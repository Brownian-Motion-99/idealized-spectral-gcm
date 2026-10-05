"""Reproducible optional comparison; normal Julia tests need no Fortran compiler."""
import argparse
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
CP = 287.04 / (2 / 7)


def run(argv, cwd):
    result = subprocess.run([str(arg) for arg in argv], cwd=cwd,
                            capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"Command failed: {argv}\n{result.stdout}\n{result.stderr}")
    return result.stdout


def compile_reference(isca, folder, exact=False):
    folder.mkdir(parents=True, exist_ok=True)
    kernel_path = isca / "src/shared/sat_vapor_pres/sat_vapor_pres_k.F90"
    kernel = kernel_path.read_text()
    function = re.search(r"(?is)\bfunction compute_es_k\b.*?end function compute_es_k", kernel).group()
    wrapper = """module sat_vapor_pres_mod
implicit none
contains
subroutine escomp(t,es)
real,intent(in)::t
real,intent(out)::es
real::v(1)
v=compute_es_k([t],273.16)
es=v(1)
end subroutine
subroutine descomp(t,des)
real,intent(in)::t
real,intent(out)::des
real::ep,em
call escomp(t+0.00001,ep)
call escomp(t-0.00001,em)
des=(ep-em)/0.00002
end subroutine
"""
    (folder / "saturation.f90").write_text(wrapper + function + "\nend module\n")
    source = (isca / "src/atmos_param/qe_moist_convection/qe_moist_convection.F90").read_text()
    if exact:
        # A control experiment, not the unmodified Isca result: replace only
        # the LCL lookup and conversion with a log-pressure bisection.
        pattern = (r"call get_lcl_temp\(lcl_temp_table, value, val_min, val_max, TLCL\)"
                   r"\s*pLCL = pref \* \(TLCL/theta0\)\*\*\(1\./kappa\)")
        source, changes = re.subn(pattern, "call validation_exact_lcl(theta0,r0,p_full(1), &\n"
                                "               p_full(k_surface),pLCL,TLCL)", source)
        assert changes == 1, "Unexpected Isca LCL source; review the control patch"
        helper = """
  subroutine validation_exact_lcl(theta,r0,ptop,ps,plcl,tlcl)
    real,intent(in)::theta,r0,ptop,ps
    real,intent(out)::plcl,tlcl
    real::lo,hi,mid,p,t,es
    integer::i
    lo=log(ptop)
    hi=log(ps)
    do i=1,80
      mid=0.5*(lo+hi)
      p=exp(mid)
      t=theta*(p/pref)**kappa
      call escomp(t,es)
      if (mixing_ratio(es,p)>r0) then
        hi=mid
      else
        lo=mid
      endif
    enddo
    plcl=exp(0.5*(lo+hi))
    tlcl=theta*(plcl/pref)**kappa
  end subroutine
"""
        source = source.replace("end module qe_moist_convection_mod", helper + "\nend module qe_moist_convection_mod")
    (folder / "qe_moist_convection.F90").write_text(source)
    (folder / "input.nml").write_text("&qe_moist_convection_nml\n tau_bm=7200., rhbm=.8, Tmin=173.16\n/\n")
    executable = folder / "compare_columns"
    run(["gfortran", "-fdefault-real-8", "-fdefault-double-8", "-ffree-line-length-none",
         "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", HERE / "isca_stubs.f90",
         folder / "saturation.f90", folder / "qe_moist_convection.F90",
         HERE / "isca_columns.f90", "-o", executable], folder)
    return executable


def read_columns(folder):
    lines = iter((folder / "columns.txt").read_text().splitlines())
    count, nd = map(int, next(lines).split())
    inputs = []
    for _ in range(count):
        enabled = int(next(lines))
        inputs.append((enabled, [list(map(float, next(lines).split())) for _ in range(nd)]))
    lines = iter((folder / "julia.txt").read_text().splitlines())
    julia = []
    for _ in range(count):
        header = next(lines).split()
        julia.append((header, [list(map(float, next(lines).split())) for _ in range(nd)]))
    return inputs, julia, nd


def read_reference(path, nd):
    lines = iter(path.read_text().splitlines())
    results = {}
    for line in lines:
        header = line.split()
        results[int(header[0]) - 1] = (header, [list(map(float, next(lines).split())) for _ in range(nd)])
    return results


def compare(inputs, julia, references):
    records = []
    for index, (header, values) in references.items():
        jhead, jvalues = julia[index]
        regime = "deep" if float(header[6]) > 0 else "shallow" if any(row[0] or row[1] for row in values) else "none"
        mass = [(row[2] - row[1]) / 9.8 for row in inputs[index][1]]
        energy = sum((CP * row[0] + 2.5e6 * row[1]) * m for row, m in zip(jvalues, mass))
        water = float(jhead[8]) + sum(row[1] * m for row, m in zip(jvalues, mass))
        deep_branch = None
        if regime == "deep":
            q_tolerance = 128 * len(mass) * math.ulp(1.0) * max(row[4] for row in inputs[index][1]) / 7200
            deep_branch = "humidity_timescale" if any(
                abs(row[1] - (row[3] - source[4]) / 7200) > q_tolerance
                for row, source in zip(values, inputs[index][1])) else "temperature_shift"
        records.append(dict(name=jhead[0], julia_regime=jhead[1], isca_regime=regime,
            deep_branch=deep_branch,
            lcl_match=int(jhead[2]) == int(header[2]), lzb_match=int(jhead[3]) == int(header[3]),
            temperature_error=max(abs(j[0] - r[0]) for j, r in zip(jvalues, values)),
            humidity_error=max(abs(j[1] - r[1]) for j, r in zip(jvalues, values)),
            precipitation_error=abs(float(jhead[8]) - float(header[6])),
            cape_error=abs(float(jhead[6]) - float(header[4])),
            cin_error=abs(float(jhead[7]) - float(header[5])),
            julia_energy_residual=energy, julia_water_residual=water))
    return records


def verify(inputs, julia, original, control):
    failures, exceptions = [], []
    epsilon = math.ulp(1.0)
    for index, (chead, cvalues) in control.items():
        jhead, jvalues = julia[index]
        ohead, ovalues = original[index]
        rows = inputs[index][1]
        factor = 128 * len(rows) * epsilon
        t_tol = factor * max(row[3] for row in rows) / 7200
        q_tol = factor * max(row[4] for row in rows) / 7200
        mass = [(row[2] - row[1]) / 9.8 for row in rows]
        rain_tol = q_tol * sum(mass)
        cape_tol = factor * 287.04 * max(row[3] for row in rows) * sum(math.log(row[2]/row[1]) for row in rows)
        for j, o, c in zip(jvalues, ovalues, cvalues):
            for field, tol in ((0, t_tol), (1, q_tol)):
                if abs(j[field] - c[field]) > tol:
                    failures.append(f"{jhead[0]}: exact-LCL rate disagreement")
                if abs(j[field] - o[field]) > abs(c[field] - o[field]) + tol:
                    failures.append(f"{jhead[0]}: rate error exceeds measured lookup difference")
        for jf, rf, tol in ((8, 6, rain_tol), (6, 4, cape_tol), (7, 5, cape_tol)):
            if abs(float(jhead[jf]) - float(chead[rf])) > tol:
                failures.append(f"{jhead[0]}: exact-LCL scalar disagreement")
            if abs(float(jhead[jf]) - float(ohead[rf])) > abs(float(chead[rf]) - float(ohead[rf])) + tol:
                failures.append(f"{jhead[0]}: scalar error exceeds measured lookup difference")
        # Reject spurious nonzero labels caused by floating-point residue,
        # using input-scaled operation tolerances rather than a CAPE cutoff.
        regime = "deep" if float(chead[6]) > rain_tol else "shallow" if any(
            abs(row[0]) > t_tol or abs(row[1]) > q_tol for row in cvalues) else "none"
        if regime != jhead[1]:
            failures.append(f"{jhead[0]}: meaningful exact-LCL regime disagreement")
        if int(jhead[3]) != int(chead[3]):
            failures.append(f"{jhead[0]}: exact-LCL buoyancy limit disagreement")
        lcl_difference = abs(int(jhead[2]) - int(chead[2]))
        origin_r = rows[-1][4] / (1 - rows[-1][4])
        neutral_origin = abs(float(chead[7])) <= 128 * epsilon * abs(origin_r)
        if lcl_difference:
            # A saturated input can differ by one ULP in r between equivalent
            # q<->r formulas. Both LCLs must lie in the surface cell pair.
            if not neutral_origin or max(int(jhead[2]), int(chead[2])) != len(rows) or lcl_difference > 1:
                failures.append(f"{jhead[0]}: unexplained exact-LCL LCL disagreement")
            else:
                exceptions.append(f"{jhead[0]}: surface-coincident LCL rounding")
        orig_regime = "deep" if float(ohead[6]) > rain_tol else "shallow" if any(
            abs(row[0]) > t_tol or abs(row[1]) > q_tol for row in ovalues) else "none"
        if orig_regime != regime:
            if not neutral_origin or int(jhead[2]) != len(rows) or regime != "none" or float(jhead[6]) != 0:
                failures.append(f"{jhead[0]}: unexplained lookup regime transition")
            else:
                exceptions.append(f"{jhead[0]}: lookup creates CAPE at a neutral saturated origin")
    return sorted(set(failures)), exceptions


def write_fixture(path, inputs, julia, original, control):
    lines = ["# Independent Isca outputs: unmodified SBM and an exact-LCL control.",
             "# name level pf ph_upper ph_lower T q; then dT dq rain CAPE CIN LCL LZB for each reference."]
    for index in sorted(original):
        for level, (row, o, c) in enumerate(zip(inputs[index][1], original[index][1], control[index][1]), 1):
            oh, ch = original[index][0], control[index][0]
            values = row + o[:2] + [float(oh[k]) for k in (6,4,5,2,3)] + c[:2] + [float(ch[k]) for k in (6,4,5,2,3)]
            lines.append(julia[index][0][0] + "\t" + str(level) + "\t" + "\t".join(format(x,".17g") for x in values))
    path.write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--isca-root", type=Path, required=True)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--write-fixture", type=Path)
    args = parser.parse_args()
    isca, output = args.isca_root.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    run([args.julia, f"--project={ROOT}", HERE / "export_columns.jl", output], ROOT)
    inputs, julia, nd = read_columns(output)
    report = dict(isca_revision=run(["git", "rev-parse", "HEAD"], isca).strip(),
        julia_revision=run(["git", "rev-parse", "HEAD"], ROOT).strip(),
        isca_source_sha256=hashlib.sha256((isca / "src/atmos_param/qe_moist_convection/qe_moist_convection.F90").read_bytes()).hexdigest(),
        julia_source_sha256=hashlib.sha256((ROOT / "src/Physics/Betts_Miller.jl").read_bytes()).hexdigest(),
        saturation_kernel_sha256=hashlib.sha256((isca / "src/shared/sat_vapor_pres/sat_vapor_pres_k.F90").read_bytes()).hexdigest(),
        julia_version=run([args.julia, "--version"], ROOT).strip(),
        compiler_version=run(["gfortran", "--version"], ROOT).splitlines()[0],
        constants=dict(rd=287.04, rv=461.5, cp=CP, lv=2.5e6, gravity=9.8, kappa=2/7),
        options=dict(tau=7200, rh=0.8, minimum_temperature=173.16, precision="float64",
                     saturation="direct Isca Smithsonian kernel", virtual_temperature=True),
        input_count=len(inputs), excluded=[julia[i][0][0] for i, (enabled, _) in enumerate(inputs) if not enabled])
    references = {}
    for label, exact in (("unmodified", False), ("exact_lcl_control", True)):
        build = output / label
        executable = compile_reference(isca, build, exact)
        path = build / "results.txt"
        run([executable, output / "columns.txt", path], build)
        references[label] = read_reference(path, nd)
        records = compare(inputs, julia, references[label])
        report[label] = records
        maxima = {key: max(row[key] for row in records) for key in
            ("temperature_error", "humidity_error", "precipitation_error", "cape_error", "cin_error")}
        print(label, json.dumps(dict(count=len(records), maxima=maxima,
            regime_mismatches=[r["name"] for r in records if r["julia_regime"] != r["isca_regime"]],
            index_mismatches=[r["name"] for r in records if not r["lcl_match"] or not r["lzb_match"]])))
    failures, exceptions = verify(inputs, julia, references["unmodified"], references["exact_lcl_control"])
    branches = {row["deep_branch"] for row in report["exact_lcl_control"] if row["deep_branch"]}
    if branches != {"humidity_timescale", "temperature_shift"}:
        failures.append("Reference ensemble does not cover both deep branches")
    report["validation"] = dict(passed=not failures, failures=failures, explained_boundaries=exceptions)
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print("validation", json.dumps(report["validation"]))
    if failures:
        raise SystemExit(1)
    if args.write_fixture:
        write_fixture(args.write_fixture.resolve(), inputs, julia, references["unmodified"], references["exact_lcl_control"])


if __name__ == "__main__":
    main()
