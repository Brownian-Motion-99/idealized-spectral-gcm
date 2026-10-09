"""Build experiment-centered, latitude-specific logarithmic LW LRF artifacts.

Run with the existing climlab environment. Inputs are read only. Coefficients,
validation data, and Julia export inputs are prepared in a fresh staging folder.
The established external LRF utilities supply RRTMG differencing and metrics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import h5py
import numpy as np
from netCDF4 import Dataset


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def select_archive(directory, end_day, calibration_days, heldout_end, heldout_days, finite):
    paths = sorted(
        (p for p in directory.glob("output_t*.nc") if not p.name.endswith("_plev.nc")),
        key=lambda p: int(re.search(r"output_t(\d+)", p.name)[1]),
    )
    if not paths:
        raise ValueError(f"No native output in {directory}")
    windows = {"calibration": (end_day - calibration_days, end_day),
               "heldout": (heldout_end - heldout_days, heldout_end)}
    selections = {name: [] for name in windows}
    timestamps = {name: [] for name in windows}
    coords, units, fingerprints = None, None, []
    for path in paths:
        with Dataset(path) as ds:
            times = finite(ds["time"])
            if not times.size or np.any(np.diff(times) <= 0):
                raise ValueError(f"Empty/nonmonotonic time coordinate: {path}")
            selected = {name: np.flatnonzero((times > start) & (times <= end))
                        for name, (start, end) in windows.items()}
            if not any(indices.size for indices in selected.values()):
                continue
            current = {name: finite(ds[name]) for name in ("pfull", "phalf", "lat", "lon", "pk", "bk")}
            if coords is None:
                coords, units = current, ds["time"].units
            if ds["time"].units != units or not units.startswith("days since "):
                raise ValueError(f"Time units differ or are invalid: {path}")
            for name, values in current.items():
                np.testing.assert_array_equal(values, coords[name])
            if any(ds[name].units != "hPa" for name in ("pfull", "phalf")):
                raise ValueError(f"Expected nominal pressure in hPa: {path}")
            for name, allowed in (("ta", ("K",)), ("hus", ("1", "kg/kg", "kg kg-1"))):
                if ds[name].dimensions != ("time", "pfull", "lat", "lon") or ds[name].units not in allowed:
                    raise ValueError(f"Wrong dimensions/units for {name}: {path}")
            stat = path.stat()
            fingerprints.append({"path": str(path.resolve()), "size_bytes": stat.st_size,
                                 "mtime_ns": stat.st_mtime_ns,
                                 "file_time_start": float(times[0]), "file_time_end": float(times[-1])})
            for name, indices in selected.items():
                if indices.size:
                    selections[name].append((path, indices))
                    timestamps[name].extend(times[indices])
    for name, (start, end) in windows.items():
        times = np.asarray(timestamps[name])
        dt = np.diff(times)
        if (not dt.size or np.any(dt <= 0) or not np.allclose(dt, 0.5, atol=1e-10, rtol=0)
                or not np.isclose(times[0], start + 0.5) or not np.isclose(times[-1], end)):
            raise ValueError(f"{directory.name}: incomplete or nonuniform {name} window")
        timestamps[name] = times
    return selections, timestamps, coords, units, fingerprints


def average_window(selections, times, coords, chunk_size, q0, finite, sample_times=0):
    nd, nlat, nlon = len(coords["pfull"]), len(coords["lat"]), len(coords["lon"])
    shape = (nd, nlat, nlon)
    totals = {name: np.zeros(shape) for name in ("q", "T", "chi")}
    snapshots = []
    chosen = set(np.linspace(0, len(times)-1, sample_times, dtype=int)) if sample_times else set()
    count, q_min, q_max, bm_active = 0, np.inf, -np.inf, False
    for path, indices in selections:
        before = path.stat()
        with Dataset(path) as ds:
            for first in range(indices[0], indices[-1]+1, chunk_size):
                last = min(first+chunk_size, indices[-1]+1)
                q, temperature = finite(ds["hus"], slice(first,last)), finite(ds["ta"], slice(first,last))
                if np.any((q < 0) | (q >= 1)) or np.any(temperature <= 0):
                    raise ValueError(f"Nonphysical calibration T/q: {path}")
                if np.any(finite(ds["lrf_dta_dt"], slice(first,last)) != 0):
                    raise ValueError(f"Source baseline already has active LRF: {path}")
                bm_active |= bool(np.any(finite(ds["bm_dta_dt"], slice(first,last)) != 0))
                q_min, q_max = min(q_min,float(q.min())), max(q_max,float(q.max()))
                totals["q"] += q.sum(axis=0)
                totals["T"] += temperature.sum(axis=0)
                totals["chi"] += np.log(q+q0).sum(axis=0)
                for offset in range(last-first):
                    if count+offset in chosen:
                        snapshots.append((float(times[count+offset]), temperature[offset].copy(), q[offset].copy()))
                count += last-first
        after = path.stat()
        if (before.st_size,before.st_mtime_ns) != (after.st_size,after.st_mtime_ns):
            raise ValueError(f"Archive changed during reading: {path}")
        print(f"  {path.parent.name}: averaged {path.name}, {len(indices)} records", flush=True)
    if count != len(times) or not bm_active:
        raise ValueError("Sample count mismatch or no active BM in baseline")
    return {name: total/count for name,total in totals.items()}, snapshots, q_min, q_max


def validate_samples(samples, coords, absorbers, B, q0, alpha, profile_metrics):
    import climlab
    from climlab.domain import Axis
    from validate_fixed_lrf import CASES
    state = climlab.column_state(lev=Axis(axis_type="lev",points=coords["pfull"],bounds=coords["phalf"]))
    radiation = climlab.radiation.RRTMG_LW(state=state,specific_humidity=samples[0]["q"].copy(),
                                         absorber_vmr=absorbers,icld=0)
    def heating(q):
        radiation.specific_humidity = q.copy()
        radiation.compute_diagnostics()
        values = np.asarray(radiation.TdotLW).reshape(-1).copy()
        if not np.isfinite(values).all():
            raise ValueError("Nonfinite RRTMG validation heating")
        return values
    predictions, truths = [], []
    p = coords["pfull"]
    for index,sample in enumerate(samples):
        state["Tatm"][:] = sample["T"]
        state["Ts"][:] = sample["Ts"]
        q = sample["q"]
        base = heating(q)
        lower, free = 0.1*q*(p>=700), 0.1*q*((p>=200)&(p<700))
        backgrounds = [q,q,q,q,sample["q_mean"]]
        targets = [q+lower,q-lower,q+free,q-free,q]
        truth = [heating(target)-base for target in targets[:4]]
        truth.append(base-heating(sample["q_mean"]))
        latitude_index = int(np.argmin(np.abs(coords["lat"]-sample["lat"])))
        anomalies = np.log(np.asarray(targets)+q0)-np.log(np.asarray(backgrounds)+q0)
        prediction = alpha*anomalies@B[:,:,latitude_index].T
        predictions.append(prediction); truths.append(truth)
        if (index+1)%256==0 or index+1==len(samples):
            print(f"  held-out RRTMG pairs: {index+1}/{len(samples)}", flush=True)
    predictions, truths = np.asarray(predictions),np.asarray(truths)
    metrics = profile_metrics(predictions,truths,np.diff(coords["phalf"])*100,
                              climlab.constants.cp,9.8)
    return CASES,predictions,truths,metrics


def calibration_residual(selections, coords, chunk_size, B, chi_reference, q0, finite):
    # Independently reread humidity and accumulate evaluated heating, instead of
    # asserting the tautology B @ (chi_reference - chi_reference) == 0.
    total = np.zeros_like(chi_reference)
    count = 0
    for path,indices in selections:
        with Dataset(path) as ds:
            for first in range(indices[0],indices[-1]+1,chunk_size):
                last = min(first+chunk_size,indices[-1]+1)
                q = finite(ds["hus"],slice(first,last))
                anomalies = np.log(q+q0).sum(axis=(0,3))
                samples = q.shape[0]*q.shape[3]
                anomalies -= samples*chi_reference
                total += np.einsum("ijx,jx->ix",B,anomalies,optimize=True)
                count += samples
    residual = total/count
    if np.max(np.abs(residual)) >= 1e-10:
        raise ValueError("Nonzero calibration-mean heating")
    return residual


def write_export_inputs(directory, experiment, arrays, info, hdf_path):
    binary = directory/"export_inputs"/experiment
    binary.mkdir(parents=True)
    level = np.arange(1,arrays["chi_reference"].shape[0]+1)[:,None]
    latitude = np.arange(1,arrays["chi_reference"].shape[1]+1)[None,:]
    delta = 0.02*np.sin(level*latitude)
    q = np.exp(arrays["chi_reference"]+delta)-info["q0_kg_kg"]
    if np.any(q<0):
        raise ValueError("Negative cross-language validation humidity")
    rate = info["alpha"]*np.einsum("ijx,jx->ix",arrays["B_by_lat"],delta)/86400
    export = dict(arrays,validation_q=q.T[None,:,:],validation_tendency_K_s=rate.T[None,:,:])
    for name,values in export.items():
        np.asarray(values,dtype="<f8").ravel(order="F").tofile(binary/f"{name}.f64")
    metadata = {"experiment":experiment,"nlevel":arrays["B_by_lat"].shape[0],
                "nlatitude":arrays["B_by_lat"].shape[2],"q0_kg_kg":info["q0_kg_kg"],
                "alpha":info["alpha"],"source_sha256":sha256(hdf_path),
                "metadata_json":json.dumps(info,sort_keys=True)}
    (binary/"metadata.toml").write_text("\n".join(f"{k} = {json.dumps(v)}" for k,v in metadata.items())+"\n")


def build(experiment, amplitude, args, utilities):
    finite, latitude_weights, put, calculate_latitude_kernels, sample_columns, profile_metrics, mean_statistics = utilities
    baseline = args.data_root/experiment
    selections,times,coords,units,fingerprints = select_archive(
        baseline,args.end_day,args.calibration_days,args.heldout_end,args.heldout_days,finite)
    if (len(coords["pfull"]),len(coords["lat"]),len(coords["lon"])) != (20,64,128):
        raise ValueError("Expected T42L20 archive")
    print(f"\nBuilding {experiment}: calibration {times['calibration'][0]}–{times['calibration'][-1]}",flush=True)
    calibration,_,q_min,q_max = average_window(selections["calibration"],times["calibration"],coords,
        args.chunk_size,args.q0,finite)
    q_star = calibration["q"].mean(axis=2)
    temperature = calibration["T"].mean(axis=2)
    chi_reference = calibration["chi"].mean(axis=2)
    with h5py.File(args.absorber_source) as old:
        for name in ("pfull","phalf","lat"):
            np.testing.assert_array_equal(old[name][:],coords[name])
        absorbers = {name:value[()] for name,value in old["absorber_vmr"].items()}
    intermediate = args.output/f"{experiment}.calibration.h5"
    with h5py.File(intermediate,"w") as source:
        for name in ("pfull","phalf","lat"): source[name]=coords[name]
        source["q_star"]=q_star; source["T_zonal_mean"]=temperature
        gases=source.create_group("absorber_vmr")
        for name,value in absorbers.items(): gases[name]=value
    with h5py.File(intermediate) as source:
        kernels,reference_heating,relative,fractions,sst = calculate_latitude_kernels(source,args.fraction)
    worst=float(relative.max())
    if worst>args.tolerance:
        raise ValueError(f"{experiment}: step sensitivity {worst:.3%} exceeds {args.tolerance:.3%}")
    B=kernels*(q_star+args.q0)[None,:,:]
    residual=calibration_residual(selections["calibration"],coords,args.chunk_size,B,chi_reference,args.q0,finite)
    held, snapshots, held_min, held_max = average_window(selections["heldout"],times["heldout"],coords,
        args.chunk_size,args.q0,finite,sample_times=4)
    dp=np.diff(coords["phalf"])*100
    samples=sample_columns(snapshots,held["q"],coords["lat"],coords["lon"],dp,amplitude,
                           latitude_stratified=True)
    cases,predictions,truths,metrics=validate_samples(samples,coords,absorbers,B,args.q0,args.alpha,profile_metrics)
    mean_heating=args.alpha*np.einsum("ijx,jxy->ixy",B,held["chi"]-chi_reference[:,:,None],optimize=True)
    mean_stats=mean_statistics(mean_heating,dp,latitude_weights(coords["lat"]))
    statistics={case:{"median_relative_error":float(np.nanmedian(metrics["relative_rmse"][:,index])),
                      "p90_relative_error":float(np.nanpercentile(metrics["relative_rmse"][:,index],90)),
                      "median_profile_cosine":float(np.nanmedian(metrics["profile_cosine"][:,index]))}
                for index,case in enumerate(cases)}
    info={"experiment":experiment,"baseline":str(baseline.resolve()),
          "scheme":"regularized_log_latitude_v1","q0_kg_kg":args.q0,"alpha":args.alpha,
          "calibration_time_start":float(times["calibration"][0]),"calibration_time_end":float(times["calibration"][-1]),
          "calibration_sample_count":len(times["calibration"]),"time_units":units,
          "heldout_time_start":float(times["heldout"][0]),"heldout_time_end":float(times["heldout"][-1]),
          "heldout_sample_count":len(times["heldout"]),"heldout_response_columns":len(samples),
          "selected_fraction":args.fraction,"convergence_tolerance":args.tolerance,
          "maximum_step_sensitivity":worst,"convergence_passed":True,
          "calibration_max_abs_mean_heating_K_day":float(np.abs(residual).max()),
          "calibration_humidity_range_kg_kg":[q_min,q_max],"heldout_humidity_range_kg_kg":[held_min,held_max],
          "heldout_mean_heating_rms_K_day":float(mean_stats[0]),
          "heldout_global_column_heating_W_m2":float(mean_stats[1]),
          "response_statistics":statistics,"radiation":"climlab RRTMG_LW clear sky; fixed non-water absorbers",
          "absorber_source":str(args.absorber_source.resolve()),"absorber_source_sha256":sha256(args.absorber_source),
          "sst_amplitude_K":amplitude,"calibration_sst":"271 + 29 exp(-0.5 (latitude/26)^2) K, zonal mean",
          "validation_sst":"zonal mean + amplitude exp(-0.5 (latitude/15)^2) sin(longitude)",
          "pressure_assumption":"native nominal pfull/phalf, neglecting surface-pressure variations",
          "archived_field_convention":"12-hour archived mean T/q; chi reference averages log(q_archived + q0)",
          "reference_policy":"each experiment centered separately; no shared control reference",
          "vertical_order":"top_to_bottom","matrix_dimensions":"response_level,perturbation_level,latitude",
          "source_files":fingerprints,"builder_sha256":sha256(__file__)}
    hdf_path=args.output/f"{experiment}.h5"
    with h5py.File(hdf_path,"w") as output:
        output.attrs.update(scheme=info["scheme"],complete=False,convergence_passed=True,
                            q0_kg_kg=args.q0,alpha=args.alpha,experiment=experiment,
                            baseline=info["baseline"],metadata_json=json.dumps(info,sort_keys=True))
        for name in ("pfull","phalf","lat","lon","pk","bk"): output[name]=coords[name]
        for name,values in (("B_by_lat",B),("K_by_lat",kernels),("chi_reference",chi_reference),
                            ("q_star",q_star),("T_zonal_mean",temperature),("surface_temperature",sst),
                            ("reference_heating_LW",reference_heating),("relative_kernel_change",relative),
                            ("perturbation_fractions",fractions),("calibration_mean_residual",residual)):
            output[name]=values
        gases=output.create_group("absorber_vmr")
        for name,value in absorbers.items(): gases[name]=value
        response=output.create_group("heldout_response")
        response.create_dataset("case_names",data=cases,dtype=h5py.string_dtype())
        response["prediction"]=predictions; response["truth"]=truths
        for name,values in metrics.items(): response[name]=values
        sampled=response.create_group("samples")
        for name in ("time","lat","lon","Ts","T","q","q_mean"):
            sampled[name]=np.asarray([sample[name] for sample in samples])
        output["heldout_mean_heating"]=mean_heating
        output.attrs["complete"]=True
    info["hdf5_sha256"]=sha256(hdf_path)
    (args.output/f"{experiment}.json").write_text(json.dumps(info,indent=2)+"\n")
    arrays={"B_by_lat":B,"chi_reference":chi_reference,"latitude":coords["lat"],
            "pfull":coords["pfull"],"phalf":coords["phalf"]}
    write_export_inputs(args.output,experiment,arrays,info,hdf_path)
    print(f"Completed {experiment}: step sensitivity {worst:.3%}, calibration residual "
          f"{info['calibration_max_abs_mean_heating_K_day']:.3e} K/day; "
          f"observed-anomaly median error {statistics['observed_anomaly']['median_relative_error']:.1%}",flush=True)
    return info


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root",type=Path,default=Path("/data92/garywu/undergrad_proposal"))
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--utilities",type=Path,default=Path("/home/garywu/undergrad_proposal/LRF"))
    parser.add_argument("--absorber-source",type=Path,default=Path("/home/garywu/undergrad_proposal/LRF/data/fixed_lrf.h5"))
    parser.add_argument("--calibration-days",type=float,default=1000.)
    parser.add_argument("--end-day",type=float,default=3650.)
    parser.add_argument("--heldout-days",type=float,default=500.)
    parser.add_argument("--heldout-end",type=float,default=2500.)
    parser.add_argument("--chunk-size",type=int,default=10)
    parser.add_argument("--q0",type=float,default=1e-8)
    parser.add_argument("--alpha",type=float,default=1.)
    parser.add_argument("--fraction",type=float,default=.1)
    parser.add_argument("--tolerance",type=float,default=.015)
    args=parser.parse_args()
    if args.output.exists(): parser.error("Output staging folder must be new")
    if (args.heldout_end > args.end_day-args.calibration_days or args.q0<=0 or args.alpha<0
            or not 0<args.fraction<.5 or args.chunk_size<=0):
        parser.error("Invalid or overlapping calibration/validation settings")
    sys.path.insert(0,str(args.utilities))
    from build_fixed_lrf import finite,latitude_weights,put
    from build_latitude_lrf import calculate_latitude_kernels
    from validate_fixed_lrf import sample_columns,profile_metrics
    from plot_log_mean_adjustment import mean_statistics
    args.output.mkdir(parents=True)
    utilities=(finite,latitude_weights,put,calculate_latitude_kernels,sample_columns,profile_metrics,mean_statistics)
    results=[build(name,amplitude,args,utilities)
             for name,amplitude in (("ctrl_BM",0.),("sst1.0_BM",1.),("sst2.5_BM",2.5))]
    (args.output/"build_manifest.json").write_text(json.dumps(results,indent=2)+"\n")


if __name__=="__main__": main()
