# Betts–Miller reference columns

`betts_miller_virtual_column.tsv` is a fixed, constructed near-neutral column
from the physical assessment. Its expected buoyancy diagnostics were evaluated
with the unmodified Isca `src/atmos_param/qe_moist_convection/qe_moist_convection.F90`
at revision `5c6d5ba`, in an isolated double-precision Fortran harness.

The harness used `Rd=287.04`, `Rv=461.5`, `kappa=2/7`, `cp=Rd/kappa`,
`Lv=2.5e6`, and `g=9.8`, matching JGCM defaults. Saturation vapor pressure
was evaluated directly using the Smithsonian formula from Isca's
`sat_vapor_pres_k.F90`, including the ice/liquid blend, rather than its
interpolation table. The origin is saturated, so the expected result does not
depend on LCL lookup interpolation. The Isca defaults of `tau=7200 s` and
reference RH `0.8` were used; neither affects these buoyancy diagnostics.

Expected virtual-temperature CAPE is `454.23657013397872 J/kg`, LCL is level 30,
and LZB is level 9 (top-to-bottom indexing). The column has zero CAPE with
dry-temperature buoyancy. Both configurations return zero adjustment: positive
virtual CAPE alone does not meet the thermal/drying criteria.

The fixture is read verbatim by tests. It is not regenerated from the parcel
calculation under test, and normal tests do not require Isca or a Fortran compiler.

`betts_miller_deep_column.tsv` stores a second fixed assessment column and the
expected deep-convection outputs from the same unmodified Isca revision,
constants, saturation evaluation, and options. Its columns are full pressure,
upper/lower interface pressure (Pa), environmental temperature (K), specific
humidity (kg/kg), followed by Isca temperature and humidity rates, temperature
reference, and humidity reference. The harness timestep was 1 s, so its
increments are numerically equal to rates per second.

This origin is unsaturated: Isca uses a 501 by 301 point LCL lookup table, while
Julia solves the LCL by bisection. The comparison allows this measured
interpolation difference, rather than requiring identical ascent. Expected
LCL is level 28 and LZB is level 9. Julia and Isca differ by approximately
0.016 J/kg in CAPE and 0.0015 mm/day in rain for this column; the focused tests
also check temperature and humidity profiles. The two deep closure modes can
have different heating despite agreeing on precipitation.
