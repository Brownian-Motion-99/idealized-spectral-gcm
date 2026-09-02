# Notation

This page defines the notation used throughout the documentation. The symbols
are chosen to be unambiguous across the [Dynamical core](@ref) and
[Physical parameterizations](@ref) pages. Code identifiers are shown separately
when they differ from the mathematical notation.

## Coordinates and indices

| Symbol | Meaning | Code representation |
|:---|:---|:---|
| $\lambda$ | longitude in radians | `mesh.λc`, index `i` |
| $\phi$ | latitude in radians | `mesh.θc`, index `j` |
| $\mu = \sin\phi$ | Gaussian latitude coordinate | `mesh.sinθ` |
| $a$ | planetary radius | `config.radius`, `atmo_data.radius` |
| $k$ | full-level index, increasing downward | third array index, `1:nd` |
| $k+\tfrac12$ | interface below full level $k$ | interface index `k + 1` |
| $p$ | pressure at a full level | `grid_p_full` |
| $p_s$ | surface pressure | `grid_ps_*` |
| $p_{k+1/2}$ | interface pressure | `grid_p_half[:, :, k + 1]` |
| $\Delta p_k$ | layer pressure thickness, $p_{k+1/2}-p_{k-1/2}$ | `grid_Δp` |
| $\sigma=p/p_s$ | nondimensional pressure used by the Held--Suarez forcing | computed locally |

The source uses `θ` in several grid names for latitude. In these docs,
$\phi$ always means latitude and $\theta$ is reserved for potential
temperature.

## State variables

| Symbol | Meaning | SI units | Code representation |
|:---|:---|:---|:---|
| $u$ | eastward wind | m s$^{-1}$ | `grid_u_*` |
| $v$ | northward wind | m s$^{-1}$ | `grid_v_*` |
| $\mathbf{v}=(u,v)$ | horizontal wind | m s$^{-1}$ | — |
| $\zeta$ | vertical component of relative vorticity | s$^{-1}$ | `spe_vor_*`, `grid_vor` |
| $f=2\Omega\sin\phi$ | Coriolis parameter | s$^{-1}$ | `atmo_data.coriolis` |
| $\eta=\zeta+f$ | absolute vorticity | s$^{-1}$ | `grid_absvor` |
| $\delta=\nabla_h\!\cdot\mathbf{v}$ | horizontal divergence | s$^{-1}$ | `spe_div_*`, `grid_div` |
| $T$ | air temperature | K | `spe_t_*`, `grid_t_*` |
| $T_v$ | virtual temperature | K | `vert_coord.virtual_temperature` |
| $\Pi=\ln p_s$ | log surface pressure | 1 | `spe_lnps_*`, `grid_lnps` |
| $q$ | water-vapor specific humidity per unit moist-air mass | kg kg$^{-1}$ | `grid_q_*` |
| $\Phi$ | geopotential | m$^2$ s$^{-2}$ | `grid_geopot_*` |
| $z=\Phi/g$ | geopotential height | m | `grid_z_*` |
| $K=\tfrac12(u^2+v^2)$ | specific kinetic energy | m$^2$ s$^{-2}$ | temporary energy fields |
| $\omega=Dp/Dt$ | pressure vertical velocity | Pa s$^{-1}$ | `grid_w_full`, NetCDF `wap` |
| $M_{k+1/2}$ | downward pressure-coordinate mass flux through an interface | Pa s$^{-1}$ | `grid_M_half` |

The suffixes `_p`, `_c`, and `_n` denote the previous, current, and next
leapfrog time levels. The prefix `spe_` denotes spherical-harmonic
coefficients and `grid_` denotes Gaussian-grid values.

## Thermodynamic and parameterization symbols

| Symbol | Meaning |
|:---|:---|
| $R_d$, $R_v$ | gas constants of dry air and water vapor |
| $c_p$ | specific heat of dry air at constant pressure |
| $\kappa=R_d/c_p$ | Poisson exponent |
| $L_v$ | latent heat of vaporization |
| $\epsilon=R_d/R_v$ | molecular-weight ratio used in moisture conversions |
| $e_s(T)$ | saturation vapor pressure |
| $r_s(T,p)$ | saturation mixing ratio per unit dry-air mass |
| $q_s(T,p)$ | saturation specific humidity per unit moist-air mass |
| $T_s(\lambda,\phi)$ | prescribed lower-boundary temperature |
| $K_E$ | boundary-layer eddy diffusivity |
| $P$ | precipitation mass flux, positive downward to the surface |
| $\tau_{BM}$ | Betts--Miller relaxation time scale |

Tendencies are written with a process label, for example
$\left.\partial T/\partial t\right|_{BM}$. A finite update over a physics
substep $\Delta t_p$ is written $\Delta T$ or $\Delta q$.

## Discrete spectral notation

| Symbol | Meaning |
|:---|:---|
| $m$ | zonal Fourier wavenumber |
| $n$ | total spherical-harmonic wavenumber |
| $M$ | configured triangular truncation, `num_fourier` |
| $Y_n^m$ | spherical harmonic basis function |
| $\widehat{\chi}_n^m$ | spectral coefficient of a field $\chi$ |
| $\Lambda_n=-n(n+1)/a^2$ | horizontal Laplacian eigenvalue |
| $\Delta t$ | base model timestep, `config.Δt` |
| $\alpha$ | semi-implicit weighting, `implicit_coef` |
| $\nu_{RA}$ | Robert--Asselin filter coefficient, `robert_coef` |

