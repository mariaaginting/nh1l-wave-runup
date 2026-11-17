# One-layer non-hydrostatic 2DH wave model (MATLAB)

Depth-integrated, one-layer non-hydrostatic wave model on a staggered two-dimensional grid, with wetting and drying.

Each time step consists of:

1. a shallow-water (hydrostatic) update, then
2. a non-hydrostatic pressure correction obtained from a sparse pentadiagonal system.

The code is intended for research and teaching. It is not a general-purpose coastal model.

## Requirements

- MATLAB R2016b or newer (local functions in `.m` files)
- No extra toolboxes
- Run scripts from this folder so that helper functions and data files are on the MATLAB path

## Quick start

Open MATLAB in this directory and run a case, for example:

```matlab
case31
```

Scripts plot during the run (`pause` each step). Comment out the plotting block to run faster.

| Script | Runtime (full grid, typical) |
|---|---|
| `case31.m`, `case32.m` | minutes |
| `case33.m`, `case34_conical_island.m`, `case35_inclined_bay.m` | up to several hours |

`case31.m` and `case35_inclined_bay.m` have a switch `use_nonhydrostatic` (default `true`). Set it to `false` for a hydrostatic shallow-water run.

## Test cases

| Script | Description | Output |
|---|---|---|
| `case31.m` | Standing wave in a closed flat basin vs. the linear dispersive solution | figures |
| `case32.m` | Solitary wave over a flat bottom vs. the analytical solitary wave | figures |
| `case33.m` | Periodic wave run-up on a 1:25 plane beach (Carrier–Greenspan), embedded influx | `results_carrier_greenspan.mat` |
| `case34_conical_island.m` | Solitary wave run-up around a conical island ([Briggs et al., 1995](https://doi.org/10.1029/95JC02192)). Set `VARIANT = 'A'`, `'B'` or `'C'` | `results_conical_island_<A,B,C>.mat` |
| `case35_inclined_bay.m` | Solitary wave run-up in an inclined bay (Özkan-Haller & Kirby, 1997; Hubbard & Dodd, 2002) | `results_inclined_bay.mat` |

On a moving shoreline (`case35`), the pressure solver uses wet–wet couplings only (`options.wet_wet_coupling = true`) so that dry cells do not make the Poisson system singular.

## Post-processing

After the case scripts have written `results_*.mat` files, run:

```matlab
generate_summary_tables
```

This prints:

1. error metrics against reference or analytical data (RMSE, NRMSE, correlation)
2. mass conservation and positivity
3. computational cost (wall time, factorization/solve)

Inclined-bay reference curves are not bundled; that row is skipped until a reference is supplied in `generate_summary_tables.m`.

## Repository layout

| File | Role |
|---|---|
| `SWE2D_engine_v2.m` | Hydrostatic momentum update (upwind advection, wet/dry mask) |
| `Pnonhydros1Nwd_new.m` | Non-hydrostatic pressure system (`options.wet_wet_coupling = true` on a moving shoreline) |
| `dispr_swash1L.m` | Dispersion relation and group velocity of the one-layer model (influx, case 3) |
| `compute_error_metrics.m` | RMSE, NRMSE, correlation |
| `closest.m`, `inv_label.m`, `freqspace.m` | Grid and index helpers |

`freqspace.m` in this folder shadows MATLAB’s built-in `freqspace` (Image Processing Toolbox). Keep this copy on the path when running case 3.

## Data

| File | Content |
|---|---|
| `Data_NOAAori.xlsx` (sheet `TS2A`) | Conical-island laboratory gauges, case A |
| `dataB.txt`, `dataC.txt` | Conical-island laboratory gauges, cases B and C |
| `cg_eksak.mat` | Analytical Carrier–Greenspan shoreline series (variable `A`: columns time, elevation) |

Gauge records follow the conical-island experiments of Briggs et al. (1995).

## Notes

- Cells are numbered *x*-fastest; see `inv_label.m`.
- Velocities `U` and `V` live on cell faces; elevation, depth, and pressure live at cell centres.
- Negative water depth is clipped in the run-up cases. Dry-face velocities are zeroed before the pressure step in `case35_inclined_bay.m`.
- Generated `results_*.mat` files and MATLAB figures are not required to reproduce the source; they are outputs of the case scripts.

## References

- Briggs, M. J., Synolakis, C. E., Harkins, G. S., & Green, D. R. (1995). Laboratory experiments of tsunami runup on a circular island. *Pure and Applied Geophysics*, 144, 569–593.
- Carrier, G. F., & Greenspan, H. P. (1958). Water waves of finite amplitude on a sloping beach. *Journal of Fluid Mechanics*, 4(1), 97–109.
- Hubbard, M. E., & Dodd, N. (2002). A 2D numerical model of wave run-up and overtopping. *Coastal Engineering*, 47(1), 1–26.
- Özkan-Haller, H. T., & Kirby, J. T. (1997). Nonlinear evolution of shear instabilities of the longshore current: A comparison of observations and computations. *Journal of Geophysical Research: Oceans*, 102(C8), 17771–17776.
- Stelling, G., & Zijlema, M. (2003). An accurate and efficient finite-difference algorithm for non-hydrostatic free-surface flow with application to wave propagation. *International Journal for Numerical Methods in Fluids*, 43(1), 1–23.

## License

This source is provided as-is for research and educational use, without warranty of any kind. You may copy and modify it for non-commercial scientific work. Please cite the laboratory and analytical sources listed above when you use the bundled data or test cases.
