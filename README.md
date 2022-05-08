# VoronoiRT

## 3D Radiative Transfer simulations on irregular grids

Code related to my master thesis, titled "3D Radiative Transfer on Irregular Grids". I simulate radiation in the solar atmosphere on an irregular grid using a Voronoi tesselation. 

The voronoi tesselation is calculated with the open source [`voro++`](https://github.com/chr1shr/voro) library. `C` files for computing and outputting grid statistics for the Voronoi tesselation are located in the [rt_preprocessing/](https://github.com/meudnaes/VoronoiRT/tree/master/rt_preprocessing) folder.

The radiative transfer simulation is written in [`julia`](https://julialang.org/), using its native thread parallelism to speed up computations. Code for the RT simulations are in the [src/](https://github.com/meudnaes/VoronoiRT/tree/master/src) directory. Results are written to file in an HDF5 format. Analysing results and making plots is done with [`python`](https://www.python.org/), located in the [python/](https://github.com/meudnaes/VoronoiRT/tree/master/python) folder.

To run the code, an atmospheric model is needed. The model has to be stored in an `HDF5` format. For my calculations, I used a *Bifrost* model astmosphere, provided by the Institute of Theoretical Astrophysics at the University of Oslo. 

Rasiative transer on irregular grids is performed and tested against regular grids with:
- **searchlight beam test** with the `compare_searchlight test`
- **Continuum calculations** in `compare_continuum.jl`
- **NLTE spectral synthesis** are performed in `compare_line.jl`

Data are written to file with `plot_utils.jl`, and plotted/analysed with `python`. Figures generated from the results are located in the image folder [img/](https://github.com/meudnaes/VoronoiRT/tree/master/img).

Any questions regarding the repository can be directed to: m.e.r.udnas@astro.uio.no

[![DOI](https://zenodo.org/badge/411663522.svg)](https://zenodo.org/badge/latestdoi/411663522)
