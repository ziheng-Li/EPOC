# PALEOexamples/src/OOEOAE_base/P_O_romglb_shelf_2023118

Connecting romglb model predictions via shelf area to idealized column model

Updated from PALEOexamples/src/OOEOAE_base/P_O_romglb_shelf_20231116

- target Corg burial reduced to 4.5e12 mol yr-1, modify C:P to target 4e10 mol P yr-1
- small changes to target Corg vs O2 utilisation to match Fig 2 in paper
- calculate "true" O2 utilisation by running a model with zero nutrients (not a major difference?)
- additional oxic shelf (now 3x shelves connected to boxes 14, 15 with different exchange rates to create range in
  low-nutrient shelves)
  Uses config transportromglbshelf7.yaml
- reduce tshelfexch to 0.5 yr from 1 yr (~halves shelf area required)



Based on configuration:

PALEOexamples/src/OOEOAE_base/3_P_O_romglb_test/P_O_romglb5_shelf.jl.
experiment ("burial_shelves4_BioProdPrest_Martin_Ozaki2011", ["../transportromglbshelf4.yaml","../P_O_romglb.yaml"], "model_BioProdPrest_Martin_Ozaki2011",  [("biopumpCorg_Martin", 0.858, 100.0)]),

- calculate oxygen utilization and Corg burial for each shelf in the catalog of ~18
- use this to derive shelf areas to map to any prescribed Corg burial distribution

Requires updated model configuration:

"../transportromglbshelf5.yaml"
- as ../transportromglbshelf4.yaml, with
  ReactionOceanTransportRomanielloShelf option to create additional oceanhlatgyresurface: Domain for nutrient input
(initial version, not used in runs below)

"../transportromglbshelf7.yaml"
- 21 shelves, using updated ReactionOceanTransportRomanielloShelf that allows shelves to be connected to two ocean boxes
  so that additional oxic shelves can be created with smaller steps in nutrient level
- ReactionOceanTransportRomanielloShelf option to create additional oceanhlatgyresurface: Domain for nutrient input

"../P_O_romglb_hlatgyre.yaml"
- Additional oceanhlatgyresurface: Domain, so nutrients (for restoring) can be added only to hlat, gyre surface boxes
 (requires updated ReactionOceanTransportRomanielloShelf)

## Installation

These examples use `ReactionOceanTransportRomanielloShelf` to read Matlab data files with the
3-column box model geometry from [Romaniello2010](@cite).

The Matlab datafiles are available as a zip file from <https://github.com/PALEOtoolkit/PALEOocean.jl/releases>,
generated from the Matlab model code available as Supplementary Information to [Romaniello2010](@cite).

The examples assume the zip file has been downloaded and unpacked to subfolder `romaniello2010_transport`, the script
`download_romaniello2010_files.jl` provides a function to do this:

    include("download_romaniello2010_files.jl")

    download_romaniello2010_files()  # download and unzip

## Tuning shelf areas to prescribed linear Corg burial vs O2 utilisation function

    julia> include("P_O_romglb_baseline2_20231218.jl")

Calculates shelf areas that would be needed to match a prescribed linear function of Corg vs shelf floor fractional O2 utilisation

Output saved in Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/P_O_romglb_shelf_20231116/P_O_romglb_baseline_20231218/

    Accumulated_Sum_Corg_burial.svg - Corg vs O2 utilization, adjusted
    Shelf_Area_Summary.svg - adjusted shelf areas    

## Calculate contours of P burial in (P, O) plane for adjusted shelf areas

    julia> include("P_O_romglb_save_grid.jl")

Runs grid over P, O of ocean-only models with prescribed shelf areas, and plots contours of constant P burial
(ie to show the dP/dt=0 nullcline aka critical manifold for a P-O model)

NB: adjusted shelf areas need to be cut-and-pasted from output of "P_O_romglb_sjd_baseline_20231116.jl"

Output (netcdf data files and plots) saved to Dropbox/BACE_OOEOAE_large/DainesLiOverleaf/OOEOAE_base/P_O_romglb_shelf_20231116/P_O_romglb_save_grid_20231218


    julia> include("P_O_romglb_plot_grid.jl")

 -> Pnullcline_Summary.svg (aka contours of constant P burial) in P-O plane.
(saved to Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/P_O_romglb_shelf_20231116/P_O_romglb_baseline_20231218)

## Make comparison with GLODAP data cloud

    julia> include("P_O_romglb_vs_GLODAP_20240102.jl")

Use the adjusted area of shelves, restoring ocean.P and atm.O2 to current values.

"../transportromglbshelf8.yaml"

the same as the transportromglbshelf7.yaml, only change the restoring target of the P_norm for hlat, and upw column, as a test