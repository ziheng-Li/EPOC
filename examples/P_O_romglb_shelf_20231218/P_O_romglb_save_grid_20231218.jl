using Logging
import DataFrames

using Plots 

import PALEOboxes as PB
import PALEOmodel
# import PALEOreactions
import PALEOcopse
import PALEOocean
import Sundials
import SciMLBase
import Interpolations 
import NCDatasets

global_logger(ConsoleLogger(stderr,Logging.Info))

# include("../ReactionsOOEOAE_dev.jl")
include("OceanTransportRomanielloShelf.jl")
include("SedimentationRate_dev.jl")
include("AtmReservoirs.jl")

include("Grid_output_P_O2.jl")
include("romglb_expts.jl")
include("romglb_plots.jl")

output_folder_name = "P_O_romglb_save_grid_20231218"
# dropbox_output_dir = "/home/sd336/Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/P_O_romglb_shelf_20231116" # will need to create this manually
dropbox_output_dir = joinpath(@__DIR__, "../../figures/P_O_romglb_shelf_20231116")
output_figures_dir = joinpath(
    dropbox_output_dir, 
    output_folder_name,
)
isdir(output_figures_dir) || mkdir(output_figures_dir)

# netcdf files are large so don't put them in shared Dropbox
netcdf_dir = joinpath(
    # "/home/sd336/Dropbox/BACE_OOEOAE_large/DainesLiOverleaf/OOEOAE_base/P_O_romglb_shelf_20231116", 
    dropbox_output_dir,
    output_folder_name,
)
isdir(netcdf_dir) || mkdir(netcdf_dir)

#####################################################
# Romglb ocean, P_O opened system burial + restoring_P + restoring O2
# No extra P_weathering or O_weathering
#####################################################

expts_table = [
    (
        "burial_shelves6_BioProdPrest_Martin_Ozaki2011_0p05_0p75", ["transportromglbshelf7.yaml","P_O_romglb_hlatgyre.yaml"], "model_BioProdPrest_Martin_Ozaki2011", 
        (0.05, 0.75), # O2_U min, max
        [
            ("biopumpCorg_Martin", true, 0.858, 100.0),
            ("tshelfexch", fill(0.5, 21)),
            # shelf_Afloor_adj calculated from first run with no area adjustment
            # 2nd iteration
            ("shelf_areas", [5.490676216388858e12, 3.226730287822806e12, 2.3645125636752827e12, 3.2948203530068857e12, 1.3761924597456716e12, 4.9589112640973816e11, 5.870241572434941e11, 8.710788191203156e11, 6.154468094969772e11, 4.6870955805937225e11, 4.041000307801425e11, 3.360071685123004e11, 3.012825399129616e11, 2.4499536213657147e11, 1.9369641585067056e11, 1.7569572020138535e11, 1.532535590623653e11, 9.173025844499614e10, 1.0e6, 1.0e6, 1.0e6]),
            ("VCI_Bergman", 0.124, 5e-3, 6e-3, 112.5*3, 112.5*3, 112.5*3, 4), # oxic P burial 4e10, Corg burial 4.5e12
        ],
    ),
    (
        "burial_shelves6_BioProdPrest_Martin_Ozaki2011_0p15_0p6", ["transportromglbshelf7.yaml","P_O_romglb_hlatgyre.yaml"], "model_BioProdPrest_Martin_Ozaki2011", 
        (0.15, 0.6), # O2_U min, max
        [
            ("biopumpCorg_Martin", true, 0.858, 100.0),
            ("tshelfexch", fill(0.5, 21)),
            # shelf_Afloor_adj calculated from first run with no area adjustment
             # 2nd iteration 
            ("shelf_areas", [1.0e6, 1.454348135011953e12, 3.6789715859647114e12, 5.1282267046866e12, 2.1221858901392666e12, 1.0989804089583306e12, 1.2925416209469277e12, 1.350977306427488e12, 9.533842075656481e11, 6.585772472951125e11, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6]),
            ("VCI_Bergman", 0.124, 5e-3, 6e-3, 112.5*3, 112.5*3, 112.5*3, 4), # oxic P burial 4e10, Corg burial 4.5e12
        ],
    ),
    (
        "burial_shelves6_BioProdPrest_Martin_Ozaki2011_0p23_0p25", ["transportromglbshelf7.yaml","P_O_romglb_hlatgyre.yaml"], "model_BioProdPrest_Martin_Ozaki2011", 
        (0.23, 0.25), # O2_U min, max
        [
            ("biopumpCorg_Martin", true, 0.858, 100.0),
            ("tshelfexch", fill(0.5, 21)),
            # 2nd iteration
            ("shelf_areas", [1.0e6, 1.0e6, 1.0e6, 1.9758955258573363e13, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6]),
            ("VCI_Bergman", 0.124, 5e-3, 6e-3, 112.5*3, 112.5*3, 112.5*3, 4), # oxic P burial 4e10, Corg burial 4.5e12
        ],
    ),
]

tspan = (0.0, 1e8) # yr


for (filenameroot, yamls, model_name, O2_U_target, vector_pars) in expts_table
# (filenameroot, yamls, model_name, O2_U_target, vector_pars) = expts_table[1]


    #####################################################
    # Create model
    #####################################################
    model = PB.create_model_from_config(
        joinpath.(@__DIR__, yamls), 
        model_name, 
        modelpars=Dict(),
    )

    romglb_expts(
        model, vector_pars
    )
    initial_state, modeldata = PALEOmodel.initialize!(model)

    #####################################################
    # call for grid_search and setup interpolation methods
    # save the grid into `output`
    #####################################################
    logfile = open(joinpath(netcdf_dir, filenameroot*".log"), "w")
    global_logger(SimpleLogger(logfile, Logging.Info))

    # # small grid for testing
    # grid_P_norm = collect(0.0:0.4:8.0); n_P_norm=length(grid_P_norm)
    # grid_pO2PAL = collect(0.0:1.1:5.5); n_O2_norm=length(grid_pO2PAL)

    # grid_P_norm = collect(0.0:0.2:8.0)      # 41
    grid_P_norm = collect(0.0:0.2:5.0)  
    # grid_pO2PAL = collect(0.0:0.22:5.06)    # 23
    # grid_pO2PAL = collect(0.0:0.11:2.0)
    grid_pO2PAL = collect(0.0:0.055:2.0) 

    # define a function to change ocean total P in model configuration
    P_total_modern=3.468251e15 # mol P for P_norm=1.0
    rct_restoreP = PB.get_reaction(model, "ocean", "restoreP")
    set_P_function = P_norm -> PB.setvalue!(rct_restoreP.pars.RequiredLevel, P_norm*P_total_modern)

    # define a function to change atmosphere O2 in model configuration
    rct_restoreO2 = PB.get_reaction(model, "atm", "restoreO2")
    set_O2_function = pO2PAL -> PB.setvalue!(rct_restoreO2.pars.RequiredLevel, pO2PAL*0.21*PB.Constants.k_moles1atm)
    # set_O2_function = pO2PAL -> PB.set_variable_attribute!(model, "atm", "O2", :initial_value, pO2PAL*0.21*PB.Constants.k_moles1atm)

    grid_output = Grid_output_P_O2(
        model,
        grid_P_norm,
        grid_pO2PAL;
        set_P_function,
        set_O2_function,
        tspan,
        filenameroot = joinpath(netcdf_dir, filenameroot),
        initial_state_use_last=false, # true to reuse last grid point as new initial state
    )

    # filenameroot = "P_O_romglb_shelf1"
    plot_NCDataset_P_O2(
        joinpath(netcdf_dir, filenameroot*"_cube.nc");
        P_burial_levels=[2e10, 4e10, 8e10], #, 12e10], # mol P yr-1
    )
    savefig(joinpath(output_figures_dir, "$(filenameroot)_vs_depth.svg"))

    global_logger(ConsoleLogger(stderr,Logging.Info))
    close(logfile)
end

