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
# include("Grid_output_P_O2.jl")
include("../GLODAPv2.2020/base_GLODAP.jl")
include("OceanTransportRomanielloShelf.jl")
include("SedimentationRate_dev.jl")
include("AtmReservoirs.jl")
include("romglb_expts.jl")
include("romglb_plots.jl")


# Archive figures location
# dropbox_output_dir = joinpath(@__DIR__, "../../figures/P_O_romglb_shelf_20231116")

# Local figures
isdir("figures") || mkdir("figures")
dropbox_output_dir = "figures"

output_folder_name = "P_O_romglb_vs_GLODAP_20240102"
output_figures_dir = joinpath(dropbox_output_dir, output_folder_name)
isdir(output_figures_dir) || mkdir(output_figures_dir)


#####################################################
# Romglb ocean, P_O opened system burial + restoring_P + restoring O2
# No extra P_weathering or O_weathering
#####################################################

# the shelf_areas are all after adjustment
expts_table = [
    # (
    #     "burial_shelves6_BioProdPrest_Martin_Ozaki2011_0p05_0p75", ["transportromglbshelf7.yaml","P_O_romglb_hlatgyre.yaml"], "model_BioProdPrest_Martin_Ozaki2011", 
    #     (0.05, 0.75), # O2_U min, max
    #     [
    #         ("biopumpCorg_Martin", true, 0.858, 100.0),
    #         ("tshelfexch", fill(0.5, 21)),
    #         # shelf_Afloor_adj calculated from first run with no area adjustment
    #         # 2nd iteration
    #         ("shelf_areas", [5.490676216388858e12, 3.226730287822806e12, 2.3645125636752827e12, 3.2948203530068857e12, 1.3761924597456716e12, 4.9589112640973816e11, 5.870241572434941e11, 8.710788191203156e11, 6.154468094969772e11, 4.6870955805937225e11, 4.041000307801425e11, 3.360071685123004e11, 3.012825399129616e11, 2.4499536213657147e11, 1.9369641585067056e11, 1.7569572020138535e11, 1.532535590623653e11, 9.173025844499614e10, 1.0e6, 1.0e6, 1.0e6]),
    #         ("VCI_Bergman", 0.124, 5e-3, 6e-3, 112.5*3, 112.5*3, 112.5*3, 4), # oxic P burial 4e10, Corg burial 4.5e12
    #     ],
    # ),
    # (
    #     "burial_shelves6_BioProdPrest_Martin_Ozaki2011_0p15_0p6", ["transportromglbshelf7.yaml","P_O_romglb_hlatgyre.yaml"], "model_BioProdPrest_Martin_Ozaki2011", 
    #     (0.15, 0.6), # O2_U min, max
    #     [
    #         ("biopumpCorg_Martin", true, 0.858, 100.0),
    #         ("tshelfexch", fill(0.5, 21)),
    #         # shelf_Afloor_adj calculated from first run with no area adjustment
    #          # 2nd iteration 
    #         ("shelf_areas", [1.0e6, 1.454348135011953e12, 3.6789715859647114e12, 5.1282267046866e12, 2.1221858901392666e12, 1.0989804089583306e12, 1.2925416209469277e12, 1.350977306427488e12, 9.533842075656481e11, 6.585772472951125e11, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6]),
    #         ("VCI_Bergman", 0.124, 5e-3, 6e-3, 112.5*3, 112.5*3, 112.5*3, 4), # oxic P burial 4e10, Corg burial 4.5e12
    #     ],
    # ),
    # (
    #     "burial_shelves6_BioProdPrest_Martin_Ozaki2011_0p23_0p25", ["transportromglbshelf7.yaml","P_O_romglb_hlatgyre.yaml"], "model_BioProdPrest_Martin_Ozaki2011", 
    #     (0.23, 0.25), # O2_U min, max
    #     [
    #         ("biopumpCorg_Martin", true, 0.858, 100.0),
    #         ("tshelfexch", fill(0.5, 21)),
    #         # 2nd iteration
    #         ("shelf_areas", [1.0e6, 1.0e6, 1.0e6, 1.9758955258573363e13, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6]),
    #         ("VCI_Bergman", 0.124, 5e-3, 6e-3, 112.5*3, 112.5*3, 112.5*3, 4), # oxic P burial 4e10, Corg burial 4.5e12
    #     ],
    # ),
    (
        "burial_shelves6_BioProdPrest_Martin_Ozaki2011_0p23_0p25_lowUpwHlatP", ["transportromglbshelf8.yaml","P_O_romglb_hlatgyre.yaml"], "model_BioProdPrest_Martin_Ozaki2011", 
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

P_O_romglb_Table5 = Dict() # all results, indexed by fileroot

df_hlat = Data_filter(["Atlantic_Ocean" "Indian_Ocean" "Pacific_Ocean"], 6000.0, "highlat", "oxygen")
df_gyre = Data_filter(["Atlantic_Ocean" "Indian_Ocean" "Pacific_Ocean"], 6000.0, "lowlat", "oxygen")
df_upw = Data_filter(["Arctic_Ocean" "Atlantic_Ocean" "Indian_Ocean" "Pacific_Ocean"], 6000.0, "all", "oxygen", 10)
df_shelf = Data_filter(["Arctic_Ocean" "Atlantic_Ocean" "Indian_Ocean" "Pacific_Ocean"], 200.0, "all", "oxygen")


for (filenameroot, yamls, model_name, O2_U_target, vector_pars) in expts_table
# (filenameroot, yamls, model_name, O2_U_target, vector_pars) = expts_table[1]


   #####################################################
    # Create model
    #####################################################
    local model = PB.create_model_from_config(
        joinpath.(@__DIR__, yamls), 
        model_name, 
        modelpars=Dict(), # default to restoring ocean.P and atm.O2
    )

    romglb_expts(
        model, vector_pars
    )
    local initial_state, modeldata = PALEOmodel.initialize!(model)

    local paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    P_O_romglb_Table5[filenameroot] = paleorun

    println("integrate, no jacobian")
    @time PALEOmodel.ODE.integrate(
        paleorun, initial_state, modeldata, tspan, 
        solvekwargs=( # https://diffeq.sciml.ai/stable/basics/common_solver_opts/
            reltol=1e-4,
        )
    )

    plot_compare_GLODAP(paleorun.output, df_hlat, df_gyre, df_upw, df_shelf, joinpath(output_figures_dir, filenameroot*".png"))
end

# for (i, exptroot) in ["burial_shelves6_BioProdPrest_Martin_Ozaki2011_0p05_0p75"]
    paleorun = P_O_romglb_Table5["burial_shelves6_BioProdPrest_Martin_Ozaki2011_0p23_0p25_lowUpwHlatP"]

    # gr(size=(650, 400))
    p1 = plot_vs_depth(paleorun.output, "ocean.Abox",                   "Ocean area (m2)", "Depth (m)", NaN, "", true)
    p2 = plot_vs_depth(paleorun.output, "ocean.bioprod/Prod_Corg",      "Bioprod (mol/yr)", "Depth (m)", NaN, "", true)
    p3 = plot_vs_depth(paleorun.output, "ocean.remin_Corg",             "Remin Corg (mol/yr)", "Depth (m)", NaN, "", true)
    # local p3 = plot_vs_depth(paleorun.output, "ocean.export_Corg", "Export Corg (mol/yr)", "Depth (m)", NaN, "", true) # export = bioprod
    p4 = plot_vs_depth(paleorun.output, "fluxOceanBurial.flux_Corg",    "Corg burial (mol/yr)", "Depth (m)", NaN, "", true)
    l = @layout[a b; c d]
    p = plot(p1, p2, p3, p4, layout = l, size=(500,500))
    display(p)
    savefig(joinpath(output_figures_dir, "TableS3_burial_shelves6_BioProdPrest_Martin_Ozaki2011_0p23_0p25_lowUpwHlatP.png"))

# data = PB.get_data(paleorun.output, "ocean.Abox")[end]