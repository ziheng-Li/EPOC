# Model configurations from Daines & Li (2024)

To install (Julia 1.10 required):

Change directory to OOEOAE_base

    julia> ] activate .
    (OOEOAE_base) pkg> instantiate
    (OOEOAE_base) pkg> add https://github.com/PALEOtoolkit/NLsolve.jl#project_region
    (OOEOAE_base) pkg> precompile


# Notes

This repo is a partial copy of PALEOexamples\src\OOEOAE_base folder from
PALEOdev.jl repo, branch OOE_2021_Ziheng_dev2

From:
    commit 74fa0878494557fd969518c637bec46e6e2e43d9 (HEAD -> OOE_2021_Ziheng_dev2, origin/OOE_2021_Ziheng_dev2)
    Author: Ziheng Li <zihengli@cug.edu.cn>
    Date:   Sun Jan 28 15:00:46 2024 +0800

## Code changes to run outside PALEOdev.jl repo

Remove import PALEOreactions (part of PALEOdev.jl repository)

Copy from PALEOdev/PALEOreactions to OOEOAE_base/PALEOreactions:
- Uranium.jl, AtmReservoirs.jl, Burial.jl, OceanTransportRomanielloShelf.jl

Copy from OOEOAE_base:
- OOEOAE_base/CarbBurial_dev.jl 
    from: PALEOcopse.jl\src\oceanfloor\CarbBurial.jl
        Update the mccb for each cell...
    and include in each of the P-O-A scripts 

- ReactionsOOEOAE_dev.jl
    - ReactionOceanBurialColumn: column ocean hierarchy 
    - ReactionOxWeathMinimal: simplified oxidw
    - ReactionUOceanfloor_dev: Uburial for each cell...

- SolverFunctionsOOEOAE2.jl
    Structure and functions to find the nullclines and other items in Phase plane.

- Isoline.jl
    src to support SolverFunctionsOOEOAE2.jl

- ooeoae_plots.jl
    Romaniello ocean plotting functions


## Figure cross-reference for PNAS submitted version 2024-01

Master table of experiment parameters is OOEOAE_base\OOEOAE_TableS1_sum_all_expts_20240103.xlsx

Figure graphics and README are in subfolders of Dropbox\BACE_OOEOAE\DainesLiOverleaf\Figures,
where the README documents that scripts used and the location and names of the intermediate png, svg
for the figure panels that are then composited in Inkscape.

Scripts are modified to reproduce output and save in ./figures/ instead of Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/


    Fig 1
        Figure1_model_schematic

        (shelf area panel - see SI Fig 4 below)

    Fig 2 (P-O phase plane)
        Figure2_PO_secular_stability_oscillation
            - Step 1 script
                OOEOAE_base/1_P_O_columns_test/P_O_columns_Table2_sum_20231226.jl ->
            - Step 2 script plot into dropbox
                Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/1_P_O_columns_test/P_O_columns_Table2_sum_20231226/PO_secular_stability_ocillation_summary_20231226.svg
            - Step 3 diff versions of plot deposit
                Dropbox/BACE_OOEOAE/DainesLiOverleaf/Figures/Figure2_PO_secular_stability_oscillation/PO_secular_stability_ocillation_summary_20231227.svg

    Fig 3 (P-O-A limit cycles)
        Figure4_POAU_3D_sum_20231106
            - Step 1 script
                OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table2_3D_sum_20231122_graphicstest.jl -> 
            - Step 2 script plot into dropbox
                Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table2_3D_sum_graphicstestP_O_A_U_Table2_3D_sum_graphicstest (all files in it) ->
            - Step 3 diff versions of plot deposit
                Dropbox/BACE_OOEOAE/DainesLiOverleaf/Figures/Figure4_POAU_3D_sum_20231106/Figure4_3D_phase_plane_sum_20240103.svg

    Fig 4 (P-O excitability and rate-dependent forcing)
        Figures/Figure5_POA_Ppulse
            - Step 1 script
                OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table3_1_20231227.jl ->
            - Step 2 script plot into dropbox
                Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table3_1/P_O_excitability_rate_dependent_oxic_20231227.svg ->
            - Step 3 diff versions of plot deposit
                Dropbox/BACE_OOEOAE/DainesLiOverleaf/Figures/Figure5_POA_Ppulse/P_O_excitability_rate_dependent_oxic_20231229.svg
            

    Fig 5 (P-O-A OAE excitation by CO2 pulse)
        Figures/Figure7_POAU_CO2pulse
            - Step 1 script
                OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table4_1_excitability_sum_20231202.jl ->
            - Step 2 script plot into dropbox
                Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table4_1_excitability_sum ->
            - Step 3 diff versions of plot deposit
                Dropbox/BACE_OOEOAE/DainesLiOverleaf/Figures/Figure7_POAU_CO2pulse/Excitability_CO2pulse_summary_20231221.svg

    Fig 6 (P-O-A OOE excitation by degassing reduction)
        Figures/Figure7_POAU_CO2pulse
            - Step 1 script
                OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table4_4_rate_sum_20231203.jl -> 
            - Step 2 script plot into dropbox
                Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table4_4_rate_sum/Rate_CO2pulse_summary_20231218.svg ->
            - Step 3 diff versions of plot deposit
                Dropbox/BACE_OOEOAE/DainesLiOverleaf/Figures/Figure7_POAU_CO2pulse/Rate_CO2pulse_summary_20231225.svg

    Fig 7 (P-O-A regime diagram)
        Figures/Figure_regimes/regimes_20231221.pdf
            - Step 1 script
                OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table6_contour_fix_intersection.jl ->
                (NB: takes ~30mins to run a grid of models)
            - Step 2 script plot into dropbox
                Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table6_contour_20231212/Fig6_Stability_diagram_Corg.svg ->
            - Step 3 diff versions of plot deposit
                Dropbox/BACE_OOEOAE/DainesLiOverleaf/Figures/Figure_regimes/regimes_20231221.svg
            

    Fig SI3 Romaniello global + shelves vs WOA evaluation
        Figures\FigureS2_romglb_transport\burial_shelves6_BioProdPrest_Martin_Ozaki2011_0p23_0p25_lowUpwHlatP.png

            OOEOAE_base/P_O_romglb_shelf_20231218/P_O_romglb_vs_GLODAP_20240102.jl
               -> figures/P_O_romglb_shelf_20231216/P_O_romglb_vs_GLODAP_20240102/burial_shelves6_BioProdPrest_Martin_Ozaki2011_0p23_0p25_lowUpwHlatP.png

            ZHL TODO: I had to uncomment L95-98 to produce a plot ?

    Fig SI4 Romaniello shelf area controls (expanded version of main paper Fig 1 panel C)

        Figures\Figure_ICBM_shelves\icbm_shelves_20231219.pdf

        OOEOAE_base/P_O_romglb_shelf_20231218/  
            see README.md, this requires a sequence of scripts:
            
            julia> include("P_O_romglb_baseline2_20231218.jl")  # tune and plot shelf areas to prescribed Corg burial vs [O2] distribution
                -> figures/P_O_romglb_shelf_20231216/P_O_romglb_baseline_20231218/
                        Shelf_Area_Summary.svg
                        Accumulated_Sum_Corg_burial.svg

            julia> include("P_O_romglb_save_grid.jl")  # calculate and save grid P vs O2
                -> figures/P_O_romglb_shelf_20231216/P_O_romglb_save_grid_20231218/
                (netcdf gridded output for use by plot script below)
            
            julia> include("P_O_romglb_plot_grid.jl")  # plot contours of constant P burial in P - O plane
                -> figures/P_O_romglb_shelf_20231216/Pnullcline_Summary.svg (aka contours of constant P burial)

    Fig SI 5 P-O hysteresis bifurcation
        Figures/FigureSI_PO_secular_stability_oscillation/PO_secular_stability_SI_fix_intersection_20240103.pdf
            OOEOAE_base/1_P_O_columns_test/P_O_columns_FigS1_20231202.jl
            -> figures\P_O_columns_FigS1_20231202
        ZHL TODO this script errors (at the end, seems to still produce the plots that were used) !!

    Fig SI P-O waveform vs O nullcline
        Figures/FigureSI_PO_secular_stability_oscillation/PO_secular_stability_SI_waveform_20240103.pdf
            OOEOAE_base/1_P_O_columns_test/P_O_columns_FigS1_20231202.jl
            -> figures\P_O_columns_FigS1_20231202

    Fig SI P-O excitability and rate-dependent forcing
        Figures/Figure5_POA_Ppulse/P_O_excitability_rate_dependent_oxic_SI_20231229.pdf
        (expanded version of main paper figure above, uses same script)
            OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table3_1_20231227.jl
            -> figures\P_O_A_U_Table3_1

    TODO P-O-A sensitivity studies

    Fig SI (OAE excitation by CO2 degassing ramp increase)
        Figures/Figure7_POAU_CO2pulse/Rate_CO2pulse_oxic_increase_summary_20240103.pdf
            OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table4_1_rate_sum_20231130.jl
            -> figures\P_O_A_U_Table4_1_rate_sum

## Additional scripts

    To find corg burial eff parameter corresponding to folds (to set start points etc):
        OOEOAE_base\2_P_O_A_U_columns_test\P_O_A_U_test_folds_nullcline.jl
        -> figures\P_O_A_U_Table6