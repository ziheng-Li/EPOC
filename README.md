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

Copy 
    OOEOAE_base/Uranium.jl from PALEOdev.jl\PALEOreactions\src\biogeochem\Uranium.jl

    and include in each of the P-O-A scripts 

## Figure cross-reference for PNAS submitted version 2024-01

Master table of experiment parameters is OOEOAE_base\OOEOAE_TableS1_sum_all_expts_20240103.xlsx

Figure graphics and README are in subfolders of Dropbox\BACE_OOEOAE\DainesLiOverleaf\Figures,
where the README documents that scripts used and the location and names of the intermediate png, svg
for the figure panels that are then composited in Inkscape.

    Fig 1
        Figure1_model_schematic

    Fig 2 (P-O phase plane)
        Figure2_PO_secular_stability_oscillation
            OOEOAE_base/1_P_O_columns_test/P_O_columns_Table2_sum_20231226.jl
            -> figures\P_O_columns_Table2_sum_20231226

    Fig 3 (P-O-A limit cycles)
        Figure4_POAU_3D_sum_20231106
            OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table2_3D_sum_20231122_graphicstest.jl
            -> figures\P_O_A_U_Table2_3D_sum_graphicstest

    Fig 4 (P-O excitability and rate-dependent forcing)
        Figures/Figure5_POA_Ppulse
            OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table3_1_20231227.jl
            -> figures\P_O_A_U_Table3_1

    Fig 5 (P-O-A OAE excitation by CO2 pulse)
        Figures/Figure7_POAU_CO2pulse
        ZHL TODO this is not in the README, the 3D plot does not show the trajectory correctly !!!
            OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table4_1_excitability_sum_20231202.jl
            -> figures\P_O_A_U_Table4_1_excitability_sum

    Fig 6 (P-O-A OOE excitation by degassing reduction)
        Figures/Figure7_POAU_CO2pulse
            OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table4_4_rate_sum_20231203.jl
            -> figures\P_O_A_U_Table4_4_rate_sum

    Fig 7 (P-O-A regime diagram)
        Figures/Figure_regimes/regimes_20231221.pdf
            ZHL TODO this svg / pdf is not mentioned in the README
            OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_Table6_contour_fix_intersection.jl
            -> figures/P_O_A_U_Table6_contour_20231212
            (NB: takes ~30mins to run a grid of models)

    Fig SI P-O hysteresis bifurcation
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