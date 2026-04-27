
import PALEOboxes as PB

import PALEOmodel
import DataFrames


function SPICE_expts(
    model, expts;
    modelpars=Dict{}(),
)
    ###############################################
    # choose an 'expt' (a delta to the base model)
    ###############################################
    for expt in expts
        if expt == "baseline"
            # baseline configuration

        elseif length(expt)==5 && expt[1] == "set_interp_forcing"
            _, domian, reaction, force_times, force_values = expt #
            force = PB.get_reaction(model, domian, reaction) 
            PB.setvalue!(PB.get_parameter(force, "force_times"), force_times)
            PB.setvalue!(PB.get_parameter(force, "force_values"), force_values)

        elseif length(expt)==5 && expt[1] == "set_perturb"
            _, domian, reaction, force_times, force_values = expt #
            force = PB.get_reaction(model, domian, reaction) 
            PB.setvalue!(PB.get_parameter(force, "perturb_times"), force_times)
            PB.setvalue!(PB.get_parameter(force, "perturb_totals"), force_values)

        elseif length(expt)==2 && expt[1] == "set_O2initial"
            _, multipler = expt #
            # force = PB.get_reaction(model, "atmocean", "O") 
            # PB.set_attribute!(PB.get_parameter(force, "initial_value"), multipler*3.7e19)    
            PB.set_variable_attribute!(model, "atmocean", "O", :initial_value, multipler*3.7e19)

        elseif length(expt)==2 && expt[1] == "set_Sweathering"
            _, multipler = expt #
            force = PB.get_reaction(model, "land", "land_fluxes") 
            PB.setvalue!(PB.get_parameter(force, "k22_gypw"), 2.0e12*multipler)
            PB.setvalue!(PB.get_parameter(force, "k21_pyrw"), 0.45e12*multipler)


        else
            error("unknown expt ", expt)
        end
    end

    # set parameters for carbonate-silicate weathering steady-state
    rct_ocean       = PB.get_reaction(model, "ocean", "oceanburial_copse")
    mtotpb = rct_ocean.pars.k2_mocb[]/rct_ocean.pars.CPsea0[]  +rct_ocean.pars.k7_capb[] +rct_ocean.pars.k6_fepb[]
    rct_degass      = PB.get_reaction(model, "sedcrust", "sedcrust_copse")
    rct_sfw         = PB.get_reaction(model, "oceanfloor", "sfw")
    rct_land_fluxes = PB.get_reaction(model, "land", "land_fluxes")
    PALEOcopse.Land.LandCOPSEReloaded.set_steady_state(
        rct_land_fluxes, 
        rct_ocean.pars.k2_mocb[],
        mtotpb,
        rct_degass.pars.k12_ccdeg[],
        rct_degass.pars.k13_ocdeg[], 
        0.0,
        rct_sfw.pars.k_sfw[],
    )

    # set parameters for Sr steady state
    rct_Sr_land       = PB.get_reaction(model, "land", "Sr_land")
    rct_Sr_oceanfloor = PB.get_reaction(model, "oceanfloor", "Sr_oceanfloor")
    PALEOcopse.BioGeoChem.Strontium.set_Sr_fluxes_steady_state!(
        rct_Sr_land,
        rct_Sr_oceanfloor,
        rct_land_fluxes.pars.k_basw[],
        rct_land_fluxes.pars.k_granw[], 
        rct_land_fluxes.pars.k14_carbw[], 
        rct_sfw.pars.k_sfw[],
    )

    return nothing
end

function copse_reloaded_reloaded_plot_summary(
    outputs; 
    pager=PALEOmodel.DefaultPlotPager(),
    extrakwargs::NamedTuple=NamedTuple()
)

    pager(
        plot(title="O2", 100*PALEOmodel.get_array.(outputs, "land.mrO2"); ylabel="O_2 (%)", extrakwargs...),

        plot(title="pCO2", 1e6*PALEOmodel.get_array.(outputs, "atm.pCO2atm"); ylabel="pCO2 (ppm)", extrakwargs...),

        # ocean volume in litres = PB.Constants.k18_oceanmass/1.027 = 1.397e21 / 1.027 = 1.360e21
        # [SO4] mM = 1e3 * total SO4 (mol) / ocean vol (litres) = 7.351e-19 * SO4 (mol)
        plot(title="[SO4]", 7.351e-19*PALEOmodel.get_array.(outputs, "ocean.S"); ylabel="[SO4] mM", extrakwargs...),

        plot(title="Carbon isotopes",  outputs, "ocean.mccb_delta"; ylabel="delta 13C (‰)", extrakwargs...),
           
        plot(title="d34S SO4", outputs, "ocean.S_delta"; ylabel="d34S (‰)", extrakwargs...),
  
        plot(title="d87/d86 Sr", outputs, ["ocean.Sr_delta"]; ylabel="d87/d86 Sr", extrakwargs...),
    
        :newpage,
    )
end

function copse_reloaded_reloaded_plot(
    output; 
    pager=PALEOmodel.DefaultPlotPager(),
    extrakwargs::NamedTuple=NamedTuple()
)
    if isa(output, AbstractVector)
        have_CO2pulse = all(PB.has_variable.(output, "global.CO2pulse"))
    else
        have_CO2pulse = PB.has_variable(output, "global.CO2pulse")
    end

    # Conservation checks
    pager(
        plot(title="Total carbon",          output, ["global.total_C"], ylabel="mol C"; extrakwargs...),
        plot(title="Total carbon moldelta", output, ["global.total_C.v_moldelta"], ylabel="mol C * per mil"; extrakwargs...),
        plot(title="Total sulphur",         output, ["global.total_S"], ylabel="mol S"; extrakwargs...),
        plot(title="Total sulphur moldelta",output, ["global.total_S.v_moldelta"], ylabel="mol S * per mil"; extrakwargs...),
        plot(title="Total redox",           output, ["global.total_redox"], ylabel="mol O2 equiv"; extrakwargs...),
        :newpage,

        # Forcings
        plot(title="Physical forcings",     output, "global.".*["DEGASS", "UPLIFT", "PG"], ylim=(0, 2.0), ylabel="normalized forcing"; extrakwargs...),
        plot(title="Land area forcings",    output, ["land.BA_AREA", "land.GRAN_AREA", "global.ORGEVAP_AREA"],  ylim=(0, 2.5), ylabel="normalized forcing"; extrakwargs...),
        plot(title="Basalt area forcings",  output, ["global.CFB_area", "land.oib_area"],  ylabel="area (km^2)",),
        plot(title="Evolutionary forcings", output, "global.".*["EVO", "W", "Bforcing", "CPland_relative", "F_EPSILON", "COAL"], ylabel="normalized forcing"; extrakwargs...),    
    
        have_CO2pulse ?
            plot(title="CO2pulse", output,  ["global.CO2pulse"],  ylabel="CO2 pulse (mol C yr-1)"; extrakwargs...) : :skip,

        # Outputs
        plot(title="pCO2",                  output, "atm.".*["pCO2PAL"],  ylabel="pCO2 (PAL)"; extrakwargs...),
        plot(title="Temperature",           output, ["global.TEMP"],  ylabel="T (K)"; extrakwargs...),
        plot(title="Oxygen",                output, ["atm.pO2PAL", "ocean.ANOX"]; extrakwargs...),
        plot(title="Carbon",                output, ["global.total_C", "sedcrust.C", "sedcrust.G", "atmocean.A"], ylabel="mol C"; extrakwargs...),
        plot(title="Sulphur",               output, ["global.total_S", "ocean.S", "sedcrust.PYR", "sedcrust.GYP"], ylabel="mol S"; extrakwargs...),
        plot(title="Silicate weathering",   output, ["land.silw", "land.granw", "land.basw", "oceanfloor.sfw_total"], ylabel="flux (mol/yr)"; extrakwargs...),
        plot(title="Phosphorus weathering", output, ["land.phosw", "land.phosw_c", "land.phosw_o", "land.phosw_s"], ylabel="flux (mol P/yr)"; extrakwargs...),
        plot(title="P ocean reservoir",     output, ["ocean.P"], ylabel="P (mol)"; extrakwargs...), 
        plot(title="organic C burial",      output, ["ocean.mocb", "land.locb"], ylabel="Corg burial (mol C yr-1)"; extrakwargs...), 
        plot(title="S burial",              output, ["fluxOceanBurial.flux_total_GYP", "fluxOceanBurial.flux_total_PYR"], ylabel="S burial (mol S yr-1)"; extrakwargs...), 
        plot(title="Land biota",            output, ["land.V_T", "land.V_o2", "land.V_co2", "land.V_npp", "land.firef", "land.VEG", "global.COAL", "atm.pO2PAL"]; extrakwargs...),
        plot(title="Sulphur isotopes",      output, ["ocean.S_delta", "sedcrust.PYR_delta", "sedcrust.GYP_delta"], ylabel="delta 34S (per mil)"; extrakwargs...),
        plot(title="Carbon isotopes",       output, ["atmocean.A_delta", "ocean.DIC_delta", "atm.CO2_delta", "ocean.mccb_delta", "sedcrust.C_delta"], ylabel="delta 13C (per mil)"; extrakwargs...),
        plot(title="Sr sed reservoir",      output, ["sedcrust.Sr_sed"], ylabel="Sr_sed (mol)"; extrakwargs...),    
        plot(title="Sr ocean reservoir",    output, ["ocean.Sr"], ylabel="Sr (mol)"; extrakwargs...), 
        plot(title="Sr fluxes",             output, ["fluxRtoOcean.flux_Sr", "fluxOceanBurial.flux_total_Sr", "fluxOceanfloor.soluteflux_Sr" ], ylabel="Sr flux (mol yr-1)"; extrakwargs...),
        plot(title="Sr isotopes",           output, ["sedcrust.Sr_mantle_delta", "sedcrust.Sr_new_ig_delta", "sedcrust.Sr_old_ig_delta", "sedcrust.Sr_sed_delta", "ocean.Sr_delta"], ylabel="87Sr"; extrakwargs...),
        plot(title="U ocean reservoir",     output, ["ocean.U"], ylabel="U (mol)"; extrakwargs...), 
        plot(title="U fluxes",              output, ["fluxRtoOcean.flux_U", "oceanfloor.U_anoxic", "oceanfloor.U_other" ], ylabel="U flux (mol yr-1)"; extrakwargs...),       
        plot(title="U isotopes",            output, ["ocean.U_delta"], ylabel="d238U/235U"; extrakwargs...),

        :newpage,  # flush partial screen
    )
    return nothing
end

"""
    plot_C_Li_Sr

    Feifei's C-Li-Sr data from 355 to 345Ma
"""
function plot_C_Li_Sr(filename)
    df_C = DataFrames.DataFrame(XLSX.readtable(filename,"C"));  sort!(df_C, :Age, rev=true)
    df_Li = DataFrames.DataFrame(XLSX.readtable(filename,"Li"));  sort!(df_Li, :Age, rev=true)
    df_Sr = DataFrames.DataFrame(XLSX.readtable(filename,"Sr"));  sort!(df_Sr, :Age, rev=true)
    df_U = DataFrames.DataFrame(XLSX.readtable(filename,"U"));  sort!(df_U, :Age, rev=true)

    # to compare Jianlin's color
    # plot_DC_C = scatter(df_C.Age, df_C.d13C, color=:tan1, alpha=0.8, markersize=2.5, label=false) #, label="δ¹³Ccarb"
    # plot_DC_U = scatter(df_U.Age, df_U.d238U, color=:gray, alpha=0.1, label = "δ²³⁸Ucarb")
    # plot_DC_Li = scatter(df_Li.Age, df_Li.d7Li_add2, color=:tan1, alpha=0.8, markersize=2.5, label=false) #, label = "δ⁷Li"
    # plot_DC_Sr = scatter(df_Li.Age, df_Sr.Sr, color=:gray, alpha=0.1, label = "87Sr/86Sr")

    plot_DC_C = scatter(df_C.Age, df_C.d13C, color=:gray, alpha=0.1, label=false) #, label="δ¹³Ccarb"
    plot_DC_U = scatter(df_U.Age, df_U.d238U, color=:gray, alpha=0.1, label = "δ²³⁸Ucarb")
    plot_DC_Li = scatter(df_Li.Age, df_Li.d7Li_add2, color=:gray, alpha=0.1, label=false) #, label = "δ⁷Li"
    plot_DC_Sr = scatter(df_Li.Age, df_Sr.Sr, color=:gray, alpha=0.1, label = "87Sr/86Sr")

    return plot_DC_C, plot_DC_Li, plot_DC_Sr, plot_DC_U, unique_rows(df_C, :Age), unique_rows(df_Li, :Age), unique_rows(df_Sr, :Age), unique_rows(df_U, :Age)
end

"""
    plot_C_U_S

    Feifei's C-U data from 502 to 490 Ma SPICE
"""
function plot_C_U_S(filename; Age_multipler=-1.0)
    df_C = DataFrames.DataFrame(XLSX.readtable(filename,"C"));  sort!(df_C, :Age, rev=true)
    df_U = DataFrames.DataFrame(XLSX.readtable(filename,"U"));  sort!(df_U, :Age, rev=true)
    df_S = DataFrames.DataFrame(XLSX.readtable(filename,"S"));  sort!(df_S, :Age, rev=true)

    plot_SPICE_C = scatter(df_C.Age*Age_multipler, df_C.d13C, color=:gray, alpha=0.1, label="δ¹³Ccarb") #, label="δ¹³Ccarb"
    # plot_SPICE_U = scatter(df_U.Age*Age_multipler, df_U.d238Usw, color=:gray, alpha=0.1, label = "δ²³⁸Usw (-0.27 offset)")
    plot_SPICE_U = scatter(df_U.Age*Age_multipler, df_U.d238U, color=:gray, alpha=0.1, label = "δ²³⁸Ucarb")
    plot_SPICE_S = scatter(df_S.Age*Age_multipler, df_S.d34SCAS, color=:gray, alpha=0.1, label = "δ³⁴SCAS")

    return plot_SPICE_C, plot_SPICE_U, plot_SPICE_S, unique_rows(df_C, :Age), unique_rows(df_U, :Age), unique_rows(df_S, :Age)
end

"""
    unique_rows

Input the Dataframe, the program delete the matrix rows that duplicate the elements in the column col
"""
function unique_rows(A::DataFrames.DataFrame, col::Symbol) 
    # find the index of unique rows
    row_indices = Int[]
    row_values = []
    for i in 1:size(A, 1)
        v = A[i, col]
        if !(v in row_values)
            push!(row_indices, i)
            push!(row_values, v)
        end
    end

    return A[row_indices,:]
end
