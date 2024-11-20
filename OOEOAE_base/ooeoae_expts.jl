import PALEOboxes as PB

import PALEOmodel

function ooeoae_expts(model, expts)

    for expt in expts
        if length(expt) == 2 && first(expt) == "OPinit"
            _, (O_init, P_init) = expt
            PB.set_variable_attribute!(model, "atmocean.O", :initial_value, O_init)
            PB.set_variable_attribute!(model, "ocean.P", :initial_value, P_init)

        elseif length(expt) == 2 && first(expt) == "OPAinit"
            _, (O_init, P_init, A_init) = expt
            PB.set_variable_attribute!(model, "atmocean.O", :initial_value, O_init)
            PB.set_variable_attribute!(model, "ocean.P", :initial_value, P_init)
            PB.set_variable_attribute!(model, "atmocean.A", :initial_value, A_init)

        elseif length(expt) == 2 && first(expt) == "k_O2_U"
            _, (k_O2_U_min, k_O2_U_max) = expt
            PB.set_parameter_value!(model, "ocean", "oceanburial_columns", "k_O2_U_min", k_O2_U_min)
            PB.set_parameter_value!(model, "ocean", "oceanburial_columns", "k_O2_U_max", k_O2_U_max)

        elseif length(expt) == 2 && first(expt) == "k_anox"
            _, k_anox = expt
            PB.set_parameter_value!(model, "ocean", "oceanburial_columns", "k_anox", k_anox)
        
        elseif length(expt) == 2 && first(expt) == "k_O2_U_rate" # linear or other interpolation methods
            _, rate_method = expt
            PB.set_parameter_value!(model, "ocean", "oceanburial_columns", "k_O2_U_rate", rate_method)

        elseif length(expt) == 2 && first(expt) == "corg_burial_fac"
            _, corg_burial_fac = expt
            PB.set_parameter_value!(model, "ocean", "oceanburial_columns", "corg_burial_fac", corg_burial_fac)

        elseif length(expt) == 2 && first(expt) == "anoxic_threshold"
            _, anoxic_threshold = expt
            PB.set_parameter_value!(model, "ocean", "oceanburial_columns", "anoxic_threshold", anoxic_threshold)

        elseif length(expt) == 3 && first(expt) == "Ppulse"
            _, Ptimes, Pfluxes = expt
            Ptimes_all = vcat([-1e30, first(Ptimes)-1.0], Ptimes, [last(Ptimes)+1.0, 1e30])
            Pfluxes_all = vcat([0.0, 0.0], Pfluxes, [0.0, 0.0])
            PB.set_parameter_value!(model, "global", "Ppulse", "perturb_times", Ptimes_all)
            PB.set_parameter_value!(model, "global", "Ppulse", "perturb_totals", Pfluxes_all)
        
        elseif length(expt) == 3 && expt[1] == "CO2pulse"
            _, CO2times, CO2fluxes = expt
            CO2times_all = vcat([-1e30, first(CO2times)-1.0], CO2times, [last(CO2times)+1.0, 1e30])
            CO2fluxes_all = vcat([0.0, 0.0], CO2fluxes, [0.0, 0.0])
            delta_all = repeat([-5.0], length(CO2times_all))
            PB.set_parameter_value!(model, "global", "CO2pulse", "perturb_times", CO2times_all)
            PB.set_parameter_value!(model, "global", "CO2pulse", "perturb_totals", CO2fluxes_all)
            PB.set_parameter_value!(model, "global", "CO2pulse", "perturb_deltas", delta_all) 
        
        elseif length(expt)==4 && expt[1] == "CO2pulse"
            _, size, pstart, duration = expt # (_, mol C, yr, yr)
            
            CO2pulse = PB.get_reaction(model, "global", "CO2pulse")
            delta = -5       
            
            @info "expt $expt apply CO2pulse size $size (mol C) delta $delta (per mil) duration $duration (yr) start time $pstart (model yr)"
            # Witches hat perturbation
            PB.setvalue!(PB.get_parameter(CO2pulse, "perturb_times"),  
                [-1e30, pstart, pstart+duration/2, pstart+duration,  1e30]) 
            PB.setvalue!(PB.get_parameter(CO2pulse, "perturb_totals"), 
                [0.0,   0.0,    2*size/duration,    0.0,             0.0])
            PB.setvalue!(PB.get_parameter(CO2pulse, "perturb_deltas"),
                delta.*[1.0,    1.0,     1.0,       1.0,             1.0])
        
        elseif length(expt) == 4 && expt[1] == "set_initial_value"
            # generic :initial_value set (set_initial_value, <domain>, <varname>, <initial_value)
            _, domname, varname, initial_value = expt            
            PB.set_variable_attribute!(model, domname, varname, :initial_value, initial_value)
            # Do not change normalized values e.g. A_norm, otherwise unpredictable results will be involved
        elseif length(expt) == 5 && expt[1] == "set_parameter"
            # generic :initial_value set (set_initial_value, <domain>, <varname>, <initial_value)
            _, domname, reaction, varname, value = expt            
            PB.set_parameter_value!(model, domname, reaction, varname, value)
        
        elseif expt == "set_forcing_EndP"
            PB.set_variable_attribute!(model, "global", "PG", :initial_value, 0.54)
            PB.set_variable_attribute!(model, "global", "DEGASS", :initial_value, 1.13)
            PB.set_variable_attribute!(model, "global", "UPLIFT", :initial_value, 0.55)
            PB.set_variable_attribute!(model, "global", "COAL", :initial_value, 4.07)

            PB.set_variable_attribute!(model, "land", "BA_AREA", :initial_value, 0.70)
            PB.set_variable_attribute!(model, "land", "GRAN_AREA", :initial_value, 1.56)
            PB.set_variable_attribute!(model, "land", "CARB_AREA", :initial_value, 1.0)
            PB.set_variable_attribute!(model, "land", "VEG", :initial_value, 1.0)

            PB.set_variable_attribute!(model, "atmocean", "A", :initial_value, 3.193e18*1.5^0.5)
            PB.set_variable_attribute!(model, "atmocean", "O", :initial_value, 3.7e19*1.34)
        
        elseif expt == "set_restoring_EndP"    
            PB.set_parameter_value!(model, "atmocean", "restoreA", "RequiredLevel", 3.193e18*1.5^0.5)
            PB.set_parameter_value!(model, "atmocean", "restoreO2", "RequiredLevel", 3.7e19*1.34)

        elseif length(expt) == 2 && expt[1] == "P_weathering"
            _, perturb_totals = expt
            PB.set_parameter_value!(model, "land", "P_weathering", "perturb_totals", [perturb_totals, perturb_totals])

        elseif length(expt) == 5 && expt[1] == "P_weathering"
            _, k10_phosw, k_Psilw, k_Pcarbw, k_Poxidw = expt
            # remove constant P_weathering input, enable land-surface P weathering
            PB.set_parameter_value!(model, "land", "P_weathering", "perturb_totals", [0.0, 0.0]) # no constant P input flux
            PB.set_parameter_value!(model, "land", "land_fluxes", "k10_phosw", k10_phosw)
            PB.set_parameter_value!(model, "land", "land_fluxes", "k_Psilw", k_Psilw)
            PB.set_parameter_value!(model, "land", "land_fluxes", "k_Pcarbw", k_Pcarbw)
            PB.set_parameter_value!(model, "land", "land_fluxes", "k_Poxidw", k_Poxidw)
        
        elseif length(expt) == 3 && expt[1] == "COPSE_locb"
            _, CPland0, k11_landfrac = expt
            PB.set_parameter_value!(model, "land", "land_fluxes", "CPland0", CPland0)
            PB.set_parameter_value!(model, "land", "land_fluxes", "k11_landfrac", k11_landfrac)
        
        elseif length(expt) == 7 && expt[1] == "COPSE_landw"
            _, k_basw, k_granw, k14_carbw, k17_oxidw, k10_phosw, k15_plantenhance = expt
            PB.set_parameter_value!(model, "land", "land_fluxes", "k_basw", k_basw)
            PB.set_parameter_value!(model, "land", "land_fluxes", "k_granw", k_granw)
            PB.set_parameter_value!(model, "land", "land_fluxes", "k14_carbw", k14_carbw)
            PB.set_parameter_value!(model, "land", "land_fluxes", "k17_oxidw", k17_oxidw)
            PB.set_parameter_value!(model, "land", "land_fluxes", "k10_phosw", k10_phosw)
            PB.set_parameter_value!(model, "land", "land_weathering", "k15_plantenhance", k15_plantenhance)
        
        elseif length(expt) == 5 && expt[1] == "COPSE_sedcrust"
            _, C_delta, G_delta, k12_ccdeg, k13_ocdeg = expt
            PB.set_variable_attribute!(model, "sedcrust", "C", :initial_delta, C_delta)
            PB.set_variable_attribute!(model, "sedcrust", "G", :initial_delta, G_delta)
            PB.set_parameter_value!(model, "sedcrust", "sedcrust_copse", "k12_ccdeg", k12_ccdeg)
            PB.set_parameter_value!(model, "sedcrust", "sedcrust_copse", "k13_ocdeg", k13_ocdeg)

        elseif length(expt) == 2 && expt[1] == "land_flux_Bergman"
            _, bool_locb = expt # constant weathering input, constant burial fluxes
            PB.set_parameter_value!(model, "land", "P_weathering", "perturb_totals", [3.9e10, 3.9e10])
            PB.set_parameter_value!(model, "land", "Corg_oxidation", "k_oxidw", 7.75e12)
            PB.set_parameter_value!(model, "land", "const_ocdeg", "perturb_totals", [-1.25e12, -1.25e12])
            PB.set_parameter_value!(model, "ocean", "oceanburial_columns", "k_mocb", 4.5e12)
            if bool_locb
                PB.set_parameter_value!(model, "land", "const_locb", "perturb_totals", [4.5e12, 4.5e12])
            else
                PB.set_parameter_value!(model, "land", "const_locb", "perturb_totals", [0.0, 0.0])
            end
        
        elseif length(expt) == 2 && expt[1] == "land_flux_COPSE_reload"
            _, bool_locb = expt # constant weathering input, constant burial fluxes
            PB.set_parameter_value!(model, "land", "P_weathering", "perturb_totals", [4.0e10, 4.0e10])
            PB.set_parameter_value!(model, "land", "Corg_oxidation", "k_oxidw", 3.75e12)
            PB.set_parameter_value!(model, "land", "const_ocdeg", "perturb_totals", [-1.25e12, -1.25e12])
            PB.set_parameter_value!(model, "ocean", "oceanburial_columns", "k_mocb", 2.5e12)
            if bool_locb
                PB.set_parameter_value!(model, "land", "const_locb", "perturb_totals", [2.5e12, 2.5e12])
            else
                PB.set_parameter_value!(model, "land", "const_locb", "perturb_totals", [0.0, 0.0])
            end

        elseif length(expt) == 3 && expt[1] == "anoxic_area"
            _, range_min, range_max = expt # , rangevar
            PB.set_parameter_value!(model, "oceanfloor", "anoxic_area", "range_min", range_min)
            PB.set_parameter_value!(model, "oceanfloor", "anoxic_area", "range_max", range_max)
    

        elseif length(expt) == 4 && expt[1] == "CPU_U"
            _, k_U_total_silw, k_U_anoxic, k_U_other = expt
            PB.set_parameter_value!(model, "land", "U_land", "k_U_total_silw", k_U_total_silw)
            PB.set_parameter_value!(model, "oceanfloor", "U_oceanfloor", "k_U_anoxic", k_U_anoxic)
            PB.set_parameter_value!(model, "oceanfloor", "U_oceanfloor", "k_U_other", k_U_other)
        
        elseif length(expt) == 3 && expt[1] == "CPsea"
            _, k_oxic, k_anoxic = expt
            PB.set_parameter_value!(model, "ocean", "oceanburial_columns", "k_oxic", k_oxic)
            PB.set_parameter_value!(model, "ocean", "oceanburial_columns", "k_anoxic", k_anoxic)

        else
            error("unrecognized expt=", expt)
        end

    end

    return nothing
end

