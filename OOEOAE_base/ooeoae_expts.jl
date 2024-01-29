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

        elseif length(expt) == 3 && expt[1] == "bioprod_lim"
            _, k_uPO4, k_KPO4 = expt
            PB.set_parameter_value!(model, "ocean", "bioprod", "k_uPO4", k_uPO4)
            PB.set_parameter_value!(model, "ocean", "bioprod", "k_KPO4", k_KPO4)
        
        elseif length(expt) == 4 && expt[1] == "biopumpCorg_Martin"
            _, transportfloor, martin_rovera, martin_depthmin = expt
            PB.set_parameter_value!(model, "ocean", "biopumpCorg", "input_frac", [1.0])
            PB.set_parameter_value!(model, "ocean", "biopumpCorg", "transportfloor", transportfloor)
            PB.set_parameter_value!(model, "ocean", "biopumpCorg", "martin_rovera", martin_rovera)
            PB.set_parameter_value!(model, "ocean", "biopumpCorg", "martin_depthmin", martin_depthmin)
        
        elseif length(expt) == 3 && expt[1] == "biopumpCorg_Sumexp"
            _, input_frac, sumexp_scale = expt
            PB.set_parameter_value!(model, "ocean", "biopumpCorg", "input_frac", input_frac)
            PB.set_parameter_value!(model, "ocean", "biopumpCorg", "sumexp_scale", sumexp_scale)
        
        elseif length(expt) == 2 && expt[1] == "volume_method"
            _, volume_method = expt
            PB.set_parameter_value!(model, "ocean", "transportromglb", "volume_method", volume_method)

        elseif length(expt) == 2 && expt[1] in ("ioceanshelf", "shelf_areas", "shelf_temps", "tshelfexch", "toceanexch", "ioceanshelf2", "toceanexch2")
            parname, parvalues = expt
            transportromglb = PB.get_reaction(model, "ocean", "transportromglb")
            PB.set_parameter_value!(model, "ocean", "transportromglb", parname, parvalues)

        elseif length(expt) == 7 && expt[1] == "shelf_setup"
            _, shelf_areas, shelf_temps, tshelfexch, toceanexch, ioceanshelf, tansportSv = expt
            transportromglb = PB.get_reaction(model, "ocean", "transportromglb")
            num_shelf = length(transportromglb.pars.shelf_names)
            PB.set_parameter_value!(model, "ocean", "transportromglb", "shelf_areas", shelf_areas)
            PB.set_parameter_value!(model, "ocean", "transportromglb", "shelf_temps", shelf_temps)
            PB.set_parameter_value!(model, "ocean", "transportromglb", "tshelfexch", tshelfexch)
            PB.set_parameter_value!(model, "ocean", "transportromglb", "toceanexch", toceanexch)
            PB.set_parameter_value!(model, "ocean", "transportromglb", "ioceanshelf", ioceanshelf)
            PB.set_parameter_value!(model, "ocean", "transportromglb", "tansportSv", tansportSv)
            
        elseif length(expt) == 2 && expt[1] == "sfw_distribution"
            _, ocean_nboxes = expt
            sfw_distribution_column = fill(1/ocean_nboxes, ocean_nboxes)
            sfw = PB.get_reaction(model, "oceanfloor", "sfw")
            PB.setvalue!(sfw.pars.sfw_distribution, sfw_distribution_column)

        elseif expt == "VCI"
            r_sedBEcorgP = PB.get_reaction(model, "oceanfloor", "sedBEcorgP")
            #                                                   5uM
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPO2"), [5e-3, 6e-3]) # mol m-3
            # P:Corg burial ratios to give 4e+10 mol P yr-1
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPorgCorg"), [1/5000, 1/250])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPFeCorg"),  [1/5000, 1/250])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPauthCorg"),[1/2500, 1/125])
        
        elseif length(expt) == 8 && expt[1] == "VCI_Bergman"
            _, BECorgNorm, BPO2l, BPO2h, k_PCorg_oxic, k_PFe_oxic, k_Pauth_oxic, factor_OtoAno = expt
            r_sedBEcorgP = PB.get_reaction(model, "oceanfloor", "sedBEcorgP")
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BECorgNorm"), BECorgNorm)
            #                                                   5uM
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPO2"), [BPO2l, BPO2h]) # mol m-3
            # P:Corg burial ratios to give 4e+10 mol P yr-1
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPorgCorg"), [1/(factor_OtoAno*k_PCorg_oxic), 1/k_PCorg_oxic])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPFeCorg"),  [1/(factor_OtoAno*k_PFe_oxic), 1/k_PFe_oxic])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPauthCorg"),[1/(factor_OtoAno*k_Pauth_oxic), 1/k_Pauth_oxic])
        
        elseif length(expt) == 5 && expt[1] == "VCI_fixCP"
            _, BECorgNorm, k_PCorg_oxic, k_PFe_oxic, k_Pauth_oxic = expt
            r_sedBEcorgP = PB.get_reaction(model, "oceanfloor", "sedBEcorgP")
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BECorgNorm"), BECorgNorm)
            #                                                   5uM
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPO2"), [1.0]) # mol m-3
            # P:Corg burial ratios to give 4e+10 mol P yr-1
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPorgCorg"), [1/k_PCorg_oxic])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPFeCorg"),  [1/k_PFe_oxic])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPauthCorg"),[1/k_Pauth_oxic])

        elseif expt == "VCIsmooth"
            r_sedBEcorgP = PB.get_reaction(model, "oceanfloor", "sedBEcorgP")
            #                                                   5uM
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPO2"), [5e-3, 10e-3]) # mol m-3
            # P:Corg burial ratios to give 4e+10 mol P yr-1
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPorgCorg"), [1/5000, 1/250])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPFeCorg"),  [1/5000, 1/250])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPauthCorg"),[1/2500, 1/125])

        elseif expt == "VCIx4"
            # factor-of-4 change in C:P oxic vs anoxic
            r_sedBEcorgP = PB.get_reaction(model, "oceanfloor", "sedBEcorgP")
            #                                                   5uM
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPO2"), [5e-3, 6e-3]) # mol m-3
            # P:Corg burial ratios to give 4e+10 mol P yr-1
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPorgCorg"), [1/1000, 1/250])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPFeCorg"),  [1/1000, 1/250])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPauthCorg"),[1/500, 1/125])

        elseif expt == "VCIx4_highBECorgNorm"
            # factor-of-4 change in C:P oxic vs anoxic
            r_sedBEcorgP = PB.get_reaction(model, "oceanfloor", "sedBEcorgP")
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BECorgNorm"), 0.25)
            #                                                   5uM
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPO2"), [5e-3, 6e-3]) # mol m-3
            # P:Corg burial ratios to give 4e+10 mol P yr-1
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPorgCorg"), [1/1000, 1/250])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPFeCorg"),  [1/1000, 1/250])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPauthCorg"),[1/500, 1/125])
    
        elseif expt == "VCIx4_lowBECorgNorm"
            # factor-of-4 change in C:P oxic vs anoxic
            r_sedBEcorgP = PB.get_reaction(model, "oceanfloor", "sedBEcorgP")
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BECorgNorm"), 0.05)
            #                                                   5uM
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPO2"), [5e-3, 6e-3]) # mol m-3
            # P:Corg burial ratios to give 4e+10 mol P yr-1
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPorgCorg"), [1/1000, 1/250])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPFeCorg"),  [1/1000, 1/250])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPauthCorg"),[1/500, 1/125])

        elseif expt == "VCIx4smooth"
            # factor-of-4 change in C:P oxic vs anoxic, less abrupt shift in C:P
            r_sedBEcorgP = PB.get_reaction(model, "oceanfloor", "sedBEcorgP")
            #                                                   5uM
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPO2"), [5e-3, 10e-3]) # mol m-3
            # P:Corg burial ratios to give 4e+10 mol P yr-1
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPorgCorg"), [1/1000, 1/250])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPFeCorg"),  [1/1000, 1/250])
            PB.setvalue!(PB.get_parameter(r_sedBEcorgP, "BPauthCorg"),[1/500, 1/125])

        else
            error("unrecognized expt=", expt)
        end

    end

    return nothing
end

function ooeoae_copse_reloaded_expts(
    basemodel, expts;
    modelpars=Dict{}(),
    )

    if basemodel == "reloaded_test"
        
        # modelpars = Dict("CIsotope"=>"ScalarData", "SIsotope"=>"ScalarData", "CIsotopeReacts"=>false)       
        model = PB.create_model_from_config(
            joinpath(@__DIR__, "cfg_OOEOAE_columns.yaml"), 
            "model_COPSE_test",
            modelpars=modelpars,
            # sort_methods_algorithm=PB.dfs_methods
        )

        # set parameters for carbonate-silicate weathering steady-state
        rct_ocean       = PB.get_reaction(model, "ocean", "oceanburial_copse")
        mtotpb = rct_ocean.pars.k2_mocb.v/rct_ocean.pars.CPsea0.v  +rct_ocean.pars.k7_capb.v +rct_ocean.pars.k6_fepb.v
        rct_degass      = PB.get_reaction(model, "sedcrust", "sedcrust_copse")
        rct_sfw         = PB.get_reaction(model, "oceanfloor", "sfw")
        rct_land_fluxes = PB.get_reaction(model, "land", "land_fluxes")
        PALEOcopse.Land.LandCOPSEReloaded.set_steady_state(
            rct_land_fluxes, 
            rct_ocean.pars.k2_mocb.v,
            mtotpb,
            rct_degass.pars.k12_ccdeg.v,
            rct_degass.pars.k13_ocdeg.v, 
            0.0,
            rct_sfw.pars.k_sfw.v
        )

        # set parameters for Sr steady state
        rct_Sr_land       = PB.get_reaction(model, "land", "Sr_land")
        rct_Sr_oceanfloor = PB.get_reaction(model, "oceanfloor", "Sr_oceanfloor")
        PALEOcopse.BioGeoChem.Strontium.set_Sr_fluxes_steady_state!(
            rct_Sr_land,
            rct_Sr_oceanfloor,
            rct_land_fluxes.pars.k_basw.v,
            rct_land_fluxes.pars.k_granw.v, 
            rct_land_fluxes.pars.k14_carbw.v, 
            rct_sfw.pars.k_sfw.v
        )

        run = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    elseif basemodel == "ooeoae_COPSE_column"
        model = PB.create_model_from_config(
            joinpath(@__DIR__, "cfg_OOEOAE_columns.yaml"), 
            "model_ooeoae_COPSE_column",
            modelpars=modelpars,
            # sort_methods_algorithm=PB.dfs_methods
        )   
        
        run = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())
        
    else
        error("unknown basemodel ", basemodel)
    end

    ###############################################
    # choose an 'expt' (a delta to the base model)
    ###############################################
    for expt in expts
        if expt == "baseline"
            

        # elseif length(expt)==2 && expt[1] == "tforce_constant"
        #     # Set constant forcing time
        #     # NB: use with 'modelpars'=Dict("tforcevar"=>"tforce_constant")        
        #     PB.set_variable_attribute!(model, "global", "tforce_constant", :initial_value, expt[2])

        # elseif  length(expt) == 2 && expt[1] == "CPsea"
        #     _, CPsea_index = expt
        #     PB.set_parameter_value!(model, "ocean", "oceanburial_copse", "f_CPsea", CPsea_index)
        else
            error("unknown expt ", expt)
        end
    end

    # # set parameters for carbonate-silicate weathering steady-state
    # rct_ocean       = PB.get_reaction(model, "ocean", "oceanburial_copse")
    # mtotpb = rct_ocean.pars.k2_mocb.v/rct_ocean.pars.CPsea0.v  +rct_ocean.pars.k7_capb.v +rct_ocean.pars.k6_fepb.v
    # rct_degass      = PB.get_reaction(model, "sedcrust", "sedcrust_copse")
    # rct_sfw         = PB.get_reaction(model, "oceanfloor", "sfw")
    # rct_land_fluxes = PB.get_reaction(model, "land", "land_fluxes")
    # PALEOcopse.Land.LandCOPSEReloaded.set_steady_state(
    #     rct_land_fluxes, 
    #     rct_ocean.pars.k2_mocb.v,
    #     mtotpb,
    #     rct_degass.pars.k12_ccdeg.v,
    #     rct_degass.pars.k13_ocdeg.v, 
    #     0.0,
    #     rct_sfw.pars.k_sfw.v
    # )

    # # set parameters for Sr steady state
    # rct_Sr_land       = PB.get_reaction(model, "land", "Sr_land")
    # rct_Sr_oceanfloor = PB.get_reaction(model, "oceanfloor", "Sr_oceanfloor")
    # PALEOcopse.BioGeoChem.Strontium.set_Sr_fluxes_steady_state!(
    #     rct_Sr_land,
    #     rct_Sr_oceanfloor,
    #     rct_land_fluxes.pars.k_basw.v,
    #     rct_land_fluxes.pars.k_granw.v, 
    #     rct_land_fluxes.pars.k14_carbw.v, 
    #     rct_sfw.pars.k_sfw.v
    # )

    return run
end

function copse_reloaded_reloaded_plot(
    output; 
    pager=PALEOmodel.DefaultPlotPager(),
    extrakwargs::NamedTuple=NamedTuple(),
    )
    # if isa(output, AbstractVector)
    #     have_CO2pulse = all(PB.has_variable.(output, "global.CO2pulse"))
    # else
    #     have_CO2pulse = PB.has_variable(output, "global.CO2pulse")
    # end

    # Conservation checks
     pager(
        plot(title="Total carbon",          output, ["global.total_C"], ylabel="mol C"; extrakwargs...),
        plot(title="Total carbon moldelta", output, ["global.total_C.v_moldelta"], ylabel="mol C * per mil"; extrakwargs...),
        plot(title="Total sulphur",         output, ["global.total_S"], ylabel="mol S"; extrakwargs...),
        plot(title="Total sulphur moldelta",output, ["global.total_S.v_moldelta"], ylabel="mol S * per mil"; extrakwargs...),
        plot(title="Total redox",           output, ["global.total_redox"], ylabel="mol O2 equiv"; extrakwargs...),
        # :newpage,

        # Forcings
        plot(title="Physical forcings",     output, "global.".*["DEGASS", "UPLIFT", "PG"], ylim=(0, 2.0), ylabel="normalized forcing"; extrakwargs...),
        plot(title="Land area forcings",    output, ["land.BA_AREA", "land.GRAN_AREA", "global.ORGEVAP_AREA"],  ylim=(0, 2.5), ylabel="normalized forcing"; extrakwargs...),
        plot(title="Basalt area forcings",  output, ["global.CFB_area", "land.oib_area"],  ylabel="area (km^2)",),
        plot(title="Evolutionary forcings", output, "global.".*["EVO", "W", "Bforcing", "CPland_relative", "F_EPSILON", "COAL"], ylabel="normalized forcing"; extrakwargs...),    
    
        # have_CO2pulse ?
        #     plot(title="CO2pulse", output,  ["global.CO2pulse"],  ylabel="CO2 pulse (mol C yr-1)"; extrakwargs...) : :skip,

        # # Outputs
        plot(title="pCO2",                  output, "atm.".*["pCO2PAL"],  ylabel="pCO2 (PAL)"; extrakwargs...),
        plot(title="Temperature",           output, ["global.TEMP"],  ylabel="T (K)"; extrakwargs...),
        plot(title="Oxygen",                output, ["atm.pO2PAL", "ocean.ANOX"]; extrakwargs...),
        # plot(title="Carbon",                output, ["global.total_C", "sedcrust.C", "sedcrust.G", "atmocean.A"], ylabel="mol C"; extrakwargs...),
        # plot(title="Sulphur",               output, ["global.total_S", "ocean.S", "sedcrust.PYR", "sedcrust.GYP"], ylabel="mol S"; extrakwargs...),
        # plot(title="Silicate weathering",   output, ["land.silw", "land.granw", "land.basw", "oceanfloor.sfw_total"], ylabel="flux (mol/yr)"; extrakwargs...),
        plot(title="Phosphorus weathering", output, ["land.phosw", "land.phosw_c", "land.phosw_o", "land.phosw_s"], ylabel="flux (mol P/yr)"; extrakwargs...),
        plot(title="P ocean reservoir",     output, ["ocean.P"], ylabel="P (mol)"; extrakwargs...), 
        # plot(title="organic C burial",      output, ["ocean.mocb", "land.locb"], ylabel="Corg burial (mol C yr-1)"; extrakwargs...), 
        # plot(title="S burial",              output, ["fluxOceanBurial.flux_total_GYP", "fluxOceanBurial.flux_total_PYR"], ylabel="S burial (mol S yr-1)"; extrakwargs...), 
        # plot(title="Land biota",            output, ["land.V_T", "land.V_o2", "land.V_co2", "land.V_npp", "land.firef", "land.VEG", "global.COAL", "atm.pO2PAL"]; extrakwargs...),
        # plot(title="Sulphur isotopes",      output, ["ocean.S_delta", "sedcrust.PYR_delta", "sedcrust.GYP_delta"], ylabel="delta 34S (per mil)"; extrakwargs...),
        # plot(title="Carbon isotopes",       output, ["atmocean.A_delta", "ocean.DIC_delta", "atm.CO2_delta", "ocean.mccb_delta", "sedcrust.C_delta"], ylabel="delta 13C (per mil)"; extrakwargs...),
        # plot(title="Sr sed reservoir",      output, ["sedcrust.Sr_sed"], ylabel="Sr_sed (mol)"; extrakwargs...),    
        # plot(title="Sr ocean reservoir",    output, ["ocean.Sr"], ylabel="Sr (mol)"; extrakwargs...), 
        # plot(title="Sr fluxes",             output, ["fluxRtoOcean.flux_Sr", "fluxOceanBurial.flux_total_Sr", "fluxOceanfloor.soluteflux_Sr" ], ylabel="Sr flux (mol yr-1)"; extrakwargs...),
        # plot(title="Sr isotopes",           output, ["sedcrust.Sr_mantle_delta", "sedcrust.Sr_new_ig_delta", "sedcrust.Sr_old_ig_delta", "sedcrust.Sr_sed_delta", "ocean.Sr_delta"], ylabel="87Sr"; extrakwargs...),
        # plot(title="U ocean reservoir",     output, ["ocean.U"], ylabel="U (mol)"; extrakwargs...), 
        # plot(title="U fluxes",              output, ["fluxRtoOcean.flux_U", "oceanfloor.U_anoxic", "oceanfloor.U_other" ], ylabel="U flux (mol yr-1)"; extrakwargs...),       
        # plot(title="U isotopes",            output, ["ocean.U_delta"], ylabel="d238U/235U"; extrakwargs...),

        :newpage,  # flush partial screen
    )
    return nothing
end

function ooeoae_copse_column_expts(
    basemodel, expts;
    modelpars::Dict = Dict(), 
    )

    ocean_nboxes = 100

    if basemodel == "ooeoae_COPSE_column"
        model = PB.create_model_from_config(
            joinpath(@__DIR__, "cfg_OOEOAE_columns.yaml"), 
            "model_ooeoae_COPSE_column",
            modelpars=modelpars,
            # sort_methods_algorithm=PB.dfs_methods
        )
    elseif basemodel == "ooeoae_COP_column"
        model = PB.create_model_from_config(
            joinpath(@__DIR__, "cfg_OOEOAE_columns.yaml"), 
            "model_ooeoae_COP_column",
            modelpars=modelpars,
            # sort_methods_algorithm=PB.dfs_methods
        )
    elseif basemodel == "ooeoae_AOP_column"
        model = PB.create_model_from_config(
            joinpath(@__DIR__, "cfg_OOEOAE_columns.yaml"), 
            "model_ooeoae_AOP_column",
            modelpars=modelpars,
            # sort_methods_algorithm=PB.dfs_methods
        )
    else
        error("unknown basemodel ", basemodel)
    end

    for expt in expts
        if expt == "baseline"
            # baseline configuration
            # every cell has the same proportion
            # sfw_distribution_column = fill(1/ocean_nboxes, ocean_nboxes)
            # sfw = PB.get_reaction(model, "oceanfloor", "sfw")
            # PB.setvalue!(sfw.pars.sfw_distribution, sfw_distribution_column)

        # elseif length(expt) == 2 && expt[1] == "CPsea"
        #     _, CPsea_index = expt
        #     PB.set_parameter_value!(model, "ocean", "oceanburial_copse", "f_CPsea", CPsea_index)

        elseif length(expt)==2 && expt[1] == "tforce_constant"
            # Set constant forcing time
            # NB: use with 'modelpars=Dict("tforcevar"=>"global.tforce_constant")'
            PB.set_variable_attribute!(model, "global", "tforce_constant", :initial_value, expt[2])
            
        elseif  length(expt) == 3 && expt[1] == "k_O2_U"
            _, k_O2_U_min, k_O2_U_max = expt
            PB.set_parameter_value!(model, "ocean", "oceanburial_copse", "k_O2_U_min", k_O2_U_min)
            PB.set_parameter_value!(model, "ocean", "oceanburial_copse", "k_O2_U_max", k_O2_U_max)

        elseif length(expt) ==2 && first(expt) == "corg_burial_fac"
                _, corg_burial_fac = expt
                PB.set_parameter_value!(model, "ocean", "oceanburial_copse", "corg_burial_fac", corg_burial_fac)
                
        elseif length(expt) == 4 && expt[1] == "set_initial_value"
                # generic :initial_value set (set_initial_value, <domain>, <varname>, <initial_value)
                _, domname, varname, initial_value = expt            
                PB.set_variable_attribute!(model, domname, varname, :initial_value, initial_value)

        else 
            error("unknown expt ", expt)
        end
    end

    return model
end