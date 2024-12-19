import PALEOboxes as PB

import PALEOmodel

function romglb_expts(model, expts)

    for expt in expts

        if length(expt) == 3 && expt[1] == "bioprod_lim"
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
