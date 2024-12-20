module ReactionsS_dev

import PALEOboxes as PB
using PALEOboxes.DocStrings

"""
    ReactionPyriteBurial_dev

Calculate ocean pyrite burial, proportional to ocean S, per-column pyrite burial (mol S yr-1):

    pyrb = mocb * pars.k_SC_ratio[] * vars.S_norm[]

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionPyriteBurial_dev{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("k_pyrite_D", -35.0, units="per mil",
            description="pyrite burial fractionation relative to ocean SO4"),
    
        PB.ParDouble("k_SC_ratio", 1.0/4.5*29.85, units="mol/mol",
            description="Spyrite:Corg ratio at modern S (S_norm=1)"),
      
        PB.ParBool("add_oxygen", true,
            description="true to included contribution to O_sms = 2 * pyrite burial"),

        PB.ParType(PB.AbstractData, "SIsotope", PB.IsotopeLinear,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable sulphur isotopes and specify isotope type"),
    )

end


function PB.register_methods!(rj::ReactionPyriteBurial_dev)

    SIsotopeType = rj.pars.SIsotope[]
    PB.setfrozen!(rj.pars.SIsotope)

    vars = [
        # PB.VarDepScalar("ocean.ANOX", "",  "ocean anoxic fraction"),
        PB.VarDepScalar("ocean.S_norm", "mol",  "normalized ocean SO4 (modern = 1.0)"),
        PB.VarDepScalar("ocean.S_delta", "per mil",  "ocean SO4 d34S"),
        # PB.VarDep("oceanfloor.local_anoxia", "", "oceanfloor anoxia"),
        PB.VarDep("fluxOceanBurial.flux_Corg", "mol yr-1", "ocean Corg burial flux"),
           
        PB.VarProp("S_pyrite", "mol S yr-1",  "pyrite burial flux";
            attributes=(:field_data=>SIsotopeType, :calc_total=>true,)),    
    
        PB.VarContrib("solutefluxOceanfloor_S"=>"fluxOceanfloor.soluteflux_S", "mol yr-1",  "S oceanfloor solute flux";
            attributes=(:field_data=>SIsotopeType,)),
        PB.VarContribScalar("atmocean.O_sms", "mol yr-1", "atm-ocean oxygen source - sink"),
    ]

    PB.add_method_do!(
        rj, 
        do_pyriteburial,
        (PB.VarList_namedtuple(vars), ),
        p = SIsotopeType,
    )

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end

function do_pyriteburial(
    m::PB.ReactionMethod,
    pars,
    (vars, ), 
    cellrange::PB.AbstractCellRange,
    deltat
)
    # rj = m.reaction
    SIsotopeType = m.p
    
    # r_nfloorcells = 1.0/PB.get_length(rj.domain) # fraction of flux for each oceanfloor cell
    for i in cellrange.indices
        mocb = PB.get_total(vars.flux_Corg[i]) # mol Corg yr-1 in this cell

        pyrb = mocb * pars.k_SC_ratio[] * vars.S_norm[]
        vars.S_pyrite[i] = @PB.isotope_totaldelta(
            SIsotopeType, 
            pyrb,
            vars.S_delta[] + pars.k_pyrite_D[],
        )

        vars.solutefluxOceanfloor_S[i] -= vars.S_pyrite[i]

        if pars.add_oxygen[]
            vars.O_sms[] += 2*pyrb
        end
       
    end

    return nothing
end


end # module