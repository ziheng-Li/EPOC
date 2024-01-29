"""
    Uranium

Uranium 238U/235U isotope system.

A 'two sink' U model configuration based on [Zhang2020](@cite) should contain the following Reservoirs and Reactions:

|Domain name         |Reservoirs     |  Reactions                           |
|:-------------------|:--------------|:-------------------------------------|
|land                |               |`ReactionULand`  | 
|ocean               |U               |                |
|oceanfloor          |               |`ReactionUOceanfloor`  | 
"""
module Uranium

import PALEOboxes as PB
using PALEOboxes.DocStrings

# import Infiltrator # Julia debugger

"""
    ReactionULand

Calculate U weathering flux from land surface, given relative (normalized) silicate (basalt + granite) rates.

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionULand{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("k_U_total_silw", 40e6, units="mol yr-1",
            description="total U weathering rate constant from silicate rocks (granite + basalt)"),
        PB.ParDouble("k_U_riv_delta", -0.26, units="per mil",
            description="d238U/235U of riverine flux"),

        PB.ParType(PB.AbstractData, "UIsotope", PB.IsotopeLinear,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable uranium isotopes and specify isotope type"),
    )
end


function PB.register_methods!(rj::ReactionULand)

    UIsotopeType = rj.pars.UIsotope[]
    PB.setfrozen!(rj.pars.UIsotope)

    vars = [
        PB.VarDepScalar("silw_relative", "",  "Basalt + granite weathering normalized to present"),
        PB.VarContribScalar("fluxRtoOcean_U"=>"fluxRtoOcean.flux_U", "mol yr-1",  "U riverine flux",
            attributes=(:field_data=>UIsotopeType,)),
    ]
 
    PB.add_method_do!(
        rj, 
        do_U_land,
        (PB.VarList_namedtuple(vars), ),
        p = UIsotopeType,
    )

    return nothing
end

function do_U_land(
    m::PB.ReactionMethod,
    pars,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    UIsotopeType = m.p

    # U weathering flux from silicate rocks
    U_silw = @PB.isotope_totaldelta(UIsotopeType, pars.k_U_total_silw[] * vars.silw_relative[], pars.k_U_riv_delta[])

    # Riverine flux
    vars.fluxRtoOcean_U[]        += U_silw

    return nothing
end


"""
    ReactionUOceanfloor

Calculate U ocean burial in 'anoxic' and 'other' sinks with partitioning linearly dependent on 'ocean anoxia'.

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionUOceanfloor{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("k_U_anoxic_D", 0.6, units="per mil",
            description="anoxic sink fractionation"),
        PB.ParDouble("k_U_other_D", 0.005, units="per mil",
            description="other sink fractionation"),
    
        PB.ParDouble("k_U_anoxic", 6.0e6, units="mol U yr-1",
            description="present-day anoxic sink flux"),
        PB.ParDouble("k_U_other", 34.0e6, units="mol U yr-1",
            description="present-day other sink flux"),
    
        PB.ParDouble("k_anox_0", 0.0025, units="",
            description="present-day ocean anoxia fraction"),

        PB.ParType(PB.AbstractData, "UIsotope", PB.IsotopeLinear,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable uranium isotopes and specify isotope type"),
    )

end


function PB.register_methods!(rj::ReactionUOceanfloor)

    UIsotopeType = rj.pars.UIsotope[]
    PB.setfrozen!(rj.pars.UIsotope)

    vars = [
        PB.VarDepScalar("ocean.ANOX", "",  "ocean anoxic fraction"),
        PB.VarDepScalar("ocean.U_norm", "",  "normalized ocean U"),
        PB.VarDep("ocean.oceanfloor.U_delta", "",  "ocean d238U"),
    
        PB.VarPropScalar("U_anoxic", "mol U yr-1",  "anoxic sink flux"),
        PB.VarPropScalar("U_other", "mol U yr-1",  "other sink flux"),
    
        PB.VarContrib("solutefluxOceanfloor_U"=>"fluxOceanfloor.soluteflux_U", "mol yr-1",  "U oceanfloor solute flux",
            attributes=(:field_data=>UIsotopeType,))   
    ]

    PB.add_method_do!(
        rj, 
        do_U_oceanfloor,
        (PB.VarList_namedtuple(vars), ),
        p = UIsotopeType,
    )

    return nothing
end

function do_U_oceanfloor(
    m::PB.ReactionMethod,
    pars,
    (vars, ), 
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction
    UIsotopeType = m.p

    vars.U_anoxic[] = pars.k_U_anoxic[]*vars.U_norm[]*(vars.ANOX[]/pars.k_anox_0[])
    vars.U_other[]  = pars.k_U_other[]*vars.U_norm[]*(1.0 - vars.ANOX[])/(1.0 - pars.k_anox_0[])
    
    r_nfloorcells = 1.0/PB.get_length(rj.domain) # fraction of flux for each oceanfloor cell
    @inbounds for i in cellrange.indices 
        U_anoxic_sink   = @PB.isotope_totaldelta(UIsotopeType, vars.U_anoxic[], vars.U_delta[i] + pars.k_U_anoxic_D[])
        U_other_sink    = @PB.isotope_totaldelta(UIsotopeType, vars.U_other[],  vars.U_delta[i] + pars.k_U_other_D[])
        vars.solutefluxOceanfloor_U[i] -= r_nfloorcells*(U_anoxic_sink + U_other_sink) 
    end

    return nothing
end

"""
    ReactionUReduction

Work-in-progress: calculate U reduction rate in sediment as a fraction of organic carbon remineralization flux `remin_Corg`:

```math
U_{rate} = remin{\\_}Corg * \\frac{oxUreducelimit}{oxUreducelimit + [O_2]} * \\frac{[U]}{k{\\_}U{\\_}conc}
```

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionUReduction{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("k_U_conc", 100e-3, units="mol m-3",
            description="U(VI) concentration to scale reaction rate"),
        PB.ParDouble("oxUreducelimit", 1e-3, units="mol m-3", 
            description="oxygen concentration below which U(VI) reduction is inhibited"),

        PB.ParDouble("k_U_D", 1.2, units="per mil",
            description="fractionation during reduction U(VI) to U(IV)"),

        PB.ParType(PB.AbstractData, "UIsotope", PB.IsotopeLinear,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable uranium isotopes and specify isotope type"),
    )
        
    stoich_reduce_U = PB.RateStoich(
        PB.VarProp("reduce_U", "mol U yr-1", "U(VI) reduction rate",
            attributes=(:calc_total=>true,)),
        ((-1.0, "U::Isotope"), (+1.0, "UIV::Isotope")),
        deltavarname_eta = ("U_delta", pars.k_U_D),  
        sms_prefix="",
        sms_suffix="_sms",
        processname="redox"
    )
    
end


function PB.register_methods!(rj::ReactionUReduction)
 
    UIsotopeType = rj.pars.UIsotope[]
    PB.setfrozen!(rj.pars.UIsotope)
    @info "register_methods! $(PB.fullname(rj)) UIsotopeType=$(UIsotopeType)"

    vars = [
        PB.VarDep("remin_Corg", "mol yr-1", "organic carbon remineralization rate"),
        PB.VarDep("O2_conc", "mol m-3", "O2 concentration"),
        PB.VarDep("U_conc", "mol m-3", "U(VI) solute concentration"),
        rj.stoich_reduce_U.ratevartemplate,
    ]

    PB.add_method_do!(rj, do_U_reduction_rate, (PB.VarList_namedtuple(vars), ) )

    PB.add_method_do!(rj, rj.stoich_reduce_U, isotope_data=UIsotopeType)

    PB.add_method_do_totals_default!(rj)

    PB.add_method_initialize_zero_vars_default!(rj) # for total Variables

    return nothing
end



function do_U_reduction_rate(
    m::PB.ReactionMethod,
    pars,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)

    @inbounds for i in cellrange.indices
        oxUreducefac = pars.oxUreducelimit[]/(pars.oxUreducelimit[] + max(vars.O2_conc[i], 0.0))
        vars.reduce_U[i] = vars.remin_Corg[i]*oxUreducefac*max(vars.U_conc[i]/pars.k_U_conc[], 0.0)     
    end
    
    return nothing
end


end
