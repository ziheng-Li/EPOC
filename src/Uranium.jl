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
    ReactionUOceanfloor_dev

Calculate U ocean burial in 'anoxic' and 'other' sinks with partitioning linearly dependent on 'ocean anoxia' U_anoxia.

Anoxia partitioning is a linear map:

    U_ANOX = k_anox_map_0 + ocean.ANOX * k_anox_map_1

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionUOceanfloor_dev{P} <: PB.AbstractReaction
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
        PB.ParDouble("k_anox_map_0", 0.0, units="",
            description="actual anoxic area fraction U_ANOX = k_anox_map_0 + ocean.ANOX * k_anox_map_1"),
        PB.ParDouble("k_anox_map_1", 1.0, units="",
            description="actual anoxic area fraction U_ANOX = k_anox_map_0 + ocean.ANOX * k_anox_map_1"),

        PB.ParType(PB.AbstractData, "UIsotope", PB.IsotopeLinear,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable uranium isotopes and specify isotope type"),
    )

end


function PB.register_methods!(rj::ReactionUOceanfloor_dev)

    UIsotopeType = rj.pars.UIsotope[]
    PB.setfrozen!(rj.pars.UIsotope)

    vars = [
        PB.VarDepScalar("ocean.ANOX", "",  "ocean anoxic fraction"),
        PB.VarDepScalar("ocean.U_norm", "",  "normalized ocean U"),
        PB.VarDepScalar("ocean.U_delta", "",  "ocean d238U"),
    
        PB.VarPropScalar("U_ANOX", "",  "linearly mapped ocean anoxic fraction"),
        PB.VarPropScalar("U_anoxic", "mol U yr-1",  "anoxic sink flux";
            attributes=(:field_data=>UIsotopeType,)),
        PB.VarPropScalar("U_other", "mol U yr-1",  "other sink flux";
            attributes=(:field_data=>UIsotopeType,)),
    
        PB.VarContrib("solutefluxOceanfloor_U"=>"fluxOceanfloor.soluteflux_U", "mol yr-1",  "U oceanfloor solute flux";
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

    vars.U_ANOX[] = clamp(pars.k_anox_map_0[] + vars.ANOX[] * pars.k_anox_map_1[], 0.0, 1.0) # linear map

    vars.U_anoxic[] = @PB.isotope_totaldelta(
        UIsotopeType, 
        pars.k_U_anoxic[]*vars.U_norm[]*(vars.U_ANOX[]/pars.k_anox_0[]),
        vars.U_delta[] + pars.k_U_anoxic_D[]
    )

    vars.U_other[]  = @PB.isotope_totaldelta(
        UIsotopeType, 
        pars.k_U_other[]*vars.U_norm[]*(1.0 - vars.U_ANOX[])/(1.0 - pars.k_anox_0[]),
        vars.U_delta[] + pars.k_U_other_D[],
    )
    
    r_nfloorcells = 1.0/PB.get_length(rj.domain) # fraction of flux for each oceanfloor cell
    for i in cellrange.indices 
        vars.solutefluxOceanfloor_U[i] -= r_nfloorcells*(vars.U_anoxic[] + vars.U_other[]) 
    end

    return nothing
end


end
