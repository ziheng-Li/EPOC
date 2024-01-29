module ReactionsOOEOAE

import PALEOboxes as PB
using PALEOboxes.DocStrings
# import PALEOreactions
import PALEOocean
import PALEOaqchem

#####################################################
# Atmosphere
#####################################################

"""
    ReactionOceanBurialColumn

Minimal ocean model representing organic carbon and phosphorus burial
distributed over columns with varying oceanfloor oxygen level.

Assumes a scalar state Variable `O` for normalized atmosphere-ocean oxygen level `O_norm`, and a scalar state Variable `P` for 
total ocean phosphorus with normalized level `P_norm`.

Oxygen utilisation is assumed to scale linearly with `P` and oxygen availability with `O`, so normalized seafloor oxygen
`O_norm_sf` in column `i` is given by:

    O_norm_sf[i] = O_norm - P_norm * O2_U[i] 

where the normalized oxygen utilisation in each column `O2_U` is set by parameters `k_O2_U_min` and `k_O2_U_max`.

Organic carbon burial is distributed equally across the columns, and phosphorus burial is determined from an oxygen-dependent
Corg:P burial ratio.

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionOceanBurialColumn{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple([
        PB.ParInt("ncols", 100,
            description="number of ocean columns"),
        PB.ParDouble("k_mocb", 4.5e12, units="mol yr-1",
            description="marine organic carbon burial rate"),
        PB.ParDouble("corg_burial_fac", 1.0, units="",
            description="multiplier for marine organic carbon burial rate (and C:P burial ratio)"),
        PB.ParDouble("k_anox", 1000.0,
            description="anoxia function sharpness"),
        PB.ParDouble("k_oxic", 250.0,
            description="local C:P burial ratio oxic"),
        PB.ParDouble("k_anoxic", 1000.0, 
            description="local C:P burial ratio anoxic"),
        PB.ParDouble("k_O2_U_min", 0.49,
            description="minimum normalized O2 utilisation"),
        PB.ParDouble("k_O2_U_max", 0.51,
            description="maximum normalized O2 utilisation"),
        PB.ParString("k_O2_U_rate", "linear", # allowed_values=["linear", "power0.25", "power4"],
            description="calculation method of the O2_U per cell, 
                        for example the default 'linear' use the linear interpolation within the k_O2_U_min and _max"),
        PB.ParDouble("anoxic_threshold", 0.9975, # 0.997527 in COPSE
            description="anoxic threshold"),

        # # Uranium isotopes
        # # U0 see ReactionReservoirScalar U

        # PB.ParDouble("f_anoxic0",      0.0025, units="",          description="ocean anoxia present day anoxic fraction"),
        # PB.ParDouble("k_U_anoxic",     6e6,   units="mol U/yr",   description="anoxic sink"),
        # PB.ParDouble("k_U_other",      34e6,  units="mol U/yr",   description="other sinks combined"),
        # PB.ParDouble("Delta_U_anoxic", 0.6,   units="per mil",    description="Anoxic sink fractionation"),
        # PB.ParDouble("Delta_U_other",  0.005, units="per mil",    description="Other sinks fractionation"),

        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
        PB.ParType(PB.AbstractData, "UIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable uranium isotopes and specify isotope type"),
    ])
end


function PB.set_model_geometry(rj::ReactionOceanBurialColumn, model::PB.Model)
    
    ocean_cells = rj.pars.ncols[] # Number of cells (= ocean Domain size)

    # # set minimal grid 
    # # Ziheng fixed check_length issue
    # isurf=[1]
    # surfacegrid = PB.Grids.UnstructuredVectorGrid(ncells=length(isurf))
    # # surfacegrid = PB.Grids.UnstructuredVectorGrid(ncells=ocean_cells)
    # oceangrid = PB.Grids.UnstructuredVectorGrid(ncells=ocean_cells)
    # floorgrid = PB.Grids.UnstructuredVectorGrid(ncells=ocean_cells)

    # set minimal grid for 100 columns
    ocean_cells = rj.pars.ncols.v # Number of cells (= ocean Domain size)
    isurf = collect(range(1, length=ocean_cells))   # [1]
    ifloor = collect(range(1, length=ocean_cells))
    oceangrid = PB.Grids.UnstructuredVectorGrid(ncells=ocean_cells) # 
    rj.domain.grid = oceangrid
    PB.Grids.set_subdomain!(oceangrid, "oceansurface", PB.Grids.BoundarySubdomain(isurf), true)
    @info "  set ocean.oceansurface Subdomain size=$(length(isurf))"
    PB.Grids.set_subdomain!(oceangrid, "oceanfloor", PB.Grids.BoundarySubdomain(ifloor), true)
    @info "  set ocean.oceanfloor Subdomain size=$(length(ifloor))"
    
    surfacegrid = PB.Grids.UnstructuredVectorGrid(ncells=length(isurf))
    PB.Grids.set_subdomain!(surfacegrid, "ocean", PB.Grids.InteriorSubdomain(ocean_cells, isurf), true)

    floorgrid = PB.Grids.UnstructuredVectorGrid(ncells=length(ifloor))
    PB.Grids.set_subdomain!(floorgrid, "ocean", PB.Grids.InteriorSubdomain(ocean_cells, ifloor), true)

    PALEOocean.Ocean.set_model_domains(model, oceangrid, surfacegrid, floorgrid)    
    
    return nothing
end


function PB.register_methods!(rj::ReactionOceanBurialColumn)

    CIsotopeType = rj.pars.CIsotope[]
    UIsotopeType = rj.pars.UIsotope[]

    vars_setup = [
        PB.VarPropStateIndep("O2_U", "", "normalized O2 utilisation"),
        PB.VarPropStateIndep("frac_Corg_burial", "", "fraction of total Corg burial per column"),
    ]

    PB.add_method_setup!(
        rj,
        setup_ocean_burial_column,
        (PB.VarList_namedtuple(vars_setup), )
    )

    vars_do = [
        PB.VarDep.(vars_setup)...,
        PB.VarDepScalar("P_norm", "", "normalized marine phosphorus"),
        PB.VarDerivScalar("P_sms", "mol yr-1", "marine phosphorus source - sink"),
        PB.VarDerivScalar("(DIC_sms)", "mol yr-1", "marine DIC source - sink",          # NB: DIC is optional
            attributes=(:field_data=>CIsotopeType,)), 
        PB.VarDepScalar("atmocean.O_norm", "", "normalized atm-ocean oxygen"),
        PB.VarDerivScalar("atmocean.O_sms", "mol yr-1", "atm-ocean oxygen source - sink"),
        PB.VarProp("oceanfloor.local_anoxia", "", "anoxia 0 (oxic) - 1 (anoxic)"),
        PB.VarProp("oceanfloor.local_CPsea", "", "Corg:P burial ratio"),
        PB.VarPropScalar("oceanfloor.anoxia_burial_frac", "", "fraction of Corg burial in anoxic condition")
    ]

    if CIsotopeType <: PB.AbstractIsotopeScalar
        append!(vars_do, [
                # C isotopes
                PB.VarProp("ocean.mocb_delta", "per mil", "D13C of marine organic carbon"), 

                PB.VarDepScalar("ocean.DIC_delta",  "per mil",  "d13C ocean DIC"),
                PB.VarDepScalar("ocean.D_mccb_DIC",  "per mil",  "D13C marine calcite burial relative to ocean DIC"),
                PB.VarDepScalar("ocean.D_B_mccb_mocb",     "per mil",      "D13C fractionation between marine organic and calcite burial"),
            ]
        )
    end

    # if UIsotopeType <: PB.AbstractIsotopeScalar
    #     append!(vars_do, [
    #             # U isotopes
    #             PB.VarDepScalar("U",        "mol U",        "ocean uranium",
    #                 attributes=(:field_data=>UIsotopeType,)),
    #             PB.VarContribScalar("U_sms","mol U yr-1",   "ocean uranium source - sink",
    #                 attributes=(:field_data=>UIsotopeType,)),
    #             PB.VarDepScalar("U_norm",   "",             "normalized ocean uranium"),
    #             PB.VarDepScalar("U_delta",  "per mil",      "ocean uranium fractionation"),
    #             PB.VarProp("F_U_anoxic","mol U yr-1",   "Anoxic U sink",
    #                 attributes=(:field_data=>UIsotopeType,)),
    #             PB.VarProp("F_U_other", "mol U yr-1",   "Other U sinks",
    #                 attributes=(:field_data=>UIsotopeType,)),
    #         ]
    #     )
    # end

    fluxOceanBurial = PB.Fluxes.FluxContrib(
        "fluxOceanBurial.flux_", ["Corg::$CIsotopeType", "P", "U_anoxic::$UIsotopeType","U_other::$UIsotopeType"],
    )

    PB.add_method_do!(
        rj,
        do_ocean_burial_column,
        (PB.VarList_namedtuple(vars_do),  PB.VarList_namedtuple_fields(fluxOceanBurial),),
        p = (CIsotopeType, UIsotopeType), # provide isotope types here so Julia will generate specialized (fast) code
    )
    
    return nothing
end

function setup_ocean_burial_column(m::PB.ReactionMethod, (vars, ), cellrange::PB.AbstractCellRange, attribute_name)
    rj = m.reaction
    
    attribute_name == :setup || return

    k_O2_U_rate = rj.pars.k_O2_U_rate[]

    if k_O2_U_rate == "linear"
        vars.O2_U .=  range(rj.pars.k_O2_U_min[], rj.pars.k_O2_U_max[], rj.domain.grid.ncells)   
    elseif k_O2_U_rate[1:5] == "power" 
        rate = parse(Float64, k_O2_U_rate[6:end])
        vars.O2_U .= ((1:rj.domain.grid.ncells) ./ rj.domain.grid.ncells) .^ rate .* (rj.pars.k_O2_U_max[] - rj.pars.k_O2_U_min[]) .+ rj.pars.k_O2_U_min[]
    else
        @error "Wrong input of the k_O2_U_rate=$(k_O2_U_rate), allowed_values=[linear, power*]"
    end 

    vars.frac_Corg_burial .= 1/rj.domain.grid.ncells

    return nothing
end

function do_ocean_burial_column(m::PB.ReactionMethod, (vars, fluxBurial), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction
    CIsotopeType, UIsotopeType = m.p

    newp_n = newp_norm(vars.P_norm[])
    mocb_n = mocb_norm(newp_n)
    mocb = mocb_n * rj.pars.k_mocb[] * rj.pars.corg_burial_fac[]

    vars.anoxia_burial_frac[] = 0.0

    for i in cellrange.indices

        local_Corgburial = mocb*vars.frac_Corg_burial[i]

        if CIsotopeType <: PB.AbstractIsotopeScalar 
            # delta of marine organic carbon burial
            vars.mocb_delta[i]     = vars.DIC_delta[] + vars.D_mccb_DIC[] - vars.D_B_mccb_mocb[]
        else
        end       
        
        mocb_isotope = @PB.isotope_totaldelta(CIsotopeType, local_Corgburial, vars.mocb_delta[i])

        vars.local_anoxia[i] = local_anoxia(newp_n, vars.O_norm[], vars.O2_U[i], rj.pars.k_anox.v)
        vars.local_CPsea[i] = CPsea(vars.local_anoxia[i], rj.pars.k_oxic.v, rj.pars.k_anoxic.v)*rj.pars.corg_burial_fac.v

        # TODO debugging AD Jacobian
        # if !isa(vars.local_anoxia[i], Float64)
        #    @info "i: $i local_anoxia $(vars.local_anoxia[i])  local_CPsea $(vars.local_CPsea[i])"
        # end

        if vars.local_anoxia[i] > rj.pars.anoxic_threshold[]
            vars.anoxia_burial_frac[] += vars.frac_Corg_burial[i]  # per cell
        end

        vars.O_sms[] += local_Corgburial 
        fluxBurial.Corg[i] += mocb_isotope
        PB.add_if_available(vars.DIC_sms, -mocb_isotope)

        local_Pburial = local_Corgburial/vars.local_CPsea[i]
        fluxBurial.P[i] += local_Pburial
        vars.P_sms[] -= local_Pburial
    end

    return nothing
end

# Total marine new production (normalized) as a function of normalized P
function newp_norm(Pnorm)
    newp_norm = Pnorm
    return newp_norm
end

# Total marine organic carbon burial as a function of new production (assume linear)
function mocb_norm(newp_norm)
    mocb_norm = newp_norm
    return mocb_norm
end

# Local marine anoxia function 
# Switch 0 - 1, sharpness controlled by k_anox  (from COPSE Reloaded)
function local_anoxia(newp_norm, Onorm, k_U, k_anox)
    # control variable for anoxia
    v_anox = k_U*newp_norm - Onorm
    
    # anox = 1.0/(1.0 + exp(-k_anox*v_anox))
    # TODO limit to give small value not zero, so AD Jacobian works 
    anox = 1.0/(1.0 + min(1e80, exp(-k_anox*v_anox)))
    return anox
end

# Local marine C:P burial ratio from COPSE
function CPsea(anox, k_oxic, k_anoxic)
    CPsea = k_oxic*k_anoxic/((1.0-anox)*k_anoxic + anox*k_oxic)
    return CPsea
end



"""
    ReactionOxWeathMinimal

Minimal land surface oxidative weathering

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionOxWeathMinimal{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple([
        PB.ParDouble("k_oxidw", 3.0e12, units="mol yr-1",
            description="organic carbon oxidation rate at O_norm=1.0"),    
        PB.ParBool("A_oxidw", false, description="true to link the oxidw to atmocean.A as an input"),     

        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
    ])
end


function PB.register_methods!(rj::ReactionOxWeathMinimal)

    CIsotopeType = rj.pars.CIsotope.v

    vars = [        
        PB.VarDepScalar("atmocean.O_norm", "", "normalized atm-ocean oxygen"),        
        PB.VarProp("oxweath", "mol yr-1", "organic carbon oxidation rate")
    ]

    fluxAtoLand = PB.Fluxes.FluxContrib(
        "fluxAtoLand.flux_", ["CO2::$CIsotopeType", "O2"],
    )

    PB.add_method_do!(
        rj,
        do_oxweath_minimal,
        (PB.VarList_namedtuple(vars),  PB.VarList_namedtuple_fields(fluxAtoLand),)
    )
    
    return nothing
end

function do_oxweath_minimal(m::PB.ReactionMethod, (vars, fluxAtoLand), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction
    CIsotopeType = rj.pars.CIsotope[]

    vars.oxweath[] = rj.pars.k_oxidw.v * sqrt(max(vars.O_norm[], 0.0))
    fluxAtoLand.O2[] += vars.oxweath[]    

    if rj.pars.A_oxidw.v
        # hardcode I guess the oxidw C with -5 isotopic value, need to fix later
        fluxAtoLand.CO2[] -= @PB.isotope_totaldelta(CIsotopeType, vars.oxweath[] , -5)
    end

    return nothing
end

"""
    ReactionUOceanfloor_dev

Calculate U ocean burial in 'anoxic' and 'other' sinks with partitioning linearly dependent on 'ocean anoxia'.

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
    for i in cellrange.indices 
        U_anoxic_sink   = @PB.isotope_totaldelta(UIsotopeType, vars.U_anoxic[], vars.U_delta[] + pars.k_U_anoxic_D[])
        U_other_sink    = @PB.isotope_totaldelta(UIsotopeType, vars.U_other[],  vars.U_delta[] + pars.k_U_other_D[])
        vars.solutefluxOceanfloor_U[i] -= r_nfloorcells*(U_anoxic_sink + U_other_sink) 
    end

    return nothing
end


end # module