module ReactionsOOEOAE

import PALEOboxes as PB
using PALEOboxes.DocStrings
# import PALEOreactions
import PALEOaqchem
import PALEOocean

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
        PB.ParBool("use_shelf_area_norm", false,
            description="true to multiply Corg and P burial fluxes by a factor read from Variable global.SHELF_AREA_NORM"),
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

    # set minimal grid for ncols columns
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
        PB.VarProp("O2_U_local", "", "normalized O2 utilisation multiplied by kO2U"),
        PB.VarPropStateIndep("frac_Corg_burial", "", "fraction of total Corg burial per column"),
    ]

    PB.add_method_setup!(
        rj,
        setup_ocean_burial_column,
        (PB.VarList_namedtuple(vars_setup), )
    )

    vars_do = [
        PB.VarDep.(vars_setup)...,
        PB.VarDepScalar("corg_bf", "", "multiplier for corg_burial_fac"), # the multiplier of a mutiplier
        PB.VarDepScalar("kO2U", "", "multiplier for k_O2_U_min"), # fix the oxic fold point, change the anoxic fold point to control the sharpness
        PB.VarDepScalar("P_norm", "", "normalized marine phosphorus"),
        PB.VarContribScalar("P_sms", "mol yr-1", "marine phosphorus source - sink"),
        PB.VarContribScalar("(DIC_sms)", "mol yr-1", "marine DIC source - sink",          # NB: DIC is optional
            attributes=(:field_data=>CIsotopeType,)), 
        PB.VarDepScalar("atmocean.O_norm", "", "normalized atm-ocean oxygen"),
        PB.VarContribScalar("atmocean.O_sms", "mol yr-1", "atm-ocean oxygen source - sink"),
        PB.VarProp("oceanfloor.local_anoxia", "", "anoxia 0 (oxic) - 1 (anoxic)"),
        PB.VarProp("oceanfloor.local_CPsea", "", "Corg:P burial ratio"),
        PB.VarPropScalar("oceanfloor.anoxia_burial_frac", "", "fraction of Corg burial in anoxic condition")
    ]

    if rj.pars.use_shelf_area_norm[]
        push!(vars_do, PB.VarDepScalar("global.SHELF_AREA_NORM", "", "normalized shelf area forcing for marine Corg and P burial"))
    end
    PB.setfrozen!(rj.pars.use_shelf_area_norm) # can't be changed as needs a new Variable

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
        "fluxOceanBurial.flux_", ["Corg::$CIsotopeType", "P"],
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

function do_ocean_burial_column(m::PB.ReactionMethod, pars, (vars, fluxBurial), cellrange::PB.AbstractCellRange, deltat)
   
    CIsotopeType, UIsotopeType = m.p
    rj = m.reaction

    newp_n = newp_norm(vars.P_norm[])
    mocb_n = mocb_norm(newp_n)

    shelf_area_norm = pars.use_shelf_area_norm[] ? vars.SHELF_AREA_NORM[] : 1.0

    mocb = mocb_n * pars.k_mocb[]  * shelf_area_norm * 
                            (pars.corg_burial_fac[] * vars.corg_bf[])

    vars.anoxia_burial_frac[] = 0.0

    # Update the vars.O2_U, let the O2_U can change along the time
    # fix the oxic fold point, change the anoxic fold point to control the sharpness
    k_O2_U_min = vars.O2_U[1]*vars.kO2U[]
    k_O2_U_max = vars.O2_U[end]
    vars.O2_U_local .=  range(k_O2_U_min, k_O2_U_max, rj.domain.grid.ncells)
    # @info "k_O2_U_min=$(k_O2_U_min)"

    for i in cellrange.indices

        local_Corgburial = mocb*vars.frac_Corg_burial[i]

        if CIsotopeType <: PB.AbstractIsotopeScalar 
            # delta of marine organic carbon burial
            vars.mocb_delta[i]     = vars.DIC_delta[] + vars.D_mccb_DIC[] - vars.D_B_mccb_mocb[]
        else
        end       
        
        mocb_isotope = @PB.isotope_totaldelta(CIsotopeType, local_Corgburial, vars.mocb_delta[i])

        vars.local_anoxia[i] = local_anoxia(newp_n, vars.O_norm[], vars.O2_U_local[i], pars.k_anox[])
        vars.local_CPsea[i] = CPsea(vars.local_anoxia[i], pars.k_oxic[], pars.k_anoxic[])* (pars.corg_burial_fac[] * vars.corg_bf[])

        # TODO debugging AD Jacobian
        # if !isa(vars.local_anoxia[i], Float64)
        #    @info "i: $i local_anoxia $(vars.local_anoxia[i])  local_CPsea $(vars.local_CPsea[i])"
        # end

        if vars.local_anoxia[i] > pars.anoxic_threshold[]
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

# """
#     ReactionUOceanfloor_dev

# Calculate U ocean burial in 'anoxic' and 'other' sinks with partitioning linearly dependent on 'ocean anoxia' U_anoxia.

# Anoxia partitioning is a linear map:

#     U_ANOX = k_anox_map_0 + ocean.ANOX * k_anox_map_1

# # Parameters
# $(PARS)

# # Methods and Variables
# $(METHODS_DO)
# """
# Base.@kwdef mutable struct ReactionUOceanfloor_dev{P} <: PB.AbstractReaction
#     base::PB.ReactionBase

#     pars::P = PB.ParametersTuple(
#         PB.ParDouble("k_U_anoxic_D", 0.6, units="per mil",
#             description="anoxic sink fractionation"),
#         PB.ParDouble("k_U_other_D", 0.005, units="per mil",
#             description="other sink fractionation"),
    
#         PB.ParDouble("k_U_anoxic", 6.0e6, units="mol U yr-1",
#             description="present-day anoxic sink flux"),
#         PB.ParDouble("k_U_other", 34.0e6, units="mol U yr-1",
#             description="present-day other sink flux"),
    
#         PB.ParDouble("k_anox_0", 0.0025, units="",
#             description="present-day ocean anoxia fraction"),
#         PB.ParDouble("k_anox_map_0", 0.0, units="",
#             description="actual anoxic area fraction U_ANOX = k_anox_map_0 + ocean.ANOX * k_anox_map_1"),
#         PB.ParDouble("k_anox_map_1", 1.0, units="",
#             description="actual anoxic area fraction U_ANOX = k_anox_map_0 + ocean.ANOX * k_anox_map_1"),

#         PB.ParType(PB.AbstractData, "UIsotope", PB.IsotopeLinear,
#             external=true,
#             allowed_values=PB.IsotopeTypes,
#             description="disable / enable uranium isotopes and specify isotope type"),
#     )

# end


# function PB.register_methods!(rj::ReactionUOceanfloor_dev)

#     UIsotopeType = rj.pars.UIsotope[]
#     PB.setfrozen!(rj.pars.UIsotope)

#     vars = [
#         PB.VarDepScalar("ocean.ANOX", "",  "ocean anoxic fraction"),
#         PB.VarDepScalar("ocean.U_norm", "",  "normalized ocean U"),
#         PB.VarDepScalar("ocean.U_delta", "",  "ocean d238U"),
    
#         PB.VarPropScalar("U_ANOX", "",  "linearly mapped ocean anoxic fraction"),
#         PB.VarPropScalar("U_anoxic", "mol U yr-1",  "anoxic sink flux";
#             attributes=(:field_data=>UIsotopeType,)),
#         PB.VarPropScalar("U_other", "mol U yr-1",  "other sink flux";
#             attributes=(:field_data=>UIsotopeType,)),
    
#         PB.VarContrib("solutefluxOceanfloor_U"=>"fluxOceanfloor.soluteflux_U", "mol yr-1",  "U oceanfloor solute flux";
#             attributes=(:field_data=>UIsotopeType,))   
#     ]

#     PB.add_method_do!(
#         rj, 
#         do_U_oceanfloor,
#         (PB.VarList_namedtuple(vars), ),
#         p = UIsotopeType,
#     )

#     return nothing
# end

# function do_U_oceanfloor(
#     m::PB.ReactionMethod,
#     pars,
#     (vars, ), 
#     cellrange::PB.AbstractCellRange,
#     deltat
# )
#     rj = m.reaction
#     UIsotopeType = m.p

#     vars.U_ANOX[] = clamp(pars.k_anox_map_0[] + vars.ANOX[] * pars.k_anox_map_1[], 0.0, 1.0) # linear map

#     vars.U_anoxic[] = @PB.isotope_totaldelta(
#         UIsotopeType, 
#         pars.k_U_anoxic[]*vars.U_norm[]*(vars.U_ANOX[]/pars.k_anox_0[]),
#         vars.U_delta[] + pars.k_U_anoxic_D[]
#     )

#     vars.U_other[]  = @PB.isotope_totaldelta(
#         UIsotopeType, 
#         pars.k_U_other[]*vars.U_norm[]*(1.0 - vars.U_ANOX[])/(1.0 - pars.k_anox_0[]),
#         vars.U_delta[] + pars.k_U_other_D[],
#     )
    
#     r_nfloorcells = 1.0/PB.get_length(rj.domain) # fraction of flux for each oceanfloor cell
#     for i in cellrange.indices 
#         vars.solutefluxOceanfloor_U[i] -= r_nfloorcells*(vars.U_anoxic[] + vars.U_other[]) 
#     end

#     return nothing
# end



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


"""
    ReactionLandWeatheringFluxes_dev

COPSE Reloaded(2018) weathering fluxes.
Calculates and applies global weathering fluxes and land burial given weathering rates and land areas.

Fluxes are added to flux couplers:
- `fluxAtoLand`:  CO2 and O2 exchange with atmosphere
- `fluxRtoOcean`: riverine fluxes
- `fluxLandtoSedCrust`: sedimentary reservoir weathering and land organic carbon burial

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionLandWeatheringFluxes_dev{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        # silicate weathering
        PB.ParString("f_gran_link_u",       "original",
            allowed_values=["original", "weak", "none"],
            description="granite weathering uplift dependence"),
        PB.ParString("f_bas_link_u",       "yes",
            allowed_values=["yes", "weak", "no"],
            description="basalt weathering uplift dependence"),
        
        PB.ParDouble("k_basfrac",    0.35,   units="",
            description="fraction of silicate weathering that is basaltic at present"),

        PB.ParDouble("k_granw",      NaN,   units="mol/yr",
            description="granite weathering rate"),
        PB.ParDouble("k_basw",      NaN,   units="mol/yr",
            description="basalt weathering rate"),

        # P weathering
        PB.ParString("f_p_kinetics",       "no",
            allowed_values=["no", "yes"],
            description="apatite kinetics independent of host rock"),
        PB.ParString("f_p_apportion",       "no",
            allowed_values=["no", "yes"],
            description="P content vary between granite and basalt"),
        PB.ParDouble("k_P",        2.15,
            description="enrichment P:(Ca+Mg) in granite vs basalt "), 
        PB.ParDouble("k10_phosw",      NaN,        units="mol/yr",
            description="phosphorus weathering"),
        PB.ParDouble("k_Psilw",     2/12,
            description="fraction of P weathering from silicates"),
        PB.ParDouble("k_Pcarbw",     5/12,
            description="fraction of P weathering associated with carbonates"),
        PB.ParDouble("k_Poxidw",     5/12,
            description="fraction of P weathering from organic matter"),
        PB.ParDouble("k_Psedw",      0.0,
            description="fraction of P weathering from other sedimentary rocks"),

        # Carbonate weathering
        PB.ParString("f_carb_link_u",       "yes",
            allowed_values=["yes", "no"],
            description="carbonate weathering uplift dependence"),
        PB.ParString("f_carbwC",       "Cindep",
            allowed_values=["Cindep", "Cprop"],
            description="C (carbonate reservoir) dependence of carbonate weathering"),
        PB.ParDouble("k14_carbw",      13.35e12,   units="mol/yr",
            description="carbonate weathering"),

        # oxidative weathering
        PB.ParString("f_oxwG",       "Gprop",
            allowed_values=["Gindep", "Gprop", "forced"],
            description="G (organic carbon reservoir) dependence of oxidative weathering"),
        PB.ParDouble("k17_oxidw",      NaN,        units="mol/yr",
            description="oxidative org carbon weathering"),

        # Sulphur weathering
        PB.ParBool("enableS",          true,
            description="enable S weathering"),
        PB.ParString("f_pyrweather",    "copse_O2",
            allowed_values=["copse_O2", "copse_noO2", "forced"],
            description="functional form of dependence of pyrite weathering"),
        PB.ParDouble("k21_pyrw",       0.53e12,    units="mol S/yr",
            description="pyrite weathering"),
        PB.ParString("f_gypweather",    "original",
            allowed_values=["original", "alternative", "forced"],
            description="functional form of dependence of gypsum weathering"),
        PB.ParDouble("k22_gypw",       1e12,       units="mol S/yr",
            description="gypsum weathering"),
        PB.ParBool("SRedoxAlk",     false,
            description="true to include -TAlk from pyrite weathering"),

        # burial fluxes
        PB.ParString("f_locb",    "original",
            allowed_values=["original", "Uforced", "coal", "split", "Prescribed"],
            description="functional form of dependence of land organic carbon burial"),
        PB.ParDouble("k5_locb",        4.5e12,     units="mol/yr",
            description="land organic carbon burial (f_locb==Prescribed only)"),
        PB.ParDouble("k11_landfrac",   0.10345,    units="",
            description="fract of weath P buried on land"),
        PB.ParDouble("k_aq",           0.8,       units="",
            description="fraction of locb assumed to occur in aquatic settings (not coals)"),
        PB.ParDouble("CPland0",        1000.0,     units="",
            description="present-day C/P land burial"),

        # Isotopes
        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
        PB.ParType(PB.AbstractData, "SIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable sulphur isotopes and specify isotope type"),
    )

end


function PB.register_methods!(rj::ReactionLandWeatheringFluxes_dev)
 
    # dependencies
    vars_dep = PB.VarVector(PB.VarDepScalar,
        [
            # Forcings
            ("global.tforce",  "yr"   , "time for external forcings"),
            ("global.UPLIFT",  "",      "uplift scaling"),
            ("global.PG",       "",      "paleogeographic forcing"),        
            ("global.RHOSIL",     "",  "silicate weathering enhancement"),
            ("global.RHO",     "",  "carbonate weathering enhancement"),
            ("global.F_EPSILON",     "",  "phosphorus weathering enhancement"),
            ("global.CPland_relative","","land CP burial ratio scaling"),
            ("(global.COAL)",     "",    "coal forcing"),

            # ziheng's new forcing, koxidw, a multiplier to increase the oxidw
            ("koxidw", "", "multiplier to increase the oxidw"),

            # Land areas
            ("GRAN_AREA",     "",   "granite area"),
            ("BA_AREA",     "",     "basalt area"),
            ("CARB_AREA",     "",   "carbonate area"),
            ("(ORG_AREA)",     "",    "organic carbon area"),       
            # Isotopes
            ("(D_P_CO2_locb)", "per mil",  "d13C fractionation between terrestrial organic burial and atmospheric CO2"),
            ("(D_eqbw_CO2)",    "per mil", "d13C fractionation between atmospheric CO2 and fresh water"),
            ("(atm.CO2_delta)", "per mil",  "atmospheric pCO2 delta 13C"),
            # biota
            ("VEG",     "",           "Mass of terrestrial biosphere"),  
            # rates
            ("f_gran",          "",         "granite weathering rate"),
            ("f_bas",          "",         "basalt weathering rate"),
            ("f_ap",          "",         "apatite weathering rate"),
            ("f_carb",          "",         "carbonate weathering rate"),
            ("oxw_facO",        "",         "oxidative weathering oxygen dependence factor"),
            ("oxw_facOmax",        "",      "maximum possible oxidative weathering factor at high pO2"),
            # sedcrust
            ("sedcrust.C_norm",        "",    "Sedimentary carbonate normalized to present day"),
            ("(sedcrust.C_delta)",        "",    "Sedimentary carbonate d13C"),
            ("sedcrust.G_norm",        "",    "Sedimentary organic carbon normalized to present day"),
            ("(sedcrust.G_delta)",        "",    "Sedimentary organic carbon d13C"),
            ("(sedcrust.GYP_norm)",      "",    "Sedimentary gypsum normalized to present day"),
            ("(sedcrust.GYP_delta)",      "",    "Sedimentary gypsum d34S"),
            ("(sedcrust.PYR_norm)",      "",    "Sedimentary pyrite normalized to present day"),
            ("(sedcrust.PYR_delta)",      "",    "Sedimentary pyrite d34S")
        ]
    )

    # Properties we calculate in do_react
    vars_prop = PB.VarVector(PB.VarPropScalar,
        [
            ("granw_relative",    "",    "granite weathering normalized to present"),
            ("granw",    "molC/yr",    "granite weathering"),
            ("basw_relative",    "",    "basalt weathering normalized to present"),
            ("basw",    "molC/yr",    "basalt weathering"),
            ("silw",    "molC/yr",    "total silicate weathering"),
            ("silw_relative", "",     "total silicate weathering normalized to present"),

            ("carbw_relative", "",    "Carbonate weathering normalized to present"),
            ("carbw_fac", "",    "Carbonate weathering normalized to present, without C reservoir dependence"),
            ("carbw",   "molC/yr",    "Carbonate weathering"),
            ("silwcarbw_relative",    "",   "total silicate and carbonate weathering, normalized to present"),
        
            ("orgcw",   "molC/yr",    "organic C erosion"),       
            ("oxidw",   "molC/yr",    "Oxidative organic C weathering"),
            ("gypw",    "molS/yr",    "Gypsum weathering"),
            ("pyrw",    "molS/yr",    "Pyrite weathering"),
            ("granw_ap",    "molC/yr",  "apatite effective granite weathering flux"),      
            ("basw_ap",    "molC/yr",   "apatite effective basalt weathering flux"),
            ("silw_ap",    "molC/yr",   "apatite effective total silicate weathering flux"),
            ("sedw_relative", "",     "Other sediment weathering (sandstones etc) normalized to present"),
            ("phosw_x", "molP/yr",    "Phosphorus weathering, other sediment associated"),
            ("phosw_s", "molP/yr",    "Phosphorus weathering, silicate-associated"),
            ("phosw_c", "molP/yr",    "Phosphorus weathering, carbonate-associated"),
            ("phosw_o", "molP/yr",    "Phosphorus weathering, organic carbon-associated"),
            ("phosw",   "molP/yr",    "Phosphorus weathering, total"),
            ("pland",   "molP/yr",    "Phosphorus weathering to land"),
            ("psea",    "molP/yr",    "Phosphorus weathering to sea"),
            # burial
            ("locb",   "molC/yr",     "Land organic carbon burial"),
            ("locb_delta",   "per mil","Land organic carbon burial d13C")
        ]
    )

    # isotope Types
    CIsotopeType = rj.pars.CIsotope[]
    SIsotopeType = rj.pars.SIsotope[]

    # Add flux couplers
    fluxAtoLand = PB.Fluxes.FluxContribScalar(
        "fluxAtoLand.flux_", ["CO2::$CIsotopeType", "O2"],
        isotope_data=Dict())

    fluxRtoOcean = PB.Fluxes.FluxContribScalar(
        "fluxRtoOcean.flux_", ["DIC::$CIsotopeType", "TAlk", "Ca", "P", "SO4::$SIsotopeType"],
        isotope_data=Dict())

    fluxLandtoSedCrust = PB.Fluxes.FluxContribScalar(
        "fluxLandtoSedCrust.flux_", ["Ccarb::$CIsotopeType", "Corg::$CIsotopeType", "GYP::$SIsotopeType", "PYR::$SIsotopeType"],
        isotope_data=Dict())

    PB.add_method_do!(
        rj,
        do_land_weathering_fluxes,
        (   
            PB.VarList_namedtuple_fields(fluxAtoLand),
            PB.VarList_namedtuple_fields(fluxRtoOcean),
            PB.VarList_namedtuple_fields(fluxLandtoSedCrust),
            PB.VarList_namedtuple([vars_dep; vars_prop]),
        ),
        p=(CIsotopeType, SIsotopeType),
    )

    return nothing
end

function do_land_weathering_fluxes(
    m::PB.ReactionMethod,
    pars,
    (fluxAtoLand, fluxRtoOcean, fluxLandtoSedCrust, D), 
    cellrange::PB.AbstractCellRange,
    deltat
)
    (CIsotopeType, SIsotopeType) = m.p

    # granite weathering  
    if pars.f_gran_link_u[] == "original"
        granw_fac_u = D.UPLIFT[]
    elseif pars.f_gran_link_u[] == "weak"
        granw_fac_u = (D.UPLIFT[]^0.33)
    elseif pars.f_gran_link_u[] == "none"
        granw_fac_u = 1.0
    else
        error("unrecognized f_gran_link_u ", pars.f_gran_link_u[])
    end
    D.granw_relative[]  = granw_fac_u * D.PG[]*D.RHOSIL[] * D.GRAN_AREA[] * D.f_gran[]
    D.granw[]           = D.granw_relative[] * pars.k_granw[]

    # basalt weathering
    if pars.f_bas_link_u[] == "no"
        basw_fac_u = 1.0
    elseif pars.f_bas_link_u[] == "yes"
        basw_fac_u = D.UPLIFT[]
    elseif pars.f_bas_link_u[] == "weak"
        basw_fac_u = D.UPLIFT[]^0.33
    else
        error("unrecognized f_bas_link_u ", pars.f_bas_link_u[])
    end
    D.basw_relative[]   = basw_fac_u * D.PG[]*D.RHOSIL[]*D.BA_AREA[]*D.f_bas[]
    D.basw[]            = D.basw_relative[] * pars.k_basw[]    

    D.silw[] = D.granw[] + D.basw[]
    k_silw = pars.k_granw[] + pars.k_basw[]
    D.silw_relative[] = D.silw[] / k_silw  # define relative rate for use by other modules

    # apatite associated with silicate weathering
    # does it follow its own kinetics? (or that of the host rock)
    if pars.f_p_kinetics[] == "no"
        D.granw_ap[] = D.granw[]
        D.basw_ap[] = D.basw[]
    elseif pars.f_p_kinetics[] == "yes"
        D.granw_ap[] = D.granw[] * (D.f_ap[] / D.f_gran[])
        D.basw_ap[] = D.basw[] * (D.f_ap[] / D.f_bas[])
    else
        error("unrecognized f_p_kinetics ", pars.f_p_kinetics[])
    end
    
    # is the P content assumed to vary between granite and basalt?
    if pars.f_p_apportion[] == "no"
        D.silw_ap[] = D.granw_ap[] + D.basw_ap[]
    elseif pars.f_p_apportion[] == "yes"
        D.silw_ap[] = (pars.k_P[]*D.granw_ap[] + D.basw_ap[])/(pars.k_P[]*(1-pars.k_basfrac[])+pars.k_basfrac[])
    else
        error("unrecognized f_p_apportion ", pars.f_p_apportion[])
    end
    
    D.phosw_s[] = D.F_EPSILON[]*pars.k10_phosw[]*pars.k_Psilw[]*(D.silw_ap[]/(k_silw + eps()))  # trap 0/0 if k_silw = 0
    

    # carbonate weathering 
    if pars.f_carb_link_u[] == "yes"
        carbw_fac_u = D.UPLIFT[]
    elseif pars.f_carb_link_u[] == "no"
        carbw_fac_u = 1.0
    else
        error("Unknown f_carb_link_u ", pars.f_carb_link_u[])
    end
    D.carbw_fac[] = carbw_fac_u * D.PG[]*D.RHO[]* D.CARB_AREA[] * D.f_carb[]
  
    if pars.f_carbwC[] == "Cindep"   # Copse 5_14
        carbw_fac_C = 1.0
    elseif pars.f_carbwC[] == "Cprop"    # A generalization for varying-size C reservoir
        carbw_fac_C = D.C_norm[]
    else
        error("Unknown f_carbw ", pars.f_carbwC[])
    end
    D.carbw_relative[] = carbw_fac_C * D.carbw_fac[]
    D.carbw[] = D.carbw_relative[] * pars.k14_carbw[] 
    
    # Define a normalized weathering rate for use by tracers etc.
    D.silwcarbw_relative[] = ( ( D.silw[] + D.carbw[] ) / ( k_silw + pars.k14_carbw[] ) )

    # Oxidative weathering

    # C oxidative weathering
    # not affected by Dforce.PG in G3 paper...
    if pars.f_oxwG[] == "Gindep"
        oxw_fac_G = 1.0
    elseif pars.f_oxwG[] == "Gprop"
        oxw_fac_G = D.G_norm[]    
    elseif pars.f_oxwG[] == "forced"
        oxw_fac_G = D.G_norm[]  * D.ORG_AREA[]
    else
        error("Unknown f_oxwG ", pars.f_oxwG[])
    end
    D.orgcw[] = oxw_fac_G * pars.k17_oxidw[]*D.UPLIFT[] * D.oxw_facOmax[]
    D.oxidw[] = oxw_fac_G * pars.k17_oxidw[]*D.UPLIFT[] * D.oxw_facO[] * D.koxidw[]

    # Sulphur weathering
    if pars.enableS[]
        # Gypsum weathering
        # not tied to Dforce.PG in G3 paper...
        if pars.f_gypweather[] == "original" # Gypsum weathering tied to carbonate weathering
            D.gypw[] = pars.k22_gypw[] * D.GYP_norm[]*D.carbw_fac[]
        elseif pars.f_gypweather[] == "alternative" # independent of carbonate area
            D.gypw[] = pars.k22_gypw[]*D.GYP_norm[]*D.UPLIFT[]*D.PG[]*D.RHO[]*D.f_carb[]
        elseif pars.f_gypweather[] == "forced" # dependent on evaporite area
            D.gypw[] = pars.k22_gypw[]*D.GYP_norm[]*D.EVAP_AREA[]*D.UPLIFT[]*D.PG[]*D.RHO[]*D.f_carb[]
        else
            error("unknown f_gypweather ", pars.f_gypweather[])
        end

        # Pyrite oxidative weathering 
        # not tied to Dforce.PG in G3 paper...
        if pars.f_pyrweather[] == "copse_O2"
            # with same functional form as carbon
            D.pyrw[] = pars.k21_pyrw[]*D.UPLIFT[]*D.PYR_norm[]*D.oxw_facO[]
        elseif pars.f_pyrweather[] == "copse_noO2"    # independent of O2
            D.pyrw[] = pars.k21_pyrw[]*D.UPLIFT[]*D.PYR_norm[]
        elseif pars.f_pyrweather[] == "forced" # forced by exposed shale area
            D.pyrw[] = pars.k21_pyrw[]*D.UPLIFT[]*D.SHALE_AREA[]*D.PYR_norm[]
        else
            error("unknown f_pyrweather ", pars.f_pyrweather[])
        end
        # Isotope fractionation of S
        gypw_isotope = @PB.isotope_totaldelta(SIsotopeType, D.gypw[], D.GYP_delta[])
        pyrw_isotope = @PB.isotope_totaldelta(SIsotopeType, D.pyrw[], D.PYR_delta[])
    else
        D.gypw[] = 0.0; gypw_isotope = @PB.isotope_totaldelta(SIsotopeType, 0.0, 0.0)
        D.pyrw[] = 0.0; pyrw_isotope = @PB.isotope_totaldelta(SIsotopeType, 0.0, 0.0)
    end

    # P weathering 
    # D.phosw_s is defined above

    # Introduction of P weathering flux from sandstones etc
    D.sedw_relative[] = D.UPLIFT[]*D.PG[]*D.RHO[]*D.VEG[]
    D.phosw_x[] = D.F_EPSILON[]*pars.k10_phosw[]*pars.k_Psedw[]*D.sedw_relative[]
    
    D.phosw_c[] = D.F_EPSILON[]*pars.k10_phosw[]*pars.k_Pcarbw[]*(D.carbw[]/(pars.k14_carbw[] + eps()))
    D.phosw_o[] = D.F_EPSILON[]*pars.k10_phosw[]*pars.k_Poxidw[]*(D.oxidw[]/(pars.k17_oxidw[] + eps()))  # trap 0/0 if k17_oxidw = 0
    D.phosw[]   = D.phosw_s[] + D.phosw_c[] + D.phosw_o[] + D.phosw_x[]


    #%%%%%%% Burial
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # Land organic carbon burial
    if     pars.f_locb[] == "original"
        D.pland[] = pars.k11_landfrac[]*D.VEG[]*D.phosw[]
    elseif pars.f_locb[] == "Uforced"
        # Uplift/erosion control of locb
        D.pland[] = pars.k11_landfrac[]*D.UPLIFT[]*D.VEG[]*D.phosw[]
    elseif pars.f_locb[] == "coal"
        # Coal basin forcing of locb
        D.pland[] = pars.k11_landfrac[]*D.VEG[]*D.phosw[]*(pars.k_aq[]+(1-pars.k_aq[])*D.COAL[])
    elseif pars.f_locb[] == "split"
        # Separating aquatic and coal basin components of locb
        D.pland[] = pars.k11_landfrac[]*D.VEG[]*D.phosw[]*(pars.k_aq[]*D.UPLIFT[]+(1-pars.k_aq[])*D.COAL[])
    elseif pars.f_locb[] == "Prescribed"
        # locb forced
        D.pland[] = pars.k11_landfrac[]*D.phosw[]
    else
        error("unknown f_locb ", pars.f_locb[])
    end

    D.psea[] = D.phosw[] - D.pland[]

    # Land organic carbon burial
    if pars.f_locb[] != "Prescribed"
        D.locb[] = D.pland[]*pars.CPland0[]*D.CPland_relative[]           
    end
   

    ## C isotopes
    ############################################################################
    # Isotope fractionation of organic carbon
    if CIsotopeType <: PB.AbstractIsotopeScalar
        D.locb_delta[]      = D.CO2_delta[] - D.D_P_CO2_locb[]        
    end
    locb_isotope       =  @PB.isotope_totaldelta(CIsotopeType, D.locb[], D.locb_delta[])

    oxidw_isotope       = @PB.isotope_totaldelta(CIsotopeType, D.oxidw[], D.G_delta[])

    carbw_isotope       = @PB.isotope_totaldelta(CIsotopeType, D.carbw[], D.C_delta[])

    # DIC isotopes - need to take account of atmosphere-water fractionation
    # (not critical to get this right, as there is a 'short circuit' atm <-> river -> ocean <-> atm)
    DICrunoff           = @PB.isotope_totaldelta(CIsotopeType, D.carbw[] + 2*D.silw[], D.CO2_delta[] + D.D_eqbw_CO2[])

    # fluxes
    #########################################################################

    # Atmospheric fluxes 
    fluxAtoLand.CO2[]   += DICrunoff - oxidw_isotope + locb_isotope
    
    fluxAtoLand.O2[]    += D.oxidw[] + 2*D.pyrw[] - D.locb[]

    # Riverine fluxes
    fluxRtoOcean.DIC[]  += DICrunoff + carbw_isotope   

    fluxRtoOcean.TAlk[] += 2*D.silw[] + 2*D.carbw[]
    if pars.SRedoxAlk[]
        fluxRtoOcean.TAlk[] += -2*D.pyrw[]
    end

    fluxRtoOcean.Ca[]   += D.silw[] + D.carbw[] + D.gypw[]

    fluxRtoOcean.P[]    += D.psea[]
    
    fluxRtoOcean.SO4[]  += gypw_isotope + pyrw_isotope

    # Sedimentary fluxes
    fluxLandtoSedCrust.Ccarb[]  += -carbw_isotope

    fluxLandtoSedCrust.Corg[]   += locb_isotope -oxidw_isotope
                
    fluxLandtoSedCrust.GYP[]    += -gypw_isotope

    fluxLandtoSedCrust.PYR[]    += -pyrw_isotope
                
    return nothing
end


end # module