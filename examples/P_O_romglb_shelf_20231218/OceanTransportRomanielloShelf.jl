module OceanTransportRomanielloShelf

import SparseArrays
import LinearAlgebra
import MAT   # Matlab file access

import PALEOboxes as PB
# import PALEOreactions
import PALEOocean

import Infiltrator # Julia debugger

using SparseArrays, LinearAlgebra

const BRAILLE = split("⠀⠁⠂⠄⡀⠈⠐⠠⢀", "") .|> s -> Int(s[1])

"""
    ReactionOceanTransportRomanielloShelf

Modern global ocean transport from [Romaniello2010a](@cite) [Romaniello2010](@cite) with 2-box shelves

Reads Matlab .mat files created from the published SI with Matlab commands:

    >> circ_global_79_box = romaniellog32010.Global_79_Box_ICBM_Params
    >> save('romaniello_global79','-struct', 'circ_global_79_box',  '-v6')
"""
Base.@kwdef mutable struct ReactionOceanTransportRomanielloShelf{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple([
        PB.ParString("matdir", "romaniello2010_transport",
            description="folder with Romaniello (2010) transport and geometry data files"),
        PB.ParStringVec("shelf_names", ["shelf"],
            description="names for additional shelf columns"),
        PB.ParDoubleVec("shelf_areas", [2e13], units="m-2",
            description="area of shelf columns"),
        PB.ParDoubleVec("shelf_depths", [150.0], units="m-2",
            description="area of shelf columns"),
        PB.ParDoubleVec("shelf_temps", [25.0], units="Centigrade",
            description="shelf temperatures"),
        PB.ParDouble("global_shelf_frac", 1e-2,
            description="fraction of shelf in global model to retain"),
        PB.ParDoubleVec("tshelfexch", [1.0], units="yr",
            description="vertical exchange time for 2-box shelf"),
        PB.ParDoubleVec("toceanexch", [1.0], units="yr",
            description="ocean exchange time for lower box of 2-box shelf"),
        PB.ParIntVec("ioceanshelf", [22],
            description="indices of ocean boxes to connect shelf lower boxes to"),
        PB.ParDoubleVec("toceanexch2", Float64[], units="yr",
            description="optional ocean exchange time for lower box of 2-box shelf, 2nd set of ocean boxes"),
        PB.ParIntVec("ioceanshelf2", Int[],
            description="optional indices of 2nd set of ocean boxes to connect shelf lower boxes to"),
        PB.ParDoubleVec("tansportSv", [18.0], units="Sv (10^6 m3/s)",
            description="flux transfer from gyre boxes to bottom shelf"),

        PB.ParBool("remove_upw", false,
            description="true to remove the upw column"),
        PB.ParBool("remove_shelf", false,
            description="true to remove the shelf area from hlat and gyre columns"),
        PB.ParString("volume_method", "Read",
            description="Calculation method of volume, options (1) Read, read from romglb.mat file, 
                        (2) Trapezoid, (3) Rectangle"),
        PB.ParBool("Stuart_adv", true,
            description="True to set Stuart's advection flux (none upwelling)"),
        PB.ParBool("temp_trackglobal", false,
            description="track global temperature (apply offset of 'global temp - 15C')"),
        PB.ParBool("domain_oceanhlatgyresurface", false,
            description="true to configure a domain that just includes hlat and gyre surface boxes"),
    ])

    temp_oceanC::Vector{Float64} = Float64[] # ocean temperature (degrees C)

    oceanbase           = nothing

    grid_ocean          = nothing
    grid_oceansurface   = nothing
    grid_oceanfloor     = nothing
    grid_oceanhlatgyresurface = nothing # optional
    

    "Transport matrix: Units yr^{-1}
     so dc/dt = trspt_dtm * c [yr^-1]
     where c is column vector of tracer concentrations"
    trspt_dtm::SparseArrays.SparseMatrixCSC{Float64,Int64} = SparseArrays.spzeros(0, 0)

    "Transpose of trspt_dtm"
    trspt_dtm_tr::SparseArrays.SparseMatrixCSC{Float64,Int64} = SparseArrays.spzeros(0, 0)

    rom_data                    = nothing  # raw data from .mat file
end
   


function read_datafiles(rj::ReactionOceanTransportRomanielloShelf)

    matdir = rj.pars.matdir[] # folder containing mat files

    # rom keys:
    # rom[<colname>]                        indices of boxes for this column, ordered surface to floor
    # rom["depths"]                         mid depths (m) cells
    # rom[<colname>_bnd]                    indices of horizontal surfaces for this column (length of column + 1)
    # rom["bnd_depths"][rom[<colname>_bnd]] depths (m) of horizontal surfaces (length of column + 1)
    # rom["Hyps"][rom[<colname>_bnd]]       area (m^2) of horizontal surfaces (length of column + 1)
    # rom["T"]                              transport matrix

    
    matfilename = joinpath(matdir, "romaniello_global79.mat")
    @info "$(PB.fullname(rj)) reading transport matrix from file $(matfilename)"
    rj.rom_data = MAT.matread(matfilename)
    
    if rj.pars.remove_upw.v
        columnnames=[:hlat, :gyre]    
        Icolumns=[Int.(vec(rj.rom_data["hlat"])), Int.(vec(rj.rom_data["gyre"]))]
        ocean_cells=Int(rj.rom_data["nboxes"]) - 33 # 79 = 13 + 33 + 33
    else
        columnnames=[:hlat, :gyre, :upw]    
        Icolumns=[Int.(vec(rj.rom_data["hlat"])), Int.(vec(rj.rom_data["gyre"])), Int.(vec(rj.rom_data["upw"]))]
        ocean_cells=Int(rj.rom_data["nboxes"])   
    end
      
    # adding shelves
    for shelf_name in rj.pars.shelf_names.v
        push!(columnnames, Symbol(shelf_name))     
        push!(Icolumns, [ocean_cells+1, ocean_cells+2])
        ocean_cells += 2
    end

    rj.grid_ocean = PB.Grids.UnstructuredColumnGrid(
        # domain=rj.domain,
        ncells=ocean_cells,
        columnnames=columnnames,
        Icolumns=Icolumns
    )
    
    return nothing
end


function PB.set_model_geometry(rj::ReactionOceanTransportRomanielloShelf, model::PB.Model)

    read_datafiles(rj)

    ocean_cells = rj.grid_ocean.ncells # total cells including shelf columns

    # set subdomain mappings
    isurf=[col[1] for col in rj.grid_ocean.Icolumns]    # first index in each column is surface box    
    ifloor=collect(1:ocean_cells)                       # every box has an oceanfloor box under it
    
    PB.Grids.set_subdomain!(rj.grid_ocean, "oceansurface", PB.Grids.BoundarySubdomain(isurf), true)
    @info "  set ocean.oceansurface Subdomain size=$(length(isurf))"
    PB.Grids.set_subdomain!(rj.grid_ocean, "oceanfloor", PB.Grids.BoundarySubdomain(ifloor), true)
    @info "  set ocean.oceanfloor Subdomain size=$(length(ifloor))"

    rj.grid_oceansurface = PB.Grids.UnstructuredVectorGrid(
        # domain=PB.get_domain(model, "oceansurface"),
        ncells=length(isurf),
        cellnames=Dict(rj.grid_ocean.columnnames[i]=>i for i in 1:length(rj.grid_ocean.columnnames)),
    )
    PB.Grids.set_subdomain!(rj.grid_oceansurface, "ocean", PB.Grids.InteriorSubdomain(ocean_cells, isurf), true)

    # every box has an oceanfloor box under it    
    rj.grid_oceanfloor = PB.Grids.UnstructuredColumnGrid(
        # domain=PB.get_domain(model, "oceanfloor"),
        ncells=length(ifloor),
        columnnames=rj.grid_ocean.columnnames,
        Icolumns=rj.grid_ocean.Icolumns,
    )
    PB.Grids.set_subdomain!(rj.grid_oceanfloor, "ocean", PB.Grids.InteriorSubdomain(ocean_cells, ifloor), true)

    PALEOocean.Ocean.set_model_domains(model, rj.grid_ocean, rj.grid_oceansurface, rj.grid_oceanfloor)    
       
    # optionally define a domain that just includes the  gyre and hlat surface boxes (for nutrient etc input)
    if rj.pars.domain_oceanhlatgyresurface[]
        domname = "oceanhlatgyresurface"
        # columns 1 and 2 are hlat, gyre
        icol_hlatgyre = [1, 2]
 
        isurf=[col[1] for col in rj.grid_ocean.Icolumns[icol_hlatgyre]]    # first index in each column is surface box        
        PB.Grids.set_subdomain!(rj.grid_ocean, domname, PB.Grids.BoundarySubdomain(isurf), true)
        @info "  set ocean.$domname Subdomain size=$(length(isurf))"
        rj.grid_oceanhlatgyresurface = PB.Grids.UnstructuredVectorGrid(
            # domain=PB.get_domain(model, "oceansurface"),
            ncells=length(isurf),
            cellnames=Dict(rj.grid_ocean.columnnames[i]=>i for i in icol_hlatgyre),
        )
        PB.Grids.set_subdomain!(rj.grid_oceanhlatgyresurface, "ocean", PB.Grids.InteriorSubdomain(ocean_cells, isurf), true)
 
        dom = PB.get_domain(model, "oceanhlatgyresurface")
        if !isnothing(dom)
            dom.grid = rj.grid_oceanhlatgyresurface
            @info "  $(domname) Domain size=$(dom.grid.ncells) grid=$(dom.grid)"
        end
    end

    return nothing
end

function PB.register_methods!(rj::ReactionOceanTransportRomanielloShelf)
   
    physvars = [
        PB.VarProp("sal", "psu", "Ocean salinity"),
        PB.VarProp("rho", "kg m^-3", "physical ocean density"),
        PB.VarProp("oceansurface.open_area_fraction", "", "fraction of area open to atmosphere", attributes=(:check_length=>false,)),
        PB.VarPropScalar("number_cells","", "number of ocean cells")
    ]

    PB.add_method_setup!(
        rj, 
        do_setup_grid,
        (   
            PB.VarList_namedtuple(PALEOocean.Ocean.grid_vars_all), 
            PB.VarList_namedtuple(physvars),
        ),
    )

    tempvars = [
        PB.VarDepScalar("(global.TEMP)",            "K",        "global mean temperature"),
        PB.VarProp("temp",                          "Kelvin",   "Ocean temperature"),
    ]

    PB.add_method_do!(
        rj, 
        do_temperature,
        (PB.VarList_namedtuple(tempvars), ),
    )    

    return nothing
end


function do_setup_grid(
    m::PB.ReactionMethod,
    (grid_vars, physvars),
    cellrange::PB.AbstractCellRange,
    attribute_name,
)
    rj = m.reaction

    attribute_name == :setup || return nothing

    # rom keys:
    # rom[<colname>]                        indices of boxes for this column, ordered surface to floor
    # rom["depths"]                         mid depths (m) cells
    # rom[<colname>_bnd]                    indices of horizontal surfaces for this column (length of column + 1)
    # rom["bnd_depths"][rom[<colname>_bnd]] depths (m) of horizontal surfaces (length of column + 1)
    # rom["Hyps"][rom[<colname>_bnd]]       area (m^2) of horizontal surfaces (length of column + 1)
    # rom["T"]                              transport matrix


    ocean_cells = rj.grid_ocean.ncells              # total cells including shelf columns
    if rj.pars.remove_upw.v
        num_rom_columns = 2
        rom_cells = Int(rj.rom_data["nboxes"]) - 33 # cells in 2 column transport matrix (excluding upw and shelves)
    else
        num_rom_columns = 3
        rom_cells = Int(rj.rom_data["nboxes"])      # cells in 3 column transport matrix (excluding shelves)
    end

    for icl in 1:num_rom_columns
        iboxes = rj.grid_ocean.Icolumns[icl]        # indices for this column
        ibnd = Int.(rj.rom_data[String(rj.grid_ocean.columnnames[icl])*"_bnd"]) # ibnd is a matrix

        grid_vars.zupper[iboxes]  .= -rj.rom_data["bnd_depths"][ibnd[1:end-1]]
        grid_vars.zmid[iboxes] .= -rj.rom_data["depths"][iboxes]
        # grid_vars.zlower[iboxes]  .= grid_vars.zupper[iboxes]-grid_vars.volume[iboxes]./grid_vars.Abox[iboxes]
        grid_vars.zlower[iboxes]  .= -rj.rom_data["bnd_depths"][ibnd[2:end]]
        grid_vars.zfloor[iboxes]  .= grid_vars.zlower[iboxes]    
        
        hyps = rj.rom_data["Hyps"][ibnd[1:end-1]]
        # if icl in 1:2
        #     # remove shelf (NB: this will generate an inconsistency between unmodified volume and modified area)
        #     hyps[1] = hyps[2] + (hyps[1]-hyps[2])*rj.pars.global_shelf_frac.v
        # end
        if (rj.pars.remove_shelf.v) && (icl in 1:2)
            # remove shelf area
            indices = findall(x -> x > -200, grid_vars.zmid[iboxes])   # find the index with water depth < 200m
            @info "Test, indices=$(indices), depth=$(grid_vars.zmid[iboxes]), hyps=$(hyps)"
            # all the boxes over 200m have the same area as 200m_box
            # huge issue here if the area of two vertical_linked boxes are the same, NaN would appear when compling Dunne2007
            hyps[indices] .= hyps[maximum(indices)+1] .+ 0.1 .* collect(length(indices):-1:1) 
            @info "Test_results, indices=$(indices), depth=$(grid_vars.zmid[iboxes]), hyps=$(hyps)"
        end
        hyps_bottom = vcat(hyps[2:end], rj.rom_data["Hyps"][ibnd[end]])
        grid_vars.Abox[iboxes]  .= hyps
        grid_vars.Asurf[icl]  = hyps[1]
        # assume column closed at bottom box (Romaniello hypsometry has small term here)
        grid_vars.Afloor[iboxes[1:end-1]]  .= grid_vars.Abox[iboxes[1:end-1]] - grid_vars.Abox[iboxes[2:end]]
        grid_vars.Afloor[iboxes[end]] = grid_vars.Abox[iboxes[end]]

        # NB: there is an inconsistency in Matlab (rj.rom_data) vs PALEO conventions: volume != Abox * (zupper - zlower)
        if rj.pars.volume_method.v == "Read"
            grid_vars.volume[iboxes] .= rj.rom_data["V"][iboxes]/1000.0 # convert to m^3
        elseif rj.pars.volume_method.v == "Rectangle"
            grid_vars.volume[iboxes] .=  grid_vars.Abox[iboxes] .* (grid_vars.zupper[iboxes] .- grid_vars.zlower[iboxes])
        elseif rj.pars.volume_method.v == "Trapezoid"
            grid_vars.volume[iboxes] .=  (grid_vars.Abox[iboxes].+hyps_bottom)/2 .* (grid_vars.zupper[iboxes] .- grid_vars.zlower[iboxes])
        else
            @error "No method $(rj.pars.volume_method.v) matches volume_method in OceanTransportRomanielloShelf"
        end
    end   

    # add shelf columns    
    upper_floor_frac = 1e-3  # add a neglible floor area for upper box
    for ishelf in 1:length(rj.pars.shelf_names.v)
        icl = ishelf + num_rom_columns
        iboxes = rj.grid_ocean.Icolumns[icl]  # indices for this column

        grid_vars.Abox[iboxes] .= [1.0, 1.0-upper_floor_frac].*rj.pars.shelf_areas.v[ishelf]
        grid_vars.Asurf[icl] = rj.pars.shelf_areas.v[ishelf]       
        grid_vars.Afloor[iboxes] .= [upper_floor_frac, 1.0-upper_floor_frac].*rj.pars.shelf_areas.v[ishelf]

        grid_vars.zupper[iboxes] .= [0.0, -0.5].*rj.pars.shelf_depths.v[ishelf]
        grid_vars.zmid[iboxes]   .= [-0.25, -0.75].*rj.pars.shelf_depths.v[ishelf]
        grid_vars.zlower[iboxes] .= [-0.5, -1.0].*rj.pars.shelf_depths.v[ishelf]
        grid_vars.volume[iboxes] .= grid_vars.Abox[iboxes] .* (grid_vars.zupper[iboxes] .- grid_vars.zlower[iboxes])

        grid_vars.zfloor[iboxes] .= grid_vars.zlower[iboxes]
    end


    grid_vars.volume_total[] = sum(grid_vars.volume)
    grid_vars.Afloor_total[] = sum(grid_vars.Afloor)
    
    # attach coordinates to grid for output visualisation etc
    empty!(rj.domain.grid.z_coords)
    push!(rj.domain.grid.z_coords, PB.FixedCoord("zmid", grid_vars.zmid, PB.get_variable(m, "zmid").attributes))
    push!(rj.domain.grid.z_coords, PB.FixedCoord("zlower", grid_vars.zlower, PB.get_variable(m, "zlower").attributes))
    push!(rj.domain.grid.z_coords, PB.FixedCoord("zupper", grid_vars.zupper, PB.get_variable(m, "zupper").attributes))

    # constant density
    grid_vars.rho_ref       .= 1027
    
    # set default ocean temperature 
    rj.temp_oceanC = Vector{Float64}(undef, ocean_cells)
    rj.temp_oceanC .= 2.0 # all interior boxes at 2C
    if rj.pars.remove_upw.v
        rj.temp_oceanC[rj.domain.grid.subdomains["oceansurface"].indices[1:num_rom_columns]] .= [2.5, 25.0] # guessed
    else
        rj.temp_oceanC[rj.domain.grid.subdomains["oceansurface"].indices[1:num_rom_columns]] .= [2.5, 25.0, 25.0] # guess
    end

    # add shelf columns
    for ishelf in 1:length(rj.pars.shelf_names.v)
        icl = ishelf + num_rom_columns
        iboxes = rj.grid_ocean.Icolumns[icl]  # indices for this column        
        rj.temp_oceanC[iboxes] .= rj.pars.shelf_temps.v[ishelf]
    end

    grid_vars.pressure  .= -grid_vars.zmid   # pressure(dbar) ~ depth (m)
    
    # generate Transport Matrix
    Conv = 3.15576e16   # Conversion factor between Sverdrups and kg/yr assuming 1 m3 = 1000 kg (1Sv = 10^6 m3/s)

    rj.trspt_dtm = SparseArrays.spzeros(ocean_cells, ocean_cells)

    # # Ziheng's adding shelf, with upwelling
    # (rj.trspt_dtm, _) = calculate_rom_transport_matrix(rj)

    # Stuart's adding shelf columns only, no upwelling flux
    if rj.pars.remove_upw.v
        (_, rj.trspt_dtm[1:rom_cells, 1:rom_cells]) = calculate_rom_transport_matrix(rj, grid_vars)
    else
        rj.trspt_dtm[1:rom_cells, 1:rom_cells] = calculate_rom_transport_matrix(rj, grid_vars)
    end
    
    if rj.pars.Stuart_adv.v
        for ishelf in 1:length(rj.pars.shelf_names.v)
            iu, il = rj.grid_ocean.Icolumns[num_rom_columns+ishelf] # indices of shelf upper, lower boxes
            # shelf upper-low exchange flux 
            tshelfexch = rj.pars.tshelfexch.v[ishelf]
            rj.trspt_dtm[iu, iu] += -1.0/tshelfexch
            rj.trspt_dtm[il, il] += -1.0/tshelfexch
            rj.trspt_dtm[iu, il] += 1.0/tshelfexch
            rj.trspt_dtm[il, iu] += 1.0/tshelfexch

            # shelf lower to ocean exchange flux
            ie = rj.pars.ioceanshelf.v[ishelf] # index of ocean box to exchange with
            toceanexch = rj.pars.toceanexch.v[ishelf]
            vl = grid_vars.volume[il]
            ve = grid_vars.volume[ie]        
            
            rj.trspt_dtm[il, il] += -1.0/toceanexch
            rj.trspt_dtm[ie, ie] += -(vl/ve)/toceanexch
            rj.trspt_dtm[il, ie] += 1.0/toceanexch
            rj.trspt_dtm[ie, il] += (vl/ve)/toceanexch

            # optional 2nd set of boxes for shelf lower to ocean exchange flux
            if !isempty(rj.pars.ioceanshelf2.v)
                ie = rj.pars.ioceanshelf2.v[ishelf] # index of ocean box to exchange with
                toceanexch = rj.pars.toceanexch2.v[ishelf]
                vl = grid_vars.volume[il]
                ve = grid_vars.volume[ie]        
                
                rj.trspt_dtm[il, il] += -1.0/toceanexch
                rj.trspt_dtm[ie, ie] += -(vl/ve)/toceanexch
                rj.trspt_dtm[il, ie] += 1.0/toceanexch
                rj.trspt_dtm[ie, il] += (vl/ve)/toceanexch
            end
        end
    end

    @info "The Transport Matrix of Gyre Hlat and shelves:"
    show_any_nonzero(rj.trspt_dtm, maxw=24)
    
    rj.trspt_dtm_tr = SparseArrays.sparse(transpose(rj.trspt_dtm))

    # set salinity
    physvars.sal                .= 35.0
    # constant density
    physvars.rho                .= 1027
    
    physvars.open_area_fraction .= 1.0

    num_cells = size(rj.trspt_dtm_tr, 1)
    @info "number of cells is $(num_cells)"
    physvars.number_cells       .= float(num_cells)

    return nothing
end     
           

function do_temperature(m::PB.ReactionMethod, (tempvars, ), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction

    # Set temperature
    if rj.pars.temp_trackglobal.v
        tempvars.temp      .= rj.temp_oceanC .- 15.0 .+ tempvars.TEMP[] # temperature (K)
    else
        tempvars.temp      .= rj.temp_oceanC .+ PB.Constants.k_CtoK # temperature (K)
    end

    return nothing
end

function PB.register_dynamic_methods!(rj::ReactionOceanTransportRomanielloShelf)

    (transport_conc_vars, transport_sms_vars, transport_input_vars) =
        PALEOocean.Ocean.find_transport_vars(rj.domain, add_transport_input_vars=true)

    PB.add_method_do!(
        rj, 
        do_transport,
        (   
            PB.VarList_namedtuple(PB.VarDep.(PALEOocean.Ocean.grid_vars_ocean)),
            PB.VarList_components(transport_conc_vars),
            PB.VarList_components(transport_sms_vars),
            PB.VarList_components(transport_input_vars),
        ),
        preparefn=PALEOocean.Ocean.prepare_transport
    )

    return nothing
end


function do_transport(
    m::PB.ReactionMethod,
    (grid_vars, transport_conc_components, transport_sms_components, transport_input_components, buffer), 
    cellrange::PB.AbstractCellRange, 
    deltat
)
    rj = m.reaction

    PALEOocean.Ocean.do_transport_tr(
        grid_vars, transport_conc_components, transport_sms_components, transport_input_components, buffer,
        rj.trspt_dtm_tr, 
        cellrange
    )

    return nothing
end

"""
sub2ind

convert the 2D or nD subscripts to linear indeices
"""
function sub2ind(
    range,  # this is a e.g. range[1] * range[2] matrix
    subs,
    )
    k1 = subs[:,1]
    k2 = subs[:,2]
    ind = (k2 .- 1) * range[1] .+ k1
    return ind
end

"""
norm_exp

Forms an normalized exponential function and evaluates it at x_int

    upperbnd = scalar value where function is maximum
    lowerbnd = scalar value where function set equal to zero
    L = length scale of exponent
    x = values where function is returned

    f(x) = A*(exp(-x/L)-exp(-lowerbnd/L))
    f(upperbnd) == 1
    returns values = f(x)
"""
function norm_exp(upperbnd,lowerbnd,L,x)
    A = 1/(exp(-upperbnd/L)-exp(-lowerbnd/L))
    values = A*(exp(-x/L)-exp(-lowerbnd/L))
    return values
end

"""
remove_shelf_area!
"""
function remove_shelf_area!(Hyps)
    for i = 1:6 # 0 to ~200m for hlat column, bnd_depth[7] = 208m
        Hyps[i] = Hyps[7]
    end
    for i = 15:17 # 0 to ~200m for gyre column, bnd_depth[18] = 200m
        Hyps[i] = Hyps[18]
    end
end

"""
calculate_rom_transport_matrix

Two options: keep or remove the upw column, return a Sparse Matrix
Two options: keep or remove the shelf areas from gyre and hlat columns, return a Sparse Matrix

Update: do not read ocean's physical pars from romglb's .mat file, instead get the data from argument `grid_vars`
"""
function calculate_rom_transport_matrix(rj::ReactionOceanTransportRomanielloShelf, grid_vars)

    remove_upw = rj.pars.remove_upw.v
    remove_shelf = rj.pars.remove_shelf.v

    Ocean_SA = 4.2953e14    # m^2
    Ocean_V = 135e16        # m^3
    SA_fractions = [.0240,.8569,.1191]

    # read the model input from Steve's src
    circname = "Global_79_Box" # "Global_79_Box", "Black_Sea"
    # matdir = joinpath(PALEOreactions.srcdir(), "ocean")
    matdir = @__DIR__
    matfilename = joinpath(matdir, "romaniello_global79.mat")
    rom_data = MAT.matread(matfilename)

    # What we need are only
    # 1. boxes framework box_thickness and indices
    # 2. volume of each box, pp, piecewise polynomial, satellite-based ETOPO2 (smith & sandwell 2001), another important to explore

    hlat = Int.(vec(rom_data["hlat"]))                # vector e.g. 1:3 or 1:13, indices of boxes for this column, ordered surface to floor
    gyre = Int.(vec(rom_data["gyre"]))
    upw = Int.(vec(rom_data["upw"]))
    # rj.grid_ocean.Icolumns equal to the vec(hlat, gyre, upw)

    hlat_bnd = [vec(hlat); hlat[end]+1]             # boundary indices for each boxes (length of column + 1)
    gyre_bnd = [vec(gyre); gyre[end]+1].+1 
    upw_bnd = [vec(upw); upw[end]+1].+2  

    box_thickness = grid_vars.zupper .- grid_vars.zlower # vertical thickness of each box (m)

    if remove_upw
        bnd_depths = [0; cumsum(box_thickness[hlat]);   # depths (m) of horizontal surfaces (boundaries) (length of column + 1)
                0; cumsum(box_thickness[gyre])]
        
        depths = [(bnd_depths[hlat_bnd[1:(length(hlat_bnd)-1)]] + bnd_depths[hlat_bnd[2:length(hlat_bnd)]])/2;
                (bnd_depths[gyre_bnd[1:(length(gyre_bnd)-1)]] + bnd_depths[gyre_bnd[2:length(gyre_bnd)]])/2]
                
        Hyps = rom_data["Hyps"][1:48]                      # surface area (m^2) of horizontal surfaces (length of column + 1)
        ocean_cells = Int(rom_data["nboxes"]) - 33
        V = grid_vars.volume[1:ocean_cells] * 1000     # box volume (kg, 1m^3 = 1000 kg)
    else
        bnd_depths = [0; cumsum(box_thickness[hlat]);   # depths (m) of horizontal surfaces (boundaries) (length of column + 1)
                0; cumsum(box_thickness[gyre]);
                0; cumsum(box_thickness[upw])]
        
        depths = [(bnd_depths[hlat_bnd[1:(length(hlat_bnd)-1)]] + bnd_depths[hlat_bnd[2:length(hlat_bnd)]])/2;
                (bnd_depths[gyre_bnd[1:(length(gyre_bnd)-1)]] + bnd_depths[gyre_bnd[2:length(gyre_bnd)]])/2;
                (bnd_depths[upw_bnd[1:(length(upw_bnd)-1)]] + bnd_depths[upw_bnd[2:length(upw_bnd)]])/2;]

        Hyps = rom_data["Hyps"]                      # surface area (m^2) of horizontal surfaces (length of column + 1)
        ocean_cells = Int(rom_data["nboxes"]) 
        V = grid_vars.volume[1:ocean_cells] * 1000     # box volume (kg, 1m^3 = 1000 kg)
    end   
    
    if remove_shelf
        # remove areas belonging to the shelf from gyre and hlat columns
        remove_shelf_area!(Hyps)
    end
    ################################################## 

    # specify advective fluxes
    Conv = 3.15576e16   # Conversion factor between Sverdrups and kg/yr assuming 1 m3 = 1000 kg (1Sv = 10^6 m3/s)
    # Conv = 1e6*1e3*365.25*24*60*60
    AABW = 17           # bottom water (Sv)
    AAIW = 12           # Intermediate water (Sv)
    NADW = 17           # Deep water (Sv), NADW > AAIW
    CstUpw = 18         # Upwelling in regions of intense coastal upwelling (Sv)
    Hlat_Overturn = 0   # Convection in High Latitudes (Sv)

    adv_mixing = spzeros(ocean_cells, ocean_cells)
    # adv_mixing = Antarctic_mixing + Arctic_mixing + upw_mixing
    # Steve calculate the vertial transport together, so donot seperate adv into dividied three

    if circname == "Global_13_Box"
        index_hlat_bottom_up = 3
        index_inter_up = 3 
        index_inter_low = 3
        index_bottom_up = 5
        index_upw_bottom = 2
    elseif circname == "Global_79_Box"
        index_hlat_bottom_up = 12 # hlat's ~1000m, the uppere bracket of the bottom
        index_inter_up = 15 # gyre's ~700-800m, the upper bracket of the intermediate boxes
        index_inter_low = 21 # gyre's ~1200-1300m, the lower bracket of the intermediate boxes
        index_bottom_up = 32 # gyre's ~5000m, the upper bracket of the bottom boxes, the lower bracket of the bottom boxes is gyre[end]
        index_upw_bottom = 19 # upw's ~1000m, the lower bracket of the upwelling zone
    end

    # Antarctic_mixing
    # do the deepwater advective flux
    # hlat[1] -> hlat[end]   -> gyre[end]
    #                             | 
    # hlat[1] <- hlat[end-1] <- gyre[end-1]
    adv_mixing[hlat[end],hlat[1]] = AABW
    adv_mixing[gyre[end],hlat[end]] = AABW
    adv_mixing[hlat[end-1],gyre[index_bottom_up]] = AABW
    # do intermediate water advective flux
    # Antarctic hlat[1] -> gyre[index_inter_up:index_inter_low]      -> Arctic hlat[1]
    # Arctic hlat[1]    -> gyre[index_inter_low+1:index_bottom_up-1] -> Antarctic hlat[2:end]
    adv_mixing[gyre[index_inter_up:index_inter_low],hlat[1]] = AAIW .* (box_thickness[gyre[index_inter_up:index_inter_low]] / sum(box_thickness[gyre[index_inter_up:index_inter_low]]))
    i = sub2ind([ocean_cells,ocean_cells],[hlat[2:(index_hlat_bottom_up-1)] gyre[(index_inter_low+1):(index_bottom_up-1)]])
    adv_mixing[i] = AAIW .* (box_thickness[gyre[(index_inter_low+1):(index_bottom_up-1)]] / sum(box_thickness[gyre[(index_inter_low+1):(index_bottom_up-1)]]))

    # Arctic_mixing
    # Antarctic hlat[1] -> gyre[index_inter_up:index_inter_low]      -> Arctic hlat[1]
    # Arctic hlat[1]    -> gyre[index_inter_low+1:index_bottom_up-1] -> Antarctic hlat[2:end]
    adv_mixing[hlat[1],gyre[index_inter_up:index_inter_low]] += NADW .* (box_thickness[gyre[index_inter_up:index_inter_low]] / sum(box_thickness[gyre[index_inter_up:index_inter_low]]))
    adv_mixing[gyre[index_inter_low+1:index_bottom_up-1],hlat[1]] += NADW .* (box_thickness[gyre[index_inter_low+1:index_bottom_up-1]] / sum(box_thickness[gyre[index_inter_low+1:index_bottom_up-1]]))

    if remove_upw
        # do nothing
    else
        # upw_mixing
        adv_mixing[gyre[1],upw[1]] += CstUpw * sum(-diff(norm_exp.(50,1000,100,bnd_depths[upw_bnd[2:(index_upw_bottom+1)]])))
        # adv_mixing[upw[6],gyre[6]] += CstUpw # not used, this is ~300m as the bottom of the upw
        i = sub2ind([ocean_cells,ocean_cells],[upw[2:index_upw_bottom] gyre[2:index_upw_bottom]])
        adv_mixing[i] += CstUpw .* (-diff(norm_exp.(50,1000,100,bnd_depths[upw_bnd[2:(index_upw_bottom+1)]])))
    end

    for i in 1:ocean_cells
        adv_mixing[i,i] -= sum(adv_mixing[i,:])
    end

    # Ziheng's adding shelf columns
    if remove_upw
        V_shelf = deepcopy(V)
        upw_mixing = SparseArrays.spzeros(rj.grid_ocean.ncells , rj.grid_ocean.ncells ) # total cells including shelf columns

        # upwelling fluxes from gyre boxes to bottom shelf
        for ishelf in 1:length(rj.pars.shelf_names.v)

            # prepare to convert Sv to yr-1, do it twice for upper and lower part of shelf
            append!(V_shelf, 0.5 * 1000 * rj.pars.shelf_areas.v[ishelf]*rj.pars.shelf_depths.v[ishelf])
            append!(V_shelf, 0.5 * 1000 * rj.pars.shelf_areas.v[ishelf]*rj.pars.shelf_depths.v[ishelf])

            iu, il  = rj.grid_ocean.Icolumns[2+ishelf]    # indices of shelf upper, lower boxes
            ie      = rj.pars.ioceanshelf.v[ishelf]                     # index of ocean box to exchange with

            # shelf1[il] <- ie, ie = gyre[?] 
            upw_mixing[il, ie] += rj.pars.tansportSv.v[ishelf]
            # shelf1[iu] <- shelf1[il]
            upw_mixing[iu, il] += rj.pars.tansportSv.v[ishelf]          # do vertical for shelves
            # gyre[1] <- shelf1[iu] 
            upw_mixing[gyre[1], iu] += rj.pars.tansportSv.v[ishelf]
        end

        for i in 1:rj.grid_ocean.ncells
            upw_mixing[i,i] -= sum(upw_mixing[i,:])
        end

        # do the vertial transports of gyre
        for i in length(gyre):-1:2
            convergence = sum(upw_mixing[:,gyre[i]])
            if convergence < 0
                upw_mixing[gyre[i-1], gyre[i]] -= convergence
                upw_mixing[gyre[i-1], gyre[i-1]] += convergence
            elseif convergence > 0
                upw_mixing[gyre[i], gyre[i-1]] += convergence
                upw_mixing[gyre[i], gyre[i]] -= convergence
            end
        end
    end
    ##################################################

    # do the vertial transports of adv_mixing
    if remove_upw
        vector_range = [hlat, gyre]
    else
        vector_range = [hlat, gyre, upw]
    end

    # from the bottom of each column, calculate upwards, 
    # if the net output flux of the lower box >0 then set a input to the bottom box from its overlying box
    # if the net output flux of the lower box <0 then set a output to its overlying box
    for range in vector_range
        if length(range) > 1
            for i in length(range):-1:2
                convergence = sum(adv_mixing[:,range[i]])
                if convergence < 0
                    adv_mixing[range[i-1], range[i]] -= convergence
                    adv_mixing[range[i-1], range[i-1]] += convergence
                elseif convergence > 0
                    adv_mixing[range[i], range[i-1]] += convergence
                    adv_mixing[range[i], range[i]] -= convergence
                end
            end
        else
            @error "hlat, gyre or upw's length error"
        end
    end
    ################################################## 

    # Nonlocal mixing
    # 13 box gyre[2]            <- (gyre[1]        mix      hlat[1])
    # 79 box gyre[2:19]         <- (gyre[1]        mix      hlat[1])
    #        vantilated box     equatorial_endmember box    polar_endmember box
    # Total ventilation flux = 220 Sv, gyre's 0-1000m can be ventilated directly from the gyre[1] and hlat[1]
    # gyre's 1000-6000m should be ventilated from hlat[1:end], represent as horizontal diffusion
    NonLocal_mixing = spzeros(ocean_cells, ocean_cells)
    Subtrop_Overturn = 220 # hallow Meridional Overturning Associated with Equatorial Upwelling and Subtropical Subduction (Sv)

    if circname == "Global_13_Box"
        num_box_surf_gyre = 1
        Kh_Outcrop = [220]
        T = [11] # Mean winter Outcrop Temperature (˚C), set to let feq = 0.57

    elseif circname == "Global_79_Box"
        num_box_surf_gyre = 18 # the number or indeices of the vantilated boxes
        # gyre[num_box_surf_gyre+1] should have a boundary depth of ~1000m

        # sum Kh_Outcrop = 220 (Sv) equal to the Subtrop_Overturn
        Kh_Outcrop = [40;40;40;20;15;15;10;10;6;4;4;4;4;4;1;1;1;1]
        # Mean winter Outcrop Temperature (˚C)
        T = [14.8398, 12.5611, 10.7175, 10.3659, 9.5439, 8.3651, 7.3111, 6.1024, 5.2852, 4.4662, 3.5275, 2.4268, 1.7761, 0.9855, 0.9403, 0.5058, 0.1729, -0.0911]
        # length(T) = length(Kh_Outcrop) = num_box_surf_gyre
    end

    feq = (T.-(-1))./(20-(-1)) # hlat[1] vs. gyre[1] contribution (%), calculate based only on temperature

    for i in 1:num_box_surf_gyre
        NonLocal_mixing[gyre[i+1],gyre[1]] += feq[i] * Kh_Outcrop[i]
        NonLocal_mixing[gyre[1],gyre[i+1]] += feq[i] * Kh_Outcrop[i]

        NonLocal_mixing[gyre[i+1],hlat[1]] += (1-feq[i]) * Kh_Outcrop[i]
        NonLocal_mixing[hlat[1],gyre[i+1]] += (1-feq[i]) * Kh_Outcrop[i]

        NonLocal_mixing[gyre[1],gyre[1]] -= feq[i] * Kh_Outcrop[i]
        NonLocal_mixing[hlat[1],hlat[1]] -= (1-feq[i]) * Kh_Outcrop[i]
        NonLocal_mixing[gyre[i+1],gyre[i+1]] -= Kh_Outcrop[i]
    end

    # All the rows and columns of T should always sum to zero to within numerical round off error, 
        # so long as the model has no allochthonus inputs and outflows of water (e.g. rivers, evaporation)
    for i=1:ocean_cells
        if !(sum(NonLocal_mixing[i,:]) < 1e-6 && sum(NonLocal_mixing[:,i]) < 1e-6)
            @error "Not equal to zero error in the advective circulation field"
        end
    end
    ################################################# 

    # Vertical diffusion
    Z_mixing = spzeros(ocean_cells, ocean_cells)
    # eddy diffusion coefficients
    Kz_surf = 1.5e-5    # m2/s 
    Kz_hlat = 100e-5    # m2/s 
    Kz_deep = 20e-5     # m2/s 
    # diffusion profile shape parameters
    trans_depth = 1000  # m
    trans_length = 1000 # m

    if remove_upw
        vdp = [[[hlat]    [hlat_bnd]   Kz_surf    Kz_hlat],    
            [[gyre]    [gyre_bnd]   Kz_surf    Kz_deep]] 
    else
        vdp = [[[hlat]    [hlat_bnd]   Kz_surf    Kz_hlat],    
            [[gyre]    [gyre_bnd]   Kz_surf    Kz_deep],    
            [[upw]     [upw_bnd]    Kz_surf    Kz_deep]] 
    end  

    num=length(vdp)
    for a in 1:num
        local i = vdp[a][1]                           # = hlat, gyre, upw e.g. [1,2,3]
        local bi = Int.(vdp[a][2][2:(length(vdp[a][2])-1)]) # internal box interfaces e.g. [2,3]

        local Kz = vdp[a][3] .+ (vdp[a][4]-vdp[a][3]) .* (1 .+ tanh.((bnd_depths[bi] .- trans_depth)/trans_length))/2

        local Flux = Hyps[bi] .* Kz ./ diff(depths[vdp[a][1]],dims=1) * 1e-6 # in Sv
        local Flux_k = [Flux; Flux; -[Flux;0]-[0;Flux]]# Flux with linear indeices

        # convert the 2D subscripts to linear indeices
        local k1 = [i[1:(length(i)-1)]; i[2:length(i)]; i]
        local k2 = [i[2:length(i)]; i[1:(length(i)-1)]; i]
        local k = (k2 .- 1) * ocean_cells .+ k1

        Z_mixing[k] = Z_mixing[k] + Flux_k
    end
    for i=1:ocean_cells
        if !(sum(Z_mixing[i,:])/Conv < 1e-6 && sum(Z_mixing[:,i])/Conv < 1e-6)
            @error "Not equal to zero error in the vertical diffusion"
        end
    end
    ################################################## # Vertical diffusion
    
    # Horizontal diffusion
    H_mixing = spzeros(ocean_cells, ocean_cells)
    Kh= 1e-3;           # Sv/m
    # Isopycnal eddy diffusivity, hlat->gyre
    Kh_hlat = 1e-3
    # Isopycnal eddy diffusivity, gyre->upw
    Kh_upw = 1e-3

    if circname == "Global_13_Box"
        num_box_surf_hlat = 1
        num_box_surf_gyre = 1
        index_750m_gyre = 3

    elseif circname == "Global_79_Box"
        # The ocean surface, hlat surface water (0-100m) and gyre's 750-1200m are isodense
        # No matter how many vertial boxes we have, the depth dependtent Isopycnal zone is constant
        num_box_surf_hlat = 1   # this is only a special case of num_box_surf_hlat = 1
        num_box_surf_gyre = 7
        index_750m_gyre = 15    # the index of the gyre's box with boundary depth = 750
    end

    i = sub2ind([ocean_cells, ocean_cells], [repeat([hlat[1]],num_box_surf_gyre)    gyre[index_750m_gyre:(index_750m_gyre-1+num_box_surf_gyre)]])
    # Sv          Sv/m = m2 yr-1      *        (m                                 + m)
    H_mixing[i] = Kh_hlat*0.5 .* (box_thickness[hlat[1]] / sum(box_thickness[gyre[index_750m_gyre:(index_750m_gyre-1+num_box_surf_gyre)]]) .* box_thickness[gyre[index_750m_gyre:(index_750m_gyre-1+num_box_surf_gyre)]]
                                    .+ box_thickness[gyre[index_750m_gyre:(index_750m_gyre-1+num_box_surf_gyre)]])

    # Mid to deep, hlat's (100-6000m) and gyre's 1200-6000m are isodense
    num_box_nosurf_hlat = length(hlat) - num_box_surf_hlat
    i = sub2ind([ocean_cells, ocean_cells], [hlat[(num_box_surf_hlat+1):end]        gyre[(index_750m_gyre+num_box_surf_gyre):(index_750m_gyre+num_box_surf_gyre-1+num_box_nosurf_hlat)]]) # hlat[2:13],gyre[22:33]
    # Sv          Sv/m = m2 yr-1      *        (m                                 + m)
    H_mixing[i] = Kh_hlat*0.5 .* (box_thickness[hlat[(num_box_surf_hlat+1):end]] 
                                    .+ box_thickness[gyre[(index_750m_gyre+num_box_surf_gyre):(index_750m_gyre+num_box_surf_gyre-1+num_box_nosurf_hlat)]])
    ####################

    if remove_upw
        # do nothing
        # Do gyre by matching with hlat
        for a in gyre
            for b in [hlat]
                H_mixing[a,b] = H_mixing[a,b] + H_mixing[b,a]
            end
        end
    else
        # The gyre and upw have the same number of boxes 
        i = sub2ind([ocean_cells, ocean_cells], [upw gyre])
        H_mixing[i] = Kh_upw*0.5 .* (box_thickness[upw] .+ box_thickness[gyre])

        # Do gyre by matching with hlat and upw
        for a in gyre
            for b in [hlat, upw]
                H_mixing[a,b] = H_mixing[a,b] + H_mixing[b,a]
            end
        end
    end

    # Do diagonal terms (net efflux from a box)
    for a in 1:ocean_cells
        H_mixing[a,a] = -sum(H_mixing[a,:])
    end

    for i=1:ocean_cells
        if !(sum(H_mixing[i,:])/Conv < 1e-6 && sum(H_mixing[:,i])/Conv < 1e-6)
            @error "Not equal to zero error in the vertical diffusion"
        end
    end
    ################################################## 

    # Compile the full transport matrix
    # V is a n=ocean_cells vector, use Conv./V to create a diagonal matrix
    # adv_mixing: unit is SV = = 10^6 m3/s
    Conv = 3.15576e16   # Conversion factor between Sverdrups and kg/yr assuming 1 m3 = 1000 kg (1Sv = 10^6 m3/s)
    # Conv = 1e6*1e3*365.25*24*60*60, 
    # yr-1  =  (normalized    ./  kg) * Sv   
    
    tm_rom      =  Diagonal(Conv ./  V) * (adv_mixing+Z_mixing+H_mixing+NonLocal_mixing)

    if remove_upw
        tm_shelf    =  Diagonal(Conv ./  V_shelf) * upw_mixing
        tm_shelf[1:ocean_cells, 1:ocean_cells] += tm_rom
        return tm_shelf, tm_rom
    else
        return tm_rom
    end
end

"""
Brillient Braille print of SparseArrays from: sparse_show.jl

https://gist.github.com/maxbennedich/5d32bf7fa94f1763f5362618067aef01

"""
function show_any_nonzero(S::SparseMatrixCSC; maxw = displaysize(stdout)[2], maxh = displaysize(stdout)[1]-3)
    h,w = size(S)
    h > 4maxh && (w = max(1, (w*4maxh+h÷2)÷h); h = 4maxh)
    w > 2maxw && (h = max(1, (h*2maxw+w÷2)÷w); w = 2maxw)
    P = fill(BRAILLE[1], (w+3)÷2, (h+3)÷4)
    P[end, :].=10
    @inbounds for c = 0:w-1, r = 0:h-1
        _anynz(S, r*S.m÷h+1, c*S.n÷w+1, (r+1)*S.m÷h, (c+1)*S.n÷w) &&
            (P[c÷2+1, r÷4+1] |= BRAILLE[4c&4+r&3+2])
    end
    print(join(Char.(P)))
    @info "number of nonzero = $(nnz(S))"
end

function _anynz(S::SparseMatrixCSC, r1, c1, r2, c2)
    @inbounds for c = c1:c2
        k = searchsortedfirst(S.rowval, r1, S.colptr[c], S.colptr[c+1]-1, Base.Order.Forward)
        k != S.colptr[c+1] && S.rowval[k] ≤ r2 && return true
    end
    false
end

end # module
