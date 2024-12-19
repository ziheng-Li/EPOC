# EPOC Model

To install (Julia 1.10 required):

Change directory to examples

    julia> ] activate .
    (examples) pkg> instantiate
    (examples) pkg> add https://github.com/PALEOtoolkit/NLsolve.jl#project_region
    (examples) pkg> precompile

# Examples

## Daines & Li (2024) 'theory paper'

examples/dainesli2024/


see examples/dainesli2024/README.md for work-in-progress 

## Li etal (2024) 'Gaskiers paper'

examples/3_P_O_A_U_columns_Ziheng_Egan_Formation/

## Romaniello + shelves configurations

examples/P_O_romglb_shelf_20231218/

# EPOC model code

src/

# Notes and TODO

This repo is a partial copy of PALEOexamples\src\OOEOAE_base folder from
PALEOdev.jl repo, branch OOE_2021_Ziheng_dev2

From:
    commit 74fa0878494557fd969518c637bec46e6e2e43d9 (HEAD -> OOE_2021_Ziheng_dev2, origin/OOE_2021_Ziheng_dev2)
    Author: Ziheng Li <zihengli@cug.edu.cn>
    Date:   Sun Jan 28 15:00:46 2024 +0800

(Daines & Li 2024 PNAS submission)

with some work-in-progress partial updates to later versions

## Code changes to run outside PALEOdev.jl repo

Remove import PALEOreactions (part of PALEOdev.jl repository)

Copy from PALEOdev/PALEOreactions to OOEOAE_base/PALEOreactions:
- Uranium.jl, AtmReservoirs.jl, Burial.jl, OceanTransportRomanielloShelf.jl

Copy from OOEOAE_base:
- OOEOAE_base/CarbBurial_dev.jl 
    from: PALEOcopse.jl\src\oceanfloor\CarbBurial.jl
        Update the mccb for each cell...
    and include in each of the P-O-A scripts 

- ReactionsOOEOAE_dev.jl
    - ReactionOceanBurialColumn: column ocean hierarchy 
    - ReactionOxWeathMinimal: simplified oxidw
    - ReactionUOceanfloor_dev: Uburial for each cell...

- SolverFunctionsOOEOAE2.jl
    Structure and functions to find the nullclines and other items in Phase plane.

- Isoline.jl
    src to support SolverFunctionsOOEOAE2.jl

- ooeoae_plots.jl
    Romaniello ocean plotting functions

