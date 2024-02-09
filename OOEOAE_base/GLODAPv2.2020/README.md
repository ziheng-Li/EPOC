# Data from GLobal Ocean Data Analysis Project Version 2.2020 (GLODAP)

You can find the data from [GLODAP](https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system/oceans/GLODAPv2_2020/)

You can download the data from [GLODAP_2020](https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0210813/)

### Variables
NB: the keys or column names in the Datasets have a prefix `G2` and some of them have postfix `f` or `qc`.
- Columns with no postfix is the real data
- Columns with `f` is the code of quality flag convention: WOCE quality control flags are used: `0 = interpolated`, `2 = good value`, `9 = no value`
- Columns with `qc`, I don't know, only guess 1 and 0 in them represent Observation type (measured_data_synthesis, or calculated)
- Column `expocode` and `expocodeno` are the code and number of the cruise, length != datapoints
 - μM = μmol/kg

Table 1. The list of variable names
|Name   |Abbreviation   | Unit|
|-------|---------------|-----|
|`Temperature`|temperature|˚C
|Potential temperature|theta|˚C
|Salinity|salinity|‰ (?)
|Potential density|sigma0|kg/m3
|Potential density, ref 1000 dbar|sigma1|kg/m3
|Potential density, ref 2000 dbar|sigma2|kg/m3
|Potential density, ref 3000 dbar|sigma3|kg/m3
|Potential density, ref 4000 dbar|sigma4|kg/m3
|Neutral density|gamma|kg/m3
|`Oxygen`|oxygen|μmol/kg
|Apparent oxygen utilization|aou|μmol/kg
|Nitrate|nitrate|μmol/kg
|Nitrite|nitrite|μmol/kg
|Silicate|silicate|μmol/kg
|`Phosphate`|phosphate|μmol/kg
|TCO2|tco2|μmol/kg
|TAlk|talk|μmol/kg
|fCO2|fco2|microatmospheres
|pH at total scale, 25 ◦ C and zero dbar of pressure|phts25p0|(?)
|pH at total scale, in situ temperature and pressure|phtsinsitutp|
|CFC-11|cfc11|pmol/kg
|pCFC-11|pcfc11|ppt
|CFC-12|cfc12|pmol/kg
|pCFC-12|pcfc12|ppt
|CFC-113|cfc113|pmol/kg
|CCl4|ccl4|pmol/kg
|pCCl4|pccl4|ppt
|SF6|sf6|fmol/kg
|pSF6|psf6|ppt
|`δ13C`|c13|‰
|`δ18O`|o18|‰
|∆14C|c14|‰
|∆14C counting error|c14err|‰
|3H|h3|TU
|3H counting error|h3err|TU
|he3|δ3He|%
|3He counting error|he3err|%
|he|He|nmol/kg
|He counting error|heerr|%
|Ne|neon|nmol/kg
|Ne counting error|neonerr|nmol/kg
|Total organic carbon|toc|μmol/L
|Dissolved organic carbon|doc|μmol/L
|Dissolved organic nitrogen|don|μmol/L
|Total dissolved nitrogen|tdn|μmol/L
|Chlorophyll a|chla|ug/L

### Plot the data

    julia> include("run_GLODAP.jl") # to activate environment
    julia> plot_TEMP_O2_P_dC13(["Arctic_Ocean"]) # plot raw data from Actic ocean
    julia> plot_TEMP_O2_P_dC13(["Atlantic_Ocean"  "Indian_Ocean"  "Pacific_Ocean"], "highlat", "oxygen") # plot high latitude data, with custom data filter

![Arctic_Ocean_only](Arctic_Ocean_only.png)
![Atlantic_Ocean_only](Atlantic_Ocean_only.png)
![Indian_Ocean_only](Indian_Ocean_only.png)
![Pacific_Ocean_only](Pacific_Ocean_only.png)
#### Fig. 1a-d. P_O2_TEMP_d13C data from diffrent oceans.

### Data comparsion Model vs. data

To achieve this, you mush pre-run a top-level script to get the `paleorun.output`. Which save the output of a model run, e.g.

    julia> include("OOEOAE_base/P_O_romglb")

Then run the custom function to make the plot

    julia> plot_compare_GLODAP(paleorun.output) # src in the base_GLODAP.jl

![data_vs_model](Romglb_vs_GLODAP_GeoClim.png)
#### Fig 2. GLODAP data vs. model output, data cloud filtered following [GEOCLim]().

### Data Citation
Whenever GLODAPv2.2020 is used two papers should be cited: [Olsen et al. (2019)](https://doi.org/10.5194/essd-2019-66) and [Olsen et al. (2020)](https://essd.copernicus.org/articles/12/3653/2020/), where the latter which describes the procedures for the update.

And other citations:

Olsen, A., R. M. Key, S. van Heuven, S. K. Lauvset, A. Velo, X. Lin, C. Schirnick, A. Kozyr, T. Tanhua, M. Hoppema, S. Jutterström, R. Steinfeldt, E. Jeansson, M. Ishii, F. F. Pérez and T. Suzuki. The Global Ocean Data Analysis Project version 2 (GLODAPv2) – an internally consistent data product for the world ocean, Earth Syst. Sci. Data, 8, 297–323, 2016, doi:10.5194/essd-8-297-2016.

Lauvset, S. K, R. M. Key, A. Olsen, S. van Heuven, A. Velo, X. Lin, C. Schirnick, A. Kozyr, T. Tanhua, M. Hoppema, S. Jutterström, R. Steinfeldt, E. Jeansson, M. Ishii, F. F. Pérez, T. Suzuki and S. Watelet. A new global interior ocean mapped climatology: the 1°x1° GLODAP version 2, Earth Syst. Sci. Data, 8, 325–340, 2016, doi:10.5194/essd-8-325-2016.

# Data from World Ocean Atlas (WOA2018)

You can find the data from [WOA18](https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/)

- woa18_all_o00_01.nc, dissolved Oxygen
- woa18_all_p00_01.nc, Phosphate
- woa18_decav_t00_01.nc, averaged decades surface temperature
- woa18_decav_t13_01.nc, averaged decades (winter only) surface temperature

To read the .nc file

    julia> pathname = joinpath(@__DIR__, "woa18_decav_t00_01.nc")
    julia> ds = NCDataset(pathname)

Make the plots using this file

    julia> contourf(ds["lon"][:], ds["lat"][:], ds["t_an"][:,:,1,1]', clims=(-5,4),) # show the area with average TEMP < 4 degree

![TEMP_<_4](WOA_TEMP_4_degree.png)
#### Fig. 3. WOA data, average TEMP < 4 degree.