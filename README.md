# energy_autarky

This projects estimates different kinds of electricity supply and demand for each administrative unit (from municipality to continent) in Europe.

## Manually downloaded data

The ERA5 reanalysis data and geoboundaries data are downloaded automatically, but other data need to be retrieved manually. The detailed inventory and sources for manual data are listed in the Table 1 in the appendix in the thesis.

For Exclusive Economic Zone data,

    1. Navigate to the website https://www.marineregions.org/downloads.php.
    2. Select Maritime Boundaries v11 
    3. Select World EEZ v11 - Shapefile
    4. Decompress the downloaded file to directory "data/maritime_data/"
    5. Make sure the downloaded file name and suffix are the same as the input file in snakemake rule,
        In "rules/energy_supply.smk",
        Replace "data/maritime_data/EEZ.geojson" to the new directory

For Maritime military areas, Offshore oil and gas wells, Maritime traffic density map and Offshore pipelines data,

    1. Navigate to the website https://emodnet.ec.europa.eu/geoviewer/
    2. In "Layers", select EMODnet Human Activities
    3. For Maritime military areas, select Military Areas - Military Areas (Polygons) - i - download
       For Offshore oil and gas wells, select Oil and Gas - Offshore Installations - i - download
       For Maritime traffic density map, select Route Density - All vessels (Annual totals 2019-2022) - i - download
       For Offshore pipelines, select Pipelines - Pipelines - i - download
    4. Decompress the downloaded file to directory "data/maritime_data/"
    5. Make sure the downloaded file name and suffix are the same as the input files in snakemake rules,
        In "rules/energy_supply.smk",
        For Maritime military areas, replace "data/maritime_data/maritime_military_areas.geojson" with the new directory
        For Offshore oil and gas wells, replace "data/maritime_data/offshore_oil_and_gas.geojson" with the new directory
        For Maritime traffic density map, replace "data/maritime_data/wid6-all_traffic-all_europe-yearly-20190101000000_20191231235959-tdm-grid.tif" with the new directory
        For Offshore pipelines, replace "data/maritime_data/offshore_pipelines.geojson" with the new directory

For shoreline data,

    1. Navigate to https://www.eea.europa.eu/data-and-maps/data/eea-coastline-for-analysis-2
    2. Find EEA coastline - Polyline, click Download file
    3. Decompress the downloaded file to directory "data/maritime_data/"
    4. Make sure the downloaded file name and suffix are the same as the input files in snakemake rules,
        In "rules/energy_supply.smk",
        Replace "data/maritime_data/shoreline.geojson" with the new directory

For Protected areas data,

    1. Navigate to https://www.protectedplanet.net./en/thematic-areas/wdpa?tab=WDPA
    2. Click Download
    3. Decompress the downloaded file to directory "data/maritime_data/"
    4. Make sure the downloaded file name and suffix are the same as the input files in snakemake rules,
        In "rules/energy_supply.smk",
        Replace "data/maritime_data/protected_areas/WDPA_Oct2022_Public_shp_0/WDPA_Oct2022_Public_shp-polygons.shp", 
        "data/maritime_data/protected_areas/WDPA_Oct2022_Public_shp_1/WDPA_Oct2022_Public_shp-polygons.shp",
        "data/maritime_data/protected_areas/WDPA_Oct2022_Public_shp_2/WDPA_Oct2022_Public_shp-polygons.shp",
        with the new directory

For European settlement map data,

    1. Navigate to https://land.copernicus.eu/pan-european/GHSL/european-settlement-map/esm-2012-release-2017-urban-green?tab=download
    2. Type "10m" in the Search tab
    3. Download and decompress all the data into "data/ESM/"

For Corine Land Cover data,

    1. Naviagte to https://land.copernicus.eu/pan-european/corine-land-cover/clc2018?tab=download
    2. Download and decompress all the data into "data/corine_land_cover/"
        (The directory of the tif data should be data/corine_land_cover/DATA/U2006_CLC2000_V2020_20u1.tif)

For Bathymetry data,

    1. Navigate to https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/georeferenced_tiff/
    2. Download and decompress all the data into "data/maritime_data/"
        (The directory of the tif data should be data/maritime_data/ETOPO1_Bed_g_geotiff.tif)

For Population density map data,

    1. Navigate to https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/population-distribution-demography/geostat
    2. Download JRC_GRID_2018.zip
    3. Decompress all the data into "data/population_density/"
        (The directory of the tif data should be data/population_density/JRC_1K_POP_2018.tif)

For Prior data of the GLAES tool,

    1. Navigate to https://data.mendeley.com/datasets/trvfb3nwt2/1
    2. Download all the tif data into "data/excludePrior"

For Renewables.ninja simulation data,

    1. Navigate to https://zenodo.org/record/6559895#.ZDo38HZBxPb
    2. Download all the nc data into "data/renewables_ninja "

For JRC-IDEES data,

    1. Navigate to https://data.jrc.ec.europa.eu/dataset/jrc-10110-10001%E2%80%8C
    2. Download JRC-IDEES 2015 data into "data/JRC-IDEES"

For ENTSOE-E historical generation data,

    1. Navigate to https://transparency.entsoe.eu/generation/r2/actualGenerationPerProductionType/show
    2. Download the annual generation per type for each studied country into "data/Entsoe/(3-letter-country-code)"

    Note: For convenience, it is recommended to use the already-downloaded data

For Historical electricity generation of Albania,

    1. Navigate to https://www.instat.gov.al/al/temat/mjedisi-dhe-energjia/energjia/publikimet/2022/bilanci-i-energjis%C3%AB-elektrike-tremujori-i-2022
    2. Download the annual generation per type for Albania into "data/Entsoe/ALB/bilanci-i-energjisë-elektrike-tremujorë.xls"

    Note: For convenience, it is recommended to use the already-downloaded data

For Global power plant database,

    1. Navigate to https://datasets.wri.org/dataset/globalpowerplantdatabase
    2. Download Global Power Plant Database v1.3.0
    3. Put the downloaded file into "data/global_power_plant_database_v_1_3/global_power_plant_database.xls"

For ETS industrial emission inventory,

    1. Navigate to https://www.eea.europa.eu/data-and-maps/data/industrial-reporting-under-the-industrial-6
    2. Download and decompress the data into "data/others/industrial emission database.xlsx"

    Note This version may be deprecated, in that case, it is recommended to use the already-downloaded data

For Eurostat data,

    1. Navigate to https://ec.europa.eu/eurostat/databrowser/explore/all/all_themes
    2. For biomass data, please navigate to:
        https://ec.europa.eu/eurostat/databrowser/view/ef_lpc_fruit/default/table?lang=en, 
        https://ec.europa.eu/eurostat/databrowser/view/ef_lpc_vineyd/default/table?lang=en, 
        https://ec.europa.eu/eurostat/databrowser/view/tag00047/default/table?lang=en.
        https://ec.europa.eu/eurostat/databrowser/view/tag00049/default/table?lang=en.
        https://ec.europa.eu/eurostat/databrowser/view/tag00051/default/table?lang=en.
        https://ec.europa.eu/eurostat/databrowser/view/tag00053/default/table?lang=en.
        https://ec.europa.eu/eurostat/databrowser/view/tag00093/default/table?lang=en.
        https://ec.europa.eu/eurostat/databrowser/view/tag00110/default/table?lang=en.
        https://ec.europa.eu/eurostat/databrowser/view/tag00120/default/table?lang=en.
        https://ec.europa.eu/eurostat/databrowser/view/tag00122/default/table?lang=en.
        https://ec.europa.eu/eurostat/databrowser/view/ten00108/default/table?lang=en.
        https://ourworldindata.org/grapher/rice-production.
        
        Download and Decompress the data into "data/biomass"

    3. For House to flat ratio data,

        Navigate to https://ec.europa.eu/eurostat/databrowser/view/ilc_lvho01/default/table?lang=en
        Download and decompress the data into "data/others/ilc_lvho01_linear.csv"

    4. For Final energy consumption by sector data,

        Navigate to https://ec.europa.eu/eurostat/databrowser/view/ten00124/default/table?lang=en
        Download and decompress the data into "data/others/ten00124_linear.csv"
    
    5. For National population data,

        Navigate to https://ec.europa.eu/eurostat/databrowser/view/tps00001/default/table?lang=en.
        Download and decompress the data into "data/others/tps00001_page_linear.csv"

    Note This process may be too lengthy, it is recommended to use the already-downloaded data

For Heating threshold temperatures data,

    The data are from the paper https://doi.org/10.1016/j.enbuild.2019.07.013.
    It is recommended to use the already-downloaded data in "data/others/heating_threshold_temperatures.csv"

For When2heat parameters,

    The data are from https://github.com/oruhnau/when2heat/tree/2022-02-22/input
    It is recommended to use the already-downloaded data in "data/when2heat"


## Getting ready

You need [mamba](https://mamba.readthedocs.io/en/latest/) to run the analysis. Using mamba, you can create an environment from within you can run it:

    mamba env create -f environment.yaml --no-default-packages

## Run the analysis

    snakemake --profile profiles/default

This will run all analysis steps to reproduce results and eventually build the report.

You can also run certain parts only by using other `snakemake` rules; to get a list of all rules run `snakemake --list`.

To generate a PDF of the dependency graph of all steps `build/dag.pdf` run:

    snakemake --profile profiles/default -f dag

To obtain the detailed inventory of energy supply and demand for national results, run:

    snakemake  --profile profiles/default --rerun-incomplete preprocess_national_visualization_data

    For regional results, run:

    snakemake  --profile profiles/default --rerun-incomplete preprocess_regional_visualization_data

    For municipal results, run:

    snakemake  --profile profiles/default --rerun-incomplete preprocess_municipal_visualization_data

To obtain the first two review graphs in the thesis, run:

    snakemake  --profile profiles/default --rerun-incomplete plot_potential_demand_1 
    snakemake  --profile profiles/default --rerun-incomplete plot_potential_demand_2

To obtain the choropleth maps, run:

    snakemake  --profile profiles/default --rerun-incomplete preprocess_national_visualization_data
    snakemake  --profile profiles/default --rerun-incomplete preprocess_regional_visualization_data
    snakemake  --profile profiles/default --rerun-incomplete preprocess_municipal_visualization_data

To obtain the boxplots, run:

    snakemake  --profile profiles/default --rerun-incomplete get_municipal_boxplots
    snakemake  --profile profiles/default --rerun-incomplete get_regional_boxplots

To obtain the national energy structure, run:

    snakemake  --profile profiles/default --rerun-incomplete plot_national_energy_structure

## Be notified of build successes or fails

  As the execution of this workflow may take a while, you can be notified whenever the execution terminates either successfully or unsuccessfully. Notifications are sent by email. To activate notifications, add the email address of the recipient to the configuration key `email`. You can add the key to your configuration file, or you can run the workflow the following way to receive notifications:

      snakemake --profile profiles/default --config email=<your-email>

## Run the tests

    snakemake --profile profiles/default test

## Repo structure

* `report`: contains all files necessary to build the report; plots and result files are not in here but generated automatically
* `scripts`: contains the Python source code as scripts
* `rules`: contains Snakemake rule definitions
* `envs`: contains execution environments
* `tests`: contains the test code
* `config`: configurations used in the study
* `profiles`: Snakemake execution profiles
* `data`: place for raw data
* `build`: will contain all results (does not exist initially)

## License

The code in this repo is MIT licensed, see `./LICENSE.md`.
