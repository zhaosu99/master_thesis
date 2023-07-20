from pathlib import Path
import requests


def func_download_national_municipalities(wildcards):

    adm_list = [str(i) for i in range(5, 0, -1)]
    url_list = ["https://www.geoboundaries.org/data/geoBoundaries-3_0_0/{gid}/ADM{adm}/geoBoundaries-3_0_0-{gid}-ADM{adm}.geojson".
        format(gid=wildcards.gid, adm=x) for x in adm_list]
    for i in url_list:
        if requests.get(i).status_code == 200:
            wanted_url = i
            break
        else:
            continue

    return wanted_url


rule download_national_municipalities:
    message: "Download lowest ADM level for {wildcards.gid}"
    params: func_download_national_municipalities
    output:
        "data/geoboundaries/{gid}_municipalities.geojson"
    conda: "../envs/default.yaml"
    shell:
        "curl -o {output} '{params}'"


rule preprocess_national_municipalities:
    message: "Standardize national municipalities data for {wildcards.gid}"
    input:
        "data/geoboundaries/{gid}_municipalities.geojson"
    output:
        "temp/geoboundaries/{gid}_municipalities.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/preprocess_municipalities.py"


rule download_national_regions:
    message: "Download ADM 1 level for {wildcards.gid}"
    params: "https://www.geoboundaries.org/data/geoBoundaries-3_0_0/{gid}/ADM1/geoBoundaries-3_0_0-{gid}-ADM1.geojson"
    output:
        "data/geoboundaries/{gid}_regions.geojson"
    conda: "../envs/default.yaml"
    shell:
        "curl -o {output} '{params}'"


rule preprocess_national_regions:
    message: "Standardize national regions data for {wildcards.gid}"
    input:
        "data/geoboundaries/{gid}_regions.geojson",
        "data/geoboundaries/GBR_simplified/geoBoundaries-GBR-ADM1_simplified.geojson"
    output:
        "temp/geoboundaries/{gid}_regions.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/preprocess_regions.py"


rule split_national_regions:
    message: "Split all national regions to individuals"
    input:
        national_regions = ["temp/geoboundaries/{gid}_regions.geojson".format(gid=x) for x in config['country_codes']]
    output:
        regions = ["temp/geoboundaries/{gid}/{gid}_region_{n}.geojson".format(gid=x, n=y)
            for x in config['country_codes']
            for y in range(config['regions_amount'][x])]
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/split_country_to_regions.py"


rule create_onshore_voronoi_polygons_with_site_id:
    message: "Create a voronoi geodataframe with the points from existing renewables.ninja onshore simulations"
    input:
        "data/renewables_ninja/wind-onshore-timeseries.nc"
    output:
        "temp/renewables_ninja_data/created_onshore_voronoi_polygons_with_site_id.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/create_voronoi_polygons_with_site_id.py"


rule split_municipalities_by_site_id:
    message: "Split {wildcards.gid} municipalities according to their belonging site id"
    input:
        "temp/geoboundaries/{gid}_municipalities.geojson",
        "temp/renewables_ninja_data/created_onshore_voronoi_polygons_with_site_id.geojson"
    output:
        "temp/geoboundaries/{gid}_municipalities_with_site_id.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/split_municipalities_by_site_id.py"


rule create_offshore_voronoi_polygons_with_site_id:
    message: "Create a voronoi geodataframe with the points from existing renewables.ninja offshore simulations"
    input:
        "data/renewables_ninja/wind-offshore-timeseries.nc"
    output:
        "temp/renewables_ninja_data/created_offshore_voronoi_polygons_with_site_id.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/create_voronoi_polygons_with_site_id.py"


rule create_csv_with_site_id_and_capacity_factor:
    message: "Restructure the existing ncs to dataframes with site IDs and monthly capacity factors"
    input:
        "data/renewables_ninja/open-field-pv-timeseries.nc",
        "data/renewables_ninja/rooftop-pv-e-w-timeseries.nc",
        "data/renewables_ninja/rooftop-pv-n-timeseries.nc",
        "data/renewables_ninja/rooftop-pv-s-flat-timeseries.nc",
        "data/renewables_ninja/wind-onshore-timeseries.nc",
        "data/renewables_ninja/wind-offshore-timeseries.nc"
    output:
        "temp/renewables_ninja_data/open_field_solar_PV_site_id_with_capacity_factor.csv",
        "temp/renewables_ninja_data/rooftop_solar_PV_e_w_site_id_with_capacity_factor.csv",
        "temp/renewables_ninja_data/rooftop_solar_PV_n_site_id_with_capacity_factor.csv",
        "temp/renewables_ninja_data/rooftop_solar_PV_s_flat_site_id_with_capacity_factor.csv",
        "temp/renewables_ninja_data/onshore_wind_site_id_with_capacity_factor.csv",
        "temp/renewables_ninja_data/offshore_wind_site_id_with_capacity_factor.csv"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/create_csv_with_site_id_and_capacity_factor.py"


rule generate_synthesized_capacity_factor_csv:
    message: "Generate synthesized average monthly capacity factor for rooftop solar PV according to the share of buildings towards each direction"
    input:
        EW = "temp/renewables_ninja_data/rooftop_solar_PV_e_w_site_id_with_capacity_factor.csv",
        N = "temp/renewables_ninja_data/rooftop_solar_PV_n_site_id_with_capacity_factor.csv",
        SFLAT = "temp/renewables_ninja_data/rooftop_solar_PV_s_flat_site_id_with_capacity_factor.csv"
    output:
        "temp/renewables_ninja_data/rooftop_solar_PV_synthesized_site_id_with_capacity_factor.csv"
    params:
        EW_ratio = config["parameters"]["west_and_east_oriented_rooftops_ratio"],
        N_ratio = config["parameters"]["north_oriented_rooftops_ratio"],
        SFLAT_ratio = config["parameters"]["flat_and_south_oriented_rooftops_ratio"]
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/generate_synthesized_capacity_factor_csv.py"


def func_generate_ESM_bounds(wildcards):
    return list(Path("data/ESM").rglob("*.tif"))


rule generate_ESM_bounds:
    message: "Create a dataframe file containing the bounds of input ESM rasters"
    input:
        ESM_list = func_generate_ESM_bounds
    output:
        "temp/ESM/ESM_bounds.csv"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/generate_ESM_bounds.py"


rule preprocess_eez:
    message: "Select EEZ of studied countries and exclude too distant sea areas"
    input:
        "data/maritime_data/EEZ.geojson",
        "data/maritime_data/shoreline.geojson"
    output:
        eez_by_country = ["temp/EEZ/{gid}/{gid}_EEZ.geojson".format(gid=x) for x in config["country_codes"]]
    params:
        country_codes = config["country_codes"],
        max_distance_to_land = config['parameters']['max_offshore_turbine_distance_to_land']
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/preprocess_EEZ.py"


rule preprocess_offshore_protected_areas:
    message: "Select the marine proteced areas of studied countries"
    input:
        protected_areas = ["data/maritime_data/protected_areas/WDPA_Oct2022_Public_shp_{n}/WDPA_Oct2022_Public_shp-polygons.shp".format(n=x)
            for x in [0, 1, 2]],
        eez_by_country = ["temp/EEZ/{gid}/{gid}_EEZ.geojson".format(gid=x) for x in config["country_codes"]]
    output:
        offshore_protected_areas_by_country = ["temp/EEZ/{gid}/{gid}_offshore_protected_areas.geojson".format(gid=x) for x in config["country_codes"]]
    params:
        config["country_codes"]
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/preprocess_protected_areas.py"


rule preprocess_offshore_military_areas:
    message: "Select valid military areas"
    input:
        "data/maritime_data/maritime_military_areas.geojson",   
    output:
        "temp/maritime_data/preprocessed_offshore_military_areas.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/preprocess_offshore_military_areas.py"


rule preprocess_offshore_oil_and_gas:
    message: "Select valid offshore oil and gas wells"
    input:
        "data/maritime_data/offshore_oil_and_gas.geojson",   
    output:
        "temp/maritime_data/preprocessed_offshore_oil_and_gas.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/preprocess_offshore_oil_and_gas.py"


rule preprocess_offshore_pipelines:
    message: "Select valid offshore pipelines"
    input:
        "data/maritime_data/offshore_pipelines.geojson"      
    output:
        "temp/maritime_data/preprocessed_offshore_pipelines.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/preprocess_offshore_pipelines.py"


rule preprocess_bathymetry:
    message: "Reproject bathymetry data to EPSG: 3035"
    input:
        "data/maritime_data/ETOPO1_Bed_g_geotiff.tif"
    output:
        "temp/maritime_data/bathymetry_4326.tif",
        "temp/maritime_data/bathymetry_3035.tif"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/preprocess_bathymetry.py"


rule mask_bathymetry:
    message: "Mask bathymetry data to the EEZ of {wildcards.gid}"
    input:
        "temp/EEZ/{gid}/{gid}_EEZ.geojson",
        "temp/maritime_data/bathymetry_3035.tif"
    output:
        "temp/EEZ/{gid}/{gid}_bathymetry.tif"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/split_bathymetry.py"


rule offshore_wind_glaes:
    message: "Conduct sea exclusion calculation for {wildcards.gid} EEZ to select offshore wind turbine sites"
    input:
        eez_by_country = "temp/EEZ/{gid}/{gid}_EEZ.geojson",
        offshore_protected_areas_by_country = "temp/EEZ/{gid}/{gid}_offshore_protected_areas.geojson",
        bathymetry_by_country = "temp/EEZ/{gid}/{gid}_bathymetry.tif",
        offshore_military_areas = "temp/maritime_data/preprocessed_offshore_military_areas.geojson",
        offshore_oil_and_gas = "temp/maritime_data/preprocessed_offshore_oil_and_gas.geojson",
        offshore_pipelines = "temp/maritime_data/preprocessed_offshore_pipelines.geojson",
        offshore_traffic = "data/maritime_data/wid6-all_traffic-all_europe-yearly-20190101000000_20191231235959-tdm-grid.tif",
        shoreline = "data/maritime_data/shoreline.geojson"
    output:
        "temp/glaes/offshore_wind/offshore_turbine_locations_{gid}.csv"
    params:
        offshore_turbine_distance = config["parameters"]["offshore_turbine_separation_distance"]
    conda: "../envs/glaes.yaml"
    script:
        "../scripts/energy_supply/generate_offshore_wind_turbine_locations.py"


rule calculate_offshore_wind_potential:
    message: "Summarize the monthly average capacity factors and suitable area size of offshore wind for relevant municipals in {wildcards.gid}"
    input:
        offshore_turbines = "temp/glaes/offshore_wind/offshore_turbine_locations_{gid}.csv",
        municipal_geojson = "temp/geoboundaries/{gid}_municipalities.geojson",
        site_id_geojson = "temp/renewables_ninja_data/created_offshore_voronoi_polygons_with_site_id.geojson",
        capacity_factor_csv = "temp/renewables_ninja_data/offshore_wind_site_id_with_capacity_factor.csv"
    output:
        "temp/offshore_wind/{gid}/offshore_wind_potential_{gid}.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/calculate_offshore_wind_potential.py"


rule onshore_wind_glaes:
    message: "Conduct land exclusion calculation for each municipal in {wildcards.gid} region {wildcards.n} to select onshore wind turbine sites"
    input:
        municipal_geojson = "temp/geoboundaries/{gid}/{gid}_region_{n}.geojson"
    output:
        "temp/glaes/onshore_wind/{gid}/onshore_turbine_locations_{gid}_{n}.csv"
    params:
        prior_datasets_path = "data/excludePrior",
        onshore_turbine_distance = config["parameters"]["onshore_turbine_separation_distance"]
    conda: "../envs/glaes.yaml"
    script:
        "../scripts/energy_supply/generate_onshore_wind_turbine_locations.py"


rule generate_onshore_wind_turbine_locations_geojson:
    message: "Convert {wildcards.gid} region {wildcards.n}'s onshore wind sites csv to geojson file"
    input:
        "temp/glaes/onshore_wind/{gid}/onshore_turbine_locations_{gid}_{n}.csv"
    output:
        "temp/glaes/onshore_wind/{gid}/onshore_turbine_locations_{gid}_{n}.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/convert_csv_to_geojson.py"


rule calculate_onshore_wind_potential:
    message: "Summarize the monthly average capacity factors and onshore turbine amount for each municipal in {wildcards.gid} region {wildcards.n}"
    input:            
        onshore_turbines = "temp/glaes/onshore_wind/{gid}/onshore_turbine_locations_{gid}_{n}.geojson",
        municipal_geojson_with_site_id = "temp/geoboundaries/{gid}_municipalities_with_site_id.geojson",
        national_municipalities = "temp/geoboundaries/{gid}_municipalities.geojson",
        capacity_factor_csv = "temp/renewables_ninja_data/onshore_wind_site_id_with_capacity_factor.csv"
    output:
        "temp/onshore_wind/{gid}/onshore_wind_potential_{gid}_{n}.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/calculate_onshore_wind_potential.py"


rule open_field_solar_PV_glaes:
    message: "Conduct land exclusion calculation for each municipal in {wildcards.gid} region {wildcards.n} to select suitable open field solar PV areas"
    input:
        onshore_turbine_locations_geojson = "temp/glaes/onshore_wind/{gid}/onshore_turbine_locations_{gid}_{n}.geojson",
        region =  "temp/geoboundaries/{gid}/{gid}_region_{n}.geojson"
    output:
        "temp/glaes/open_field_solar_PV/{gid}/open_field_solar_PV_{gid}_{n}.tif"
    params:
        prior_datasets_path = "data/excludePrior"
    conda: "../envs/glaes.yaml"
    script:
        "../scripts/energy_supply/generate_open_field_solar_PV_areas.py"


rule calculate_open_field_solar_PV_potential:
    message: "Summarize the monthly average capacity factors and suitable area size of open field solar PV for each municipal in {wildcards.gid} region {wildcards.n}"
    input:
        national_municipalities_geojson = "temp/geoboundaries/{gid}_municipalities.geojson",
        municipal_geojson_with_site_id = "temp/geoboundaries/{gid}_municipalities_with_site_id.geojson",
        region =  "temp/geoboundaries/{gid}/{gid}_region_{n}.geojson",
        capacity_factor_csv = "temp/renewables_ninja_data/open_field_solar_PV_site_id_with_capacity_factor.csv",
        open_field_solar_PV_areas = "temp/glaes/open_field_solar_PV/{gid}/open_field_solar_PV_{gid}_{n}.tif",     
    output:
        "temp/open_field_solar_PV/{gid}/open_field_solar_PV_potential_{gid}_{n}.geojson"
    params: value_of_wanted_cells = [100]
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/calculate_open_field_solar_PV_potential.py"


rule calculate_rooftop_solar_PV_potential:
    message: "Summarize the monthly average capacity factors and suitable area size of rooftop solar PV for each municipal in {wildcards.gid} region {wildcards.n}"
    input:
        national_municipalities_geojson = "temp/geoboundaries/{gid}_municipalities.geojson",
        municipal_geojson_with_site_id = "temp/geoboundaries/{gid}_municipalities_with_site_id.geojson",
        region =  "temp/geoboundaries/{gid}/{gid}_region_{n}.geojson",
        ESM_bounds_csv = "temp/ESM/ESM_bounds.csv",      
        ESM_list = func_generate_ESM_bounds,       
        synthesized_capacity_factor_csv = "temp/renewables_ninja_data/rooftop_solar_PV_synthesized_site_id_with_capacity_factor.csv"
    output:
        "temp/rooftop_solar_PV/{gid}/rooftop_solar_PV_potential_{gid}_{n}.geojson"
    params: value_of_wanted_cells = [255, 250]
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/calculate_rooftop_solar_PV_potential.py"


def func_concat_national_onshore_wind(wildcards):
    return ["temp/onshore_wind/{gid}/onshore_wind_potential_{gid}_{n}.geojson".format(gid=wildcards.gid, n=x)
            for x in range(config['regions_amount'][wildcards.gid])]


def func_concat_national_open_field_solar_PV(wildcards):
    return ["temp/open_field_solar_PV/{gid}/open_field_solar_PV_potential_{gid}_{n}.geojson".format(gid=wildcards.gid, n=x)
            for x in range(config['regions_amount'][wildcards.gid])]


def func_concat_national_rooftop_solar_PV(wildcards):
    return ["temp/rooftop_solar_PV/{gid}/rooftop_solar_PV_potential_{gid}_{n}.geojson".format(gid=wildcards.gid, n=x)
            for x in range(config['regions_amount'][wildcards.gid])]


rule concat_national_onshore_wind:
    message: "Concatenate {wildcards.gid}'s regional onshore wind electricity supply to national level"
    input:
        func_concat_national_onshore_wind
    output:
        "temp/onshore_wind/{gid}/onshore_wind_potential_{gid}.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/concat_regional_results_to_national.py"


rule concat_national_open_field_solar_PV:
    message: "Concatenate {wildcards.gid}'s regional open field solar PV electricity supply to national level"
    input:
        func_concat_national_open_field_solar_PV
    output:
        "temp/open_field_solar_PV/{gid}/open_field_solar_PV_potential_{gid}.geojson",
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/concat_regional_results_to_national.py"


rule concat_national_rooftop_solar_PV:
    message: "Concatenate {wildcards.gid}'s regional rooftop solar PV electricity supply to national level"
    input:
        func_concat_national_rooftop_solar_PV
    output:
        "temp/rooftop_solar_PV/{gid}/rooftop_solar_PV_potential_{gid}.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/concat_regional_results_to_national.py"


def func_monthly_generation_input(wildcards):
    return list(Path("data/Entsoe/{gid}".format(gid=wildcards.gid)).glob("*.csv"))


rule calculate_monthly_generation_of_other_electricity:
    message: "Calculate the average monthly generation of hydro, nuclear and geothermal for {wildcards.gid}"
    input:
        monthly_generation_input = func_monthly_generation_input,
        ALB_seasonal_generation = "data/Entsoe/ALB/bilanci-i-energjisë-elektrike-tremujorë.xls"
    output:
        "temp/other_electricity/other_electricity_monthly_generation/other_electricity_monthly_generation_{gid}.csv"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/calculate_monthly_generation_of_other_electricity.py"


rule split_hydro:
    message: "Split raw dataset and create hydro power plant sites for relevant countries"
    input:
        "data/global_power_plant_database_v_1_3/global_power_plant_database.csv",
        ["temp/geoboundaries/{gid}_municipalities.geojson".format(gid=x) for x in config["country_codes"]]
    output:
        output_list = ["temp/other_electricity/hydro/hydro_power_plant_sites_geojson_{gid}.geojson".format(gid=x) for x in config["country_codes"]]
    params: "Hydro",
            config["country_codes"]
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/split_other_clean_electricity.py"


rule calculate_hydro_potential:
    message: "Summarize the installed capacity and average monthly generation of hydro power for relevant municipals in {wildcards.gid}"
    input:
        power_plant_sites_geojson = "temp/other_electricity/hydro/hydro_power_plant_sites_geojson_{gid}.geojson",
        monthly_generation = "temp/other_electricity/other_electricity_monthly_generation/other_electricity_monthly_generation_{gid}.csv",
        municipals = "temp/geoboundaries/{gid}_municipalities.geojson"
    output:
        "temp/hydro/{gid}/hydro_potential_{gid}.geojson"
    params:
        'hydro_MWh'
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/calculate_potential_of_other_electricity.py"


rule split_nuclear:
    message: "Split raw dataset and create nuclear power plant sites for relevant countries"
    input:
        "data/global_power_plant_database_v_1_3/global_power_plant_database.csv"
    output:
        output_list = ["temp/other_electricity/nuclear/nuclear_power_plant_sites_geojson_{gid}.geojson".format(gid=x) for x in config["country_codes"]]
    params: "Nuclear",
            config["country_codes"]
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/split_other_clean_electricity.py"


rule calculate_nuclear_potential:
    message: "Summarize the installed capacity and average monthly generation of nuclear power for relevant municipals in {wildcards.gid}"
    input:
        power_plant_sites_geojson = "temp/other_electricity/nuclear/nuclear_power_plant_sites_geojson_{gid}.geojson",
        monthly_generation = "temp/other_electricity/other_electricity_monthly_generation/other_electricity_monthly_generation_{gid}.csv",
        municipals = "temp/geoboundaries/{gid}_municipalities.geojson"
    output:
        "temp/nuclear/{gid}/nuclear_potential_{gid}.geojson"
    params:
        'nuclear_MWh'
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/calculate_potential_of_other_electricity.py"


rule split_geothermal:
    message: "Split raw dataset and create geothermal power plant sites for relevant countries"
    input:
        "data/global_power_plant_database_v_1_3/global_power_plant_database.csv"
    output:
        output_list = ["temp/other_electricity/geothermal/geothermal_power_plant_sites_geojson_{gid}.geojson".format(gid=x) for x in config["country_codes"]]
    params: "Geothermal",
            config["country_codes"]
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/split_other_clean_electricity.py"


rule calculate_geothermal_potential:
    message: "Summarize the installed capacity and average monthly generation of geothermal power for relevant municipals in {wildcards.gid}"
    input:
        power_plant_sites_geojson = "temp/other_electricity/geothermal/geothermal_power_plant_sites_geojson_{gid}.geojson",
        monthly_generation = "temp/other_electricity/other_electricity_monthly_generation/other_electricity_monthly_generation_{gid}.csv",
        municipals = "temp/geoboundaries/{gid}_municipalities.geojson"
    output:
        "temp/geothermal/{gid}/geothermal_potential_{gid}.geojson"
    params:
        'geothermal_MWh'
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/calculate_potential_of_other_electricity.py"


rule preprocess_pop:
    message: "Reproject population density raster to EPSG: 3035"
    input:
        "data/population_density/JRC_1K_POP_2018.tif"
    output:
        "temp/population_density/pop_3035.tif"
    params:
        "EPSG: 3035"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/reproject_raster.py"


rule preprocess_clc:
    message: "Reproject corine land cover raster to EPSG: 3035"
    input:
        "data/corine_land_cover/DATA/U2006_CLC2000_V2020_20u1.tif"
    output:
        "temp/corine_land_cover/clc_3035.tif"
    params:
        "EPSG: 3035"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/reproject_raster.py"


rule preprocess_agriculture:
    message: "Extract the national annual production of involved agricultural biomass types"
    input:
        wheat = "data/biomass/tag00047_page_linear.csv",
        rye = "data/biomass/tag00049_page_linear.csv",
        barley = "data/biomass/tag00051_page_linear.csv",
        oats = "data/biomass/tag00053_page_linear.csv",
        maize = "data/biomass/tag00093_page_linear.csv",
        rape = "data/biomass/tag00110_page_linear.csv",
        sunflower = "data/biomass/tag00120_page_linear.csv",
        rice = "data/biomass/rice-production.csv",
        fruits_and_berries = "data/biomass/ef_lpc_fruit_page_linear.csv",
        vineyards = "data/biomass/ef_lpc_vineyd_page_linear.csv",
        olives = "data/biomass/tag00122_page_linear.csv",
    output:
        "temp/biomass/agriculture_production.csv"
    params: config['country_codes']
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/preprocess_agriculture.py"


rule preprocess_forestry:
    message: "Extract the national annual production of involved forestry biomass types"
    input:
        fuelwood = "data/biomass/for_remov_page_linear_fuelwood.csv",
        roundwood = "data/biomass/for_remov_page_linear_roundwood.csv"
    output:
        "temp/biomass/forestry_production.csv"
    params: config['country_codes']
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/preprocess_forestry.py"


rule preprocess_solid_wastes:
    message: "Calculate the national annual solid wastes generation"
    input:
        solid_wastes = "data/biomass/ten00108_page_linear.csv",
    output:
        "temp/biomass/solid_wastes.csv"
    params: config['country_codes']
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/preprocess_solid_wastes.py"


rule add_population_ratio_to_municipalties:
    message: "Calculate the municipal population ratio to national population in {wildcards.gid}"
    input:
        "temp/geoboundaries/{gid}_municipalities.geojson",
        "temp/population_density/pop_3035.tif"
    output:
        "temp/geoboundaries/{gid}_municipalities_add_pop_ratio.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/add_population_ratio_to_municipalties.py"


rule add_forest_ratio_to_municipalties:
    message: "Calculate the municipal forestry land cover ratio in {wildcards.gid}"
    input:
        "temp/geoboundaries/{gid}_municipalities.geojson",
        "temp/corine_land_cover/clc_3035.tif"
    output:
        "temp/geoboundaries/{gid}_municipalities_add_forest.geojson"
    params: target_cell_values = [23, 24, 25, 27, 28, 29],
            area_name = 'forest_share'
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/add_clc_to_municipalties.py"


rule add_crop_ratio_to_municipalties:
    message: "Calculate the municipal agricultural land cover ratio in {wildcards.gid}"
    input:
        "temp/geoboundaries/{gid}_municipalities.geojson",
        "temp/corine_land_cover/clc_3035.tif",
    output:
        "temp/geoboundaries/{gid}_municipalities_add_crop.geojson"
    params: target_cell_values = [12, 13, 19, 20, 21, 22],
            area_name = 'crop_share'
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/add_clc_to_municipalties.py"


rule add_rice_ratio_to_municipalties:
    message: "Calculate the municipal rice land cover ratio in {wildcards.gid}"
    input:
        "temp/geoboundaries/{gid}_municipalities.geojson",
        "temp/corine_land_cover/clc_3035.tif",
    output:
        "temp/geoboundaries/{gid}_municipalities_add_rice.geojson"
    params: target_cell_values = [14],
            area_name = 'rice_share'
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/add_clc_to_municipalties.py"


rule add_vineyards_ratio_to_municipalties:
    message: "Calculate the municipal vineyards land cover ratio in {wildcards.gid}"
    input:
        "temp/geoboundaries/{gid}_municipalities.geojson",
        "temp/corine_land_cover/clc_3035.tif",
    output:
        "temp/geoboundaries/{gid}_municipalities_add_vineyards.geojson"
    params: target_cell_values = [15],
            area_name = 'vineyards_share'
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/add_clc_to_municipalties.py"


rule add_fruits_ratio_to_municipalties:
    message: "Calculate the municipal fruits and berries land cover ratio in {wildcards.gid}"
    input:
        "temp/geoboundaries/{gid}_municipalities.geojson",
        "temp/corine_land_cover/clc_3035.tif",
    output:
        "temp/geoboundaries/{gid}_municipalities_add_fruits.geojson"
    params: target_cell_values = [16],
            area_name = 'fruits_share'
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/add_clc_to_municipalties.py"


rule add_olives_ratio_to_municipalties:
    message: "Calculate the municipal olives land cover ratio in {wildcards.gid}"
    input:
        "temp/geoboundaries/{gid}_municipalities.geojson",
        "temp/corine_land_cover/clc_3035.tif",
    output:
        "temp/geoboundaries/{gid}_municipalities_add_olives.geojson"
    params: target_cell_values = [17],
            area_name = 'olives_share'
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/add_clc_to_municipalties.py"


def get_country_codes(wildcards):
    return f'{wildcards.gid}'


rule summarize_national_annual_biomass_potential:
    message: "Calculate the annual biomass potential for all 14 types in {wildcards.gid}"
    input:
        "temp/geoboundaries/{gid}_municipalities.geojson",
        "temp/geoboundaries/{gid}_municipalities_add_pop_ratio.geojson",
        "temp/geoboundaries/{gid}_municipalities_add_forest.geojson",
        "temp/geoboundaries/{gid}_municipalities_add_crop.geojson",
        "temp/geoboundaries/{gid}_municipalities_add_rice.geojson",
        "temp/geoboundaries/{gid}_municipalities_add_vineyards.geojson",
        "temp/geoboundaries/{gid}_municipalities_add_fruits.geojson",
        "temp/geoboundaries/{gid}_municipalities_add_olives.geojson",
        "temp/biomass/agriculture_production.csv",
        "temp/biomass/forestry_production.csv",
        "temp/biomass/solid_wastes.csv",
    output:
        "temp/biomass/{gid}/biomass_potential_{gid}.geojson"
    params: get_country_codes
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/summarize_national_annual_biomass_potential.py"


rule summarize_national_annual_energy_supply:
    message: "Summarize the annual supply of all energy supply technologies in {wildcards.gid}"
    input:
        national_municipalities = "temp/geoboundaries/{gid}_municipalities.geojson",
        onshore_wind = "temp/onshore_wind/{gid}/onshore_wind_potential_{gid}.geojson",
        offshore_wind = "temp/offshore_wind/{gid}/offshore_wind_potential_{gid}.geojson",
        open_field_solar_PV = "temp/open_field_solar_PV/{gid}/open_field_solar_PV_potential_{gid}.geojson",
        rooftop_solar_PV = "temp/rooftop_solar_PV/{gid}/rooftop_solar_PV_potential_{gid}.geojson",
        hydro = "temp/hydro/{gid}/hydro_potential_{gid}.geojson",
        nuclear = "temp/nuclear/{gid}/nuclear_potential_{gid}.geojson",
        geothermal = "temp/geothermal/{gid}/geothermal_potential_{gid}.geojson",
        biomass = "temp/biomass/{gid}/biomass_potential_{gid}.geojson"
    output:
        "build/energy_supply/municipal/{gid}_municipal_energy_supply.csv",
        "build/energy_supply/municipal/{gid}_municipal_energy_supply.geojson"
    params:
        capacity_of_onshore_wind_turbine = config['parameters']['capacity_of_onshore_wind_turbine'],
        capacity_of_offshore_wind_turbine = config['parameters']['capacity_of_offshore_wind_turbine'],
        capacity_density_of_open_field_solar_PV = config['parameters']['capacity_density_of_open_field_solar_PV'],
        capacity_density_of_rooftop_solar_PV = config['parameters']['capacity_density_of_rooftop_solar_PV'],
        ratio_of_usable_rooftop = config['parameters']['ratio_of_usable_rooftop'],
        penetration_rate_of_rooftop_solar_PV = config['parameters']['penetration_rate_of_rooftop_solar_PV'],
        biomass_energy_conversion_ratio = config['parameters']['biomass_energy_conversion_ratio']
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/summarize_national_annual_energy_supply.py"