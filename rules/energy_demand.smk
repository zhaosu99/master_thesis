from pathlib import Path 


rule add_population_amount_to_units:
    message: "Add the population to each municipality in {wildcards.gid}"
    input:
        "temp/geoboundaries/{gid}_municipalities.geojson",
        "temp/population_density/pop_3035.tif"
    output:
        "temp/geoboundaries/{gid}_municipalities_add_pop_amount.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/add_population_amount_to_units.py"


rule add_railway_length_to_municipalties:
    message: "Add the railway pixel amount to each municipality {wildcards.gid}"
    input:
        "temp/geoboundaries/{gid}_municipalities.geojson",
        "data/excludePrior/railway_proximity.1501389965.tif"
    output:
        "temp/geoboundaries/{gid}_municipalities_add_railway_length.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/add_railway_length_to_municipalties.py"


def func_jrc_transport_xlsx_list(wildcards):
    return list(Path("data/JRC-IDEES").glob("JRC-IDEES-2015_Transport_??.xlsx"))


rule get_regression_of_other_transports:
    message: "Calculate the regression function of specific energy demand and population"
    input:
        pop_csv = "data/others/tps00001_page_linear.csv",
        jrc_xlsx_list = func_jrc_transport_xlsx_list
    output:
        "temp/transport/regression_of_other_transports.csv"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/regression_of_other_transports.py"


rule calculate_transport_electricity_demand:
    message: "Calculate the transport electricity demand in {wildcards.gid}"
    input:
        pop_amount_geojson = "temp/geoboundaries/{gid}_municipalities_add_pop_amount.geojson",
        regression_csv = "temp/transport/regression_of_other_transports.csv",
        jrc_xlsx_list = func_jrc_transport_xlsx_list
    output:
        "temp/transport/{gid}/transport_electricity_demand_{gid}.geojson"
    params: config['country_codes']
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/calculate_transport_electricity_demand.py"


def func_jrc_industry_xlsx_list(wildcards):
    return list(Path("data/JRC-IDEES").glob("JRC-IDEES-2015_Industry_??.xlsx"))


rule calculate_industrial_energy_demand:
    message: "Calculate the industrial energy demand in all countries in the scope"
    input:
        municipalities = expand("temp/geoboundaries/{gid}_municipalities_add_pop_amount.geojson", gid=config["country_codes"]),
        ETS_data = "data/others/industrial emission database.xlsx",
        jrc_xlsx_list = func_jrc_industry_xlsx_list,
        backup_csv = "data/others/ten00124_linear.csv"
    output:
        expand("temp/industry/{gid}/industry_electricity_and_heat_demand_{gid}.geojson", gid=config["country_codes"])
    params:
        country_codes = config["country_codes"]
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/calculate_industry_energy_demand.py"


rule download_ERA:
    message: "Download temperature and wind speed nc data"
    output:
        wind = "data/ERA/wind.nc",
        temp = ["data/ERA/temp_{}.nc".format(x) for x in [str(i) for i in range(2000, 2019)]]
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/download_ERA.py"


rule preprocess_ERA_temperature:
    message: "Calculate annual average of temperature of existing simulation points"
    input:
        ["data/ERA/temp_{}.nc".format(x) for x in [str(i) for i in range(2000, 2019)]]
    output:
        "temp/ERA/temperature.nc"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/preprocess_ERA_temperature.py"


rule preprocess_ERA_wind:
    message: "Calculate annual average of wind speed of existing simulation points"
    input:
        "data/ERA/wind.nc"
    output:
        "temp/ERA/daily_wind.nc"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/preprocess_ERA_wind.py"


rule calculate_COP:
    message: "Calculate hourly coefficient of performance for radiator and water heat pumps"
    input:
        "temp/ERA/temperature.nc"
    output:
        "temp/ERA/radiator_cop.nc",
        "temp/ERA/water_cop.nc"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/calculate_COP.py"


rule calculate_daily_temperature:
    message: "Calculate the daily reference temperature of existing simulation points"
    input:
        "temp/ERA/temperature.nc"
    output:
        "temp/ERA/daily_reference_temperature.nc"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/calculate_daily_temperature.py"


rule get_national_nc:
    message: "Mask original netcdf file to the territory of {wildcards.gid}"
    input:
        geoboundaries_geojson = "temp/geoboundaries/{gid}_regions.geojson",
        raw_nc_list = ["temp/ERA/daily_reference_temperature.nc",
        "temp/ERA/daily_wind.nc",
        "temp/ERA/radiator_cop.nc",
        "temp/ERA/water_cop.nc"]
    output:
        output_nc_list = ["temp/ERA/{gid}/{gid}_daily_reference_temperature.nc",
        "temp/ERA/{gid}/{gid}_daily_wind.nc",
        "temp/ERA/{gid}/{gid}_radiator_cop.nc",
        "temp/ERA/{gid}/{gid}_water_cop.nc"]
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/get_national_nc.py"



rule calculate_daily_heat_demand:
    message: "Calculate the daily heat demand of single-family, multi-family and commercial house in windy and normal weather in {wildcards.gid}"
    input:
        "temp/ERA/{gid}/{gid}_daily_reference_temperature.nc",
        "data/when2heat/bgw_bdew/daily_demand.csv",     
        "data/others/heating_threshold_temperatures.csv"
    output:
        total = ["temp/ERA/{gid}/{gid}_daily_total_heat_demand_sfh.nc",
            "temp/ERA/{gid}/{gid}_daily_total_heat_demand_mfh.nc",
            "temp/ERA/{gid}/{gid}_daily_total_heat_demand_com.nc",
            "temp/ERA/{gid}/{gid}_daily_total_heat_demand_sfh_windy.nc",
            "temp/ERA/{gid}/{gid}_daily_total_heat_demand_mfh_windy.nc",
            "temp/ERA/{gid}/{gid}_daily_total_heat_demand_com_windy.nc"],
        water = ["temp/ERA/{gid}/{gid}_daily_water_heat_demand_sfh.nc",
            "temp/ERA/{gid}/{gid}_daily_water_heat_demand_mfh.nc",
            "temp/ERA/{gid}/{gid}_daily_water_heat_demand_com.nc",
            "temp/ERA/{gid}/{gid}_daily_water_heat_demand_sfh_windy.nc",
            "temp/ERA/{gid}/{gid}_daily_water_heat_demand_mfh_windy.nc",
            "temp/ERA/{gid}/{gid}_daily_water_heat_demand_com_windy.nc"],      
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/calculate_daily_heat_demand.py"


rule merge_windy_and_normal:
    message: "Merge normal and windy daily water and total heat demand series acoording to wind speed in {wildcards.gid}"
    input:
        wind = "temp/ERA/{gid}/{gid}_daily_wind.nc",
        water = ["temp/ERA/{gid}/{gid}_daily_water_heat_demand_sfh.nc",
            "temp/ERA/{gid}/{gid}_daily_water_heat_demand_mfh.nc",
            "temp/ERA/{gid}/{gid}_daily_water_heat_demand_com.nc",
            "temp/ERA/{gid}/{gid}_daily_water_heat_demand_sfh_windy.nc",
            "temp/ERA/{gid}/{gid}_daily_water_heat_demand_mfh_windy.nc",
            "temp/ERA/{gid}/{gid}_daily_water_heat_demand_com_windy.nc"],
        total = ["temp/ERA/{gid}/{gid}_daily_total_heat_demand_sfh.nc",
            "temp/ERA/{gid}/{gid}_daily_total_heat_demand_mfh.nc",
            "temp/ERA/{gid}/{gid}_daily_total_heat_demand_com.nc",
            "temp/ERA/{gid}/{gid}_daily_total_heat_demand_sfh_windy.nc",
            "temp/ERA/{gid}/{gid}_daily_total_heat_demand_mfh_windy.nc",
            "temp/ERA/{gid}/{gid}_daily_total_heat_demand_com_windy.nc"]
    output:
        water = ["temp/ERA/{gid}/{gid}_aggregated_daily_water_heat_demand_sfh.nc",
            "temp/ERA/{gid}/{gid}_aggregated_daily_water_heat_demand_mfh.nc",
            "temp/ERA/{gid}/{gid}_aggregated_daily_water_heat_demand_com.nc"],
        total = ["temp/ERA/{gid}/{gid}_aggregated_daily_total_heat_demand_sfh.nc",
            "temp/ERA/{gid}/{gid}_aggregated_daily_total_heat_demand_mfh.nc",
            "temp/ERA/{gid}/{gid}_aggregated_daily_total_heat_demand_com.nc"]
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/merge_windy_and_normal.py"


rule calculate_hourly_heat_demand:
    message: "Calculate hourly water and total heat demand series of single-family, multi-family and commercial houses in {wildcards.gid}"
    input:
        temperature = "temp/ERA/{gid}/{gid}_daily_reference_temperature.nc",
        water = ["temp/ERA/{gid}/{gid}_aggregated_daily_water_heat_demand_com.nc",
            "temp/ERA/{gid}/{gid}_aggregated_daily_water_heat_demand_sfh.nc",
            "temp/ERA/{gid}/{gid}_aggregated_daily_water_heat_demand_mfh.nc"],
        total = ["temp/ERA/{gid}/{gid}_aggregated_daily_total_heat_demand_com.nc",
            "temp/ERA/{gid}/{gid}_aggregated_daily_total_heat_demand_sfh.nc",
            "temp/ERA/{gid}/{gid}_aggregated_daily_total_heat_demand_mfh.nc"],
        hourly_factor = ["data/when2heat/bgw_bdew/hourly_factors_COM.csv",
            "data/when2heat/bgw_bdew/hourly_factors_SFH.csv",
            "data/when2heat/bgw_bdew/hourly_factors_MFH.csv"]
    output:
        water = ["temp/ERA/{gid}/{gid}_hourly_water_heat_demand_com.nc",
            "temp/ERA/{gid}/{gid}_hourly_water_heat_demand_sfh.nc",
            "temp/ERA/{gid}/{gid}_hourly_water_heat_demand_mfh.nc"],
        total = ["temp/ERA/{gid}/{gid}_hourly_total_heat_demand_com.nc",
            "temp/ERA/{gid}/{gid}_hourly_total_heat_demand_sfh.nc",
            "temp/ERA/{gid}/{gid}_hourly_total_heat_demand_mfh.nc"]
    params: building_types = ["COM", "SFH", "MFH"]
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/calculate_hourly_heat_demand.py"


rule calculate_hourly_space_heat_demand:
    message: "Calculate hourly space heat demand series of single-family, multi-family and commercial houses in {wildcards.gid}"
    input:
        water = ["temp/ERA/{gid}/{gid}_hourly_water_heat_demand_sfh.nc",
            "temp/ERA/{gid}/{gid}_hourly_water_heat_demand_mfh.nc",
            "temp/ERA/{gid}/{gid}_hourly_water_heat_demand_com.nc"],
        total = ["temp/ERA/{gid}/{gid}_hourly_total_heat_demand_sfh.nc",
            "temp/ERA/{gid}/{gid}_hourly_total_heat_demand_mfh.nc",
            "temp/ERA/{gid}/{gid}_hourly_total_heat_demand_com.nc"]
    output:
        space = ["temp/ERA/{gid}/{gid}_hourly_space_heat_demand_sfh.nc",
            "temp/ERA/{gid}/{gid}_hourly_space_heat_demand_mfh.nc",
            "temp/ERA/{gid}/{gid}_hourly_space_heat_demand_com.nc"]
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/calculate_hourly_space_heat_demand.py"


rule calculate_sfh_ratios:
    message: "Calculate the single-family house ratio in each country"
    input:
        "data/others/ilc_lvho01_linear.csv"
    output:
        "temp/others/sfh_ratios.csv"
    params: config['country_codes']
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/calculate_sfh_ratios.py"


rule calculate_aggregated_hourly_residential_heat_demand:
    message: "Aggregate single-family and multi-family heat demand to space heat demand in {wildcards.gid}"
    input:
        water = ["temp/ERA/{gid}/{gid}_hourly_water_heat_demand_sfh.nc",
            "temp/ERA/{gid}/{gid}_hourly_water_heat_demand_mfh.nc"],
        space = ["temp/ERA/{gid}/{gid}_hourly_space_heat_demand_sfh.nc",
            "temp/ERA/{gid}/{gid}_hourly_space_heat_demand_mfh.nc"],
        sfh_ratios = "temp/others/sfh_ratios.csv"
    output:
        "temp/ERA/{gid}/{gid}_hourly_water_heat_demand_res.nc",
        "temp/ERA/{gid}/{gid}_hourly_space_heat_demand_res.nc"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/calculate_aggregated_hourly_residential_heat_demand.py"


rule calculate_electricity_demand_of_heat_pumps:
    message: "Calculate electricity demand of three types of heat pumps in residential and commercial buildings in {wildcards.gid}"
    input:
        residential_space_heat_demand = "temp/ERA/{gid}/{gid}_hourly_space_heat_demand_res.nc",    
        commercial_space_heat_demand = "temp/ERA/{gid}/{gid}_hourly_space_heat_demand_com.nc",
        space_heat_pump_cop = "temp/ERA/{gid}/{gid}_radiator_cop.nc",
        residential_water_heat_demand = "temp/ERA/{gid}/{gid}_hourly_water_heat_demand_res.nc",
        commercial_water_heat_demand = "temp/ERA/{gid}/{gid}_hourly_water_heat_demand_com.nc",
        water_heat_pump_cop = "temp/ERA/{gid}/{gid}_water_cop.nc"
    output:
        residental = ["temp/ERA/final_result/{gid}_residental_radiator_heat_pump_electricity_demand.nc", 
            "temp/ERA/final_result/{gid}_residental_water_heat_pump_electricity_demand.nc"],
        commercial = ["temp/ERA/final_result/{gid}_commercial_radiator_heat_pump_electricity_demand.nc", 
            "temp/ERA/final_result/{gid}_commercial_water_heat_pump_electricity_demand.nc"]
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/calculate_electricity_demand_of_heat_pumps.py"


rule get_annual_heat_electricity_ratio:
    message: "Calculate annual average heat-elctricity conversion ratio in each existing points in {wildcards.gid}"
    input:
        ["temp/ERA/final_result/{gid}_residental_radiator_heat_pump_electricity_demand.nc", 
        "temp/ERA/final_result/{gid}_residental_water_heat_pump_electricity_demand.nc",
        "temp/ERA/final_result/{gid}_commercial_radiator_heat_pump_electricity_demand.nc", 
        "temp/ERA/final_result/{gid}_commercial_water_heat_pump_electricity_demand.nc"]
    output:
        "temp/ERA/final_result/{gid}_heat_electricity_ratio.nc",
        "temp/ERA/final_result/{gid}_heat_electricity_ratio.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/get_annual_heat_electricity_ratio.py"


def func_jrc_residential_xlsx_list(wildcards):
    return list(Path("data/JRC-IDEES").glob("JRC-IDEES-2015_Residential_??.xlsx"))


def func_jrc_commercial_xlsx_list(wildcards):
    return list(Path("data/JRC-IDEES").glob("JRC-IDEES-2015_Tertiary_??.xlsx"))


rule get_regression_of_building_energy_demand:
    message: "Calculate regression functions of population and different types of energy demands"
    input:
        pop_csv = "data/others/tps00001_page_linear.csv",
        res_jrc_xlsx_list = func_jrc_residential_xlsx_list,
        com_jrc_xlsx_list = func_jrc_commercial_xlsx_list
    output:
        "temp/others/regression_of_building_energy_demand_to_pop.csv"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/regression_of_building_energy_demand.py"


rule calculate_building_energy_demand:
    message: "Calcuate electricity demand in building sector in {wildcards.gid}"
    input:
        pop_ratio_geojson = "temp/geoboundaries/{gid}_municipalities_add_pop_ratio.geojson",
        heat_electricity_ratio_geojson = "temp/ERA/final_result/{gid}_heat_electricity_ratio.geojson",
        res_jrc_xlsx_list = func_jrc_residential_xlsx_list,
        com_jrc_xlsx_list = func_jrc_commercial_xlsx_list,
        regression_csv = "temp/others/regression_of_building_energy_demand_to_pop.csv",
        pop_csv = "data/others/tps00001_page_linear.csv"
    params: config['country_codes']
    output:
        "temp/building/{gid}/building_electricity_demand_{gid}.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/calculate_building_electricity_demand.py"


rule get_continental_electricity_demand_of_transport:
    message: "Calculate continental electricity demand in transport sector"
    input:
        expand("temp/transport/{gid}/transport_electricity_demand_{gid}.geojson", gid=config["country_codes"])
    output:
        "temp/transport/transport_electricity_demand_continental.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/concat_geojson_list.py"


rule get_continental_energy_demand_of_industry:
    message: "Calculate continental electricity and heat demand in industry sector"
    input:
        expand("temp/industry/{gid}/industry_electricity_and_heat_demand_{gid}.geojson", gid=config["country_codes"])
    output:
        "temp/industry/industry_electricity_and_heat_demand_continental.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/concat_geojson_list.py"


rule get_continental_electricity_demand_of_building:
    message: "Calculate continental electricity demand in building sector"
    input: 
        expand("temp/building/{gid}/building_electricity_demand_{gid}.geojson", gid=config["country_codes"])
    output:
        "temp/building/building_electricity_demand_continental.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/concat_geojson_list.py"


rule summarize_national_annual_energy_demand:
    message: "Summarize the annual demand of all energy supply technologies in {wildcards.gid}"
    input:
        national_municipalities = "temp/geoboundaries/{gid}_municipalities.geojson",
        transport = "temp/transport/{gid}/transport_electricity_demand_{gid}.geojson",
        industry = "temp/industry/{gid}/industry_electricity_and_heat_demand_{gid}.geojson",
        building = "temp/building/{gid}/building_electricity_demand_{gid}.geojson"
    output:
        "build/energy_demand/municipal/{gid}_municipal_energy_demand.csv",
        "build/energy_demand/municipal/{gid}_municipal_energy_demand.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/summarize_national_annual_energy_demand.py"


rule calculate_energy_inventory_for_regions:
    message: "Calculate the energy demand and supply for each region in {wildcards.gid}"
    input:
        "temp/geoboundaries/{gid}_regions.geojson",
        "build/energy_supply/municipal/{gid}_municipal_energy_supply.geojson",
        "build/energy_demand/municipal/{gid}_municipal_energy_demand.geojson"
    output:
        "build/energy_supply/regional/{gid}_regional_energy_supply.geojson",
        "build/energy_demand/regional/{gid}_regional_energy_demand.geojson",
        "build/energy_supply/regional/{gid}_regional_energy_supply.csv",
        "build/energy_demand/regional/{gid}_regional_energy_demand.csv"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/calculate_energy_inventory_for_regions.py"


rule concat_municipal_demand_and_supply:
    message: "Concatenate the energy inventory for each municipality in {wildcards.gid}"
    input:
        "temp/geoboundaries/{gid}_municipalities.geojson",
        "build/energy_demand/municipal/{gid}_municipal_energy_demand.csv",
        "build/energy_supply/municipal/{gid}_municipal_energy_supply.csv"
    output:
        "build/result/municipal/{gid}_municipal_energy_inventory.csv"
    params: "municipality"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/concat_demand_and_supply.py"


rule concat_regional_demand_and_supply:
    message: "Concatenate the energy inventory for each region in {wildcards.gid}"
    input:
        "temp/geoboundaries/{gid}_regions.geojson",
        "build/energy_supply/regional/{gid}_regional_energy_supply.csv",
        "build/energy_demand/regional/{gid}_regional_energy_demand.csv"
    output:
        "build/result/regional/{gid}_regional_energy_inventory.csv"
    params: 'region'
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/concat_demand_and_supply.py"


rule concat_continental_municipal_energy_inventory:
    input:
        expand("build/result/municipal/{gid}_municipal_energy_inventory.csv", gid=config['country_codes'])
    output:
        "build/result/municipal/continental_municipal_energy_inventory.csv"
    params: "municipal_id"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/concat_csv_list.py"


rule concat_continental_regional_energy_inventory:
    input:
        expand("build/result/regional/{gid}_regional_energy_inventory.csv", gid=config['country_codes'])
    output:
        "build/result/regional/continental_regional_energy_inventory.csv"
    params: "regional_id"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_supply/concat_csv_list.py"


