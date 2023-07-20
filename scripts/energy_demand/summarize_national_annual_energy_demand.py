import pandas as pd
import geopandas as gpd


def summarize_national_annual_energy_demand(national_municipalities_geojson, transport_geojson, industry_geojson, building_geojson,
        output_csv, output_geojson):

    # Add transport electricity demand
    transport_gdf = gpd.read_file(transport_geojson)
    transport_gdf = transport_gdf.loc[:, ['municipal_id', 'total_transport_electricity_demand']]
    transport_gdf = transport_gdf.set_index('municipal_id')

    # Add industry electricity and heat demand
    industry_gdf = gpd.read_file(industry_geojson)
    industry_gdf = industry_gdf.loc[:, ['municipal_id', 'total_industry_electricity_demand', 'total_industry_heat_demand']]
    industry_gdf = industry_gdf.set_index('municipal_id')

    # Add building electricity demand
    building_gdf = gpd.read_file(building_geojson)
    building_gdf = building_gdf.loc[:, ['municipal_id', 'total_building_electricity_demand']]
    building_gdf = building_gdf.set_index('municipal_id')

    # Aggregate generations from all demands
    concat_list = [transport_gdf, industry_gdf, building_gdf]
    result_df = pd.concat(concat_list, axis=1).fillna(0).reset_index()

    # Rename
    result_df = result_df.rename(columns={
        'total_transport_electricity_demand': 'transport_electricity_MWh',
        'total_industry_electricity_demand': 'industry_electricity_MWh',
        'total_industry_heat_demand': 'industry_heat_MWh',
        'total_building_electricity_demand': 'building_electricity_MWh'
    })

    # Output to csv and geojson
    result_df.to_csv(output_csv)
    municipalities_gdf = gpd.read_file(national_municipalities_geojson)
    result_gdf = municipalities_gdf.merge(result_df, on='municipal_id')
    result_gdf.to_file(output_geojson)


if __name__ == "__main__":
    summarize_national_annual_energy_demand(
        national_municipalities_geojson = snakemake.input.national_municipalities, 
        transport_geojson = snakemake.input.transport, 
        industry_geojson = snakemake.input.industry, 
        building_geojson = snakemake.input.building,
        output_csv = snakemake.output[0], 
        output_geojson = snakemake.output[1])