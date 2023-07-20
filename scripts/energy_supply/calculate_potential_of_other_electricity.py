import geopandas as gpd
import pandas as pd


def calculate_potential_of_other_electricity(municipals_geojson, generation_data_csv, power_plant_sites_geojson, wanted_tech, output_geojson):

    # This function calculates the monthly generation of the wanted tech within relevant municipalities
    municipals_geojson = gpd.read_file(municipals_geojson)
    power_plant_sites_geojson = gpd.read_file(power_plant_sites_geojson)
    generation_data_csv = pd.read_csv(generation_data_csv)

    # Attribute power plants to municipalities
    result = gpd.sjoin(municipals_geojson, power_plant_sites_geojson, how='inner', predicate='contains')

    if result.empty == False:
        result = result.loc[:, ['geometry', 'municipal_id', 'capacity_mw']]
        result = result.dissolve(by='municipal_id', aggfunc='sum')

        for n in range(1, 13):
            generation_sum = float(generation_data_csv.loc[n - 1, wanted_tech])
            result[str(n).zfill(2) + "_MWh"] = \
                result['capacity_mw'] / result['capacity_mw'].sum() * generation_sum
        
        result = result[result['01_MWh']!=0]

    else:
        # Beware, if the output is geojson file, then only the geometry column can be kept. other columns need to be added again later
        month_list = [str(n).zfill(2) + "_MWh" for n in range(1, 13)]
        columns = ['geometry', 'municipal_id', 'capacity_mw'] + month_list
        result = gpd.GeoDataFrame(columns=columns, geometry='geometry')    

    result.to_file(output_geojson)


if __name__ == "__main__":
    calculate_potential_of_other_electricity(
        municipals_geojson = snakemake.input.municipals, 
        generation_data_csv = snakemake.input.monthly_generation, 
        power_plant_sites_geojson = snakemake.input.power_plant_sites_geojson, 
        wanted_tech = snakemake.params[0], 
        output_geojson = snakemake.output[0]
        )