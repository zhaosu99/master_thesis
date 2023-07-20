import geopandas as gpd
import pandas as pd


def preprocess_offshore_oil_and_gas_geojsons(geojson_to_process, output_geojson):

    # This function simplifies and selects qualifed offshore oil and gas wells
    oil_and_gas = gpd.read_file(geojson_to_process)
    oil_and_gas = oil_and_gas.loc[:, ['geometry', 'COUNTRY', 'STATUS']]
    oil_and_gas = oil_and_gas[oil_and_gas['STATUS'].isin(['Operational', 'Under construction'])]
    oil_and_gas.to_file(output_geojson)


if __name__ == "__main__":
    preprocess_offshore_oil_and_gas_geojsons(
        geojson_to_process = snakemake.input[0], 
        output_geojson = snakemake.output[0]
        )
