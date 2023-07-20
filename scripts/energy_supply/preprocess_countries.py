import geopandas as gpd
from pathlib import Path


def preprocess_geoboundaries(input_geojson, output_geojson):

    # Select geometry and shapeID (renamed to municipal_id) columns for further analysis
    # For regional id, this analysis adopts the standardized form: {three-letter country code}-{ADM level}-{index}, like AUT-ADM4-1


    country_gdf = gpd.read_file(input_geojson).to_crs(3035)

    # Exclude French and Portuguese oversea territories since they are not in the study scope
    if Path(input_geojson).stem[: 3] == "FRA":
        national_regions = national_regions.loc[: 12, :]
    if Path(input_geojson).stem[: 3] == "PRT":
        national_regions = national_regions.iloc[2: , :]

    national_regions = national_regions.loc[:, ['geometry', 'shapeID']]
    national_regions = national_regions.rename(columns={'shapeID': 'regional_id'})
    if Path(input_geojson).stem[: 3] == "GBR":
        national_regions['regional_id'] = [f'GBR-ADM1-{x}' for x in range(1, 5)]
    else:
        national_regions['regional_id'] = national_regions['regional_id'].apply(lambda x: x[0: 9] + x[16: ])
    national_regions.to_file(output_geojson)


if __name__ == "__main__":
    preprocess_geoboundaries(
        input_geojson = snakemake.input[0], 
        output_geojson = snakemake.output[0]
        )