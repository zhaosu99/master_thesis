import geopandas as gpd


def preprocess_geoboundaries(input_geojson, output_geojson):

    # Select geometry and shapeID (renamed to municipal_id) columns for further analysis
    # For municipal id, this analysis adopts the standardized form: {three-letter country code}-{ADM level}-{index}, like AUT-ADM4-1
    new_col_name = 'municipal_id'
    gdf = gpd.read_file(input_geojson).to_crs(3035)
    gdf = gdf.loc[:, ['geometry', 'shapeID']]
    gdf = gdf.rename(columns={'shapeID': new_col_name})
    gdf[new_col_name] = gdf[new_col_name].apply(lambda x: x[0: 9] + x[16: ])
    gdf.to_file(output_geojson)


if __name__ == "__main__":
    preprocess_geoboundaries(
        input_geojson = snakemake.input[0], 
        output_geojson = snakemake.output[0]
        )