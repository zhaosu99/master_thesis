import geopandas as gpd


def preprocess_offshore_pipelines(geojson_to_process, output_geojson):

    # This function simplifies and selects qualified offshore pipelines
    pipelines = gpd.read_file(geojson_to_process)
    pipelines = pipelines.loc[:, ['geometry', 'COUNTRY', 'STATUS']]
    pipelines = pipelines[pipelines['STATUS'] != 'Not in use']
    pipelines.to_file(output_geojson)


if __name__ == "__main__":
    preprocess_offshore_pipelines(
        geojson_to_process = snakemake.input[0], 
        output_geojson = snakemake.output[0]
        )
