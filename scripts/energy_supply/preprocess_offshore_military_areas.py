import geopandas as gpd


def preprocess_offshore_military_areas(geojson_to_process, output_geojson):

    # This function simplifies and picks out qualified marine military areas
    military_areas = gpd.read_file(geojson_to_process)
    military_areas = military_areas.loc[:, ['geometry', 'COUNTRY', 'STATUS']]
    military_areas = military_areas[military_areas['STATUS'] != 'Deactivated']
    military_areas.to_file(output_geojson)


if __name__ == "__main__":
    preprocess_offshore_military_areas(
        geojson_to_process = snakemake.input[0], 
        output_geojson = snakemake.output[0]
        )
