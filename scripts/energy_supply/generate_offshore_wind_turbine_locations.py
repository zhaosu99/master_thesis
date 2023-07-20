import geokit as gk
import glaes as gl
import pandas as pd


def offshore_wind_turbine_locations(input_path_of_the_examined_area, maritime_military_zone_geojson, offshore_oil_and_gas_geojson, offshore_pipeline_geojson, 
                        protected_area_geojson, bathymetry_tif_path, marine_traffic_tif_path, shoreline_geojson, turbine_separation_m, output_path_of_turbine_locations_csv):

    # This function generates offshore wind turbine locations for studied area
    # In the 'except' part, an empty csv is created for countries who don't have EEZ
    try:
        suitable_areas = gl.ExclusionCalculator(input_path_of_the_examined_area, srs=3035, pixelSize=100, limitOne=False)

        # Exclude proteced areas
        suitable_areas.excludeVectorType(protected_area_geojson, buffer = 3000)

        # Exclude military zones
        suitable_areas.excludeVectorType(maritime_military_zone_geojson, buffer = 3000)

        # Exclude offshore oil and gas wells
        suitable_areas.excludeVectorType(offshore_oil_and_gas_geojson, buffer = 3000)

        # Exclude offshore pipelines
        suitable_areas.excludeVectorType(offshore_pipeline_geojson, buffer = 3000)

        # Exclude sea areas within 10 km to shorelines
        suitable_areas.excludeVectorType(shoreline_geojson, buffer = 10000)

        # Exclude overly deep sea areas
        suitable_areas.excludeRasterType(bathymetry_tif_path, value = (None, -50))

        # Exclude overly busy maritime routes
        #suitable_areas.excludeRasterType(marine_traffic_tif_path, value = (20, None))

        # Generate offshore wind turbines
        suitable_areas.distributeItems(separation = turbine_separation_m)
        turbine_coordinates = suitable_areas.itemCoords

        # CRS check
        if suitable_areas.srs.IsSame(gk.srs.EPSG3035) == 1: # this should always be the case in our analysis
            # EPSG3035's unit is meters, and sub-meter precision is not necessary.
            turbine_coordinates = turbine_coordinates.astype(int)

        # Save turbine placement
        (
            pd.DataFrame(turbine_coordinates)
            .rename(columns={0: "x_m", 1: "y_m"})
            .to_csv(output_path_of_turbine_locations_csv, index=True, header=True)
        )
    except:      
        (
        pd.DataFrame(columns=['x_m', 'y_m'])
        .to_csv(output_path_of_turbine_locations_csv, index=True, header=True)
        )
      

if __name__ == "__main__":
    offshore_wind_turbine_locations(
        input_path_of_the_examined_area = snakemake.input.eez_by_country, 
        maritime_military_zone_geojson = snakemake.input.offshore_military_areas, 
        offshore_oil_and_gas_geojson = snakemake.input.offshore_oil_and_gas, 
        offshore_pipeline_geojson = snakemake.input.offshore_pipelines, 
        protected_area_geojson = snakemake.input.offshore_protected_areas_by_country, 
        bathymetry_tif_path = snakemake.input.bathymetry_by_country, 
        marine_traffic_tif_path = snakemake.input.offshore_traffic, 
        shoreline_geojson = snakemake.input.shoreline,
        turbine_separation_m = snakemake.params.offshore_turbine_distance, 
        output_path_of_turbine_locations_csv = snakemake.output[0]
        )

