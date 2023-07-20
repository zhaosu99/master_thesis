import geokit as gk
import glaes as gl
import pandas as pd


def onshore_wind_turbine_locations(input_path_of_the_examined_area, directory_path_of_prior_datasets, 
                                   turbine_separation_m, output_path_of_turbine_locations_csv):

    # This function generates onshore wind turbine locations for each studied area
    try:
        gl.Priors.loadDirectory(directory_path_of_prior_datasets)
        suitable_areas = gl.ExclusionCalculator(input_path_of_the_examined_area, srs=3035, pixelSize=100, limitOne=False)

        # Geographical exclusions
        suitable_areas.excludePrior("elevation_threshold", value=(1750, None) ) # alpine forest line assumed at 1750 m
        suitable_areas.excludePrior("slope_threshold", value=(15, None) )
        suitable_areas.excludePrior("windspeed_100m_threshold", value=(None, 4) )

        # Natural factors exclusions
        suitable_areas.excludePrior("river_proximity", value=(None, 500) )
        suitable_areas.excludePrior("lake_proximity", value=(None, 1000) )
        suitable_areas.excludePrior("ocean_proximity", value=(None, 1000) )
        suitable_areas.excludePrior("sand_proximity", value=(None, 0) )
        suitable_areas.excludePrior("wetland_proximity", value=(None, 0) )
        suitable_areas.excludePrior("woodland_proximity", value=(None, 300) )
    
        # Human factors exclusions
        suitable_areas.excludePrior("settlement_proximity", value=(None, 2000) )
        suitable_areas.excludePrior("touristic_proximity", value=(None, 2000) )
        suitable_areas.excludePrior("industrial_proximity", value=(None, 2000) )
        suitable_areas.excludePrior("mining_proximity", value=(None, 2000) )
        suitable_areas.excludePrior("camping_proximity", value=(None, 2000) )
        suitable_areas.excludePrior("leisure_proximity", value=(None, 2000) )
        suitable_areas.excludePrior("railway_proximity", value=(None, 500) )
        suitable_areas.excludePrior("roads_main_proximity", value=(None, 500) )
        suitable_areas.excludePrior("roads_secondary_proximity", value=(None, 500) )
        suitable_areas.excludePrior("airport_proximity", value=(None, 5000) )
        suitable_areas.excludePrior("airfield_proximity", value=(None, 5000) )
        suitable_areas.excludePrior("power_line_proximity", value=(None, 300) )

        # Protected areas exclusions
        suitable_areas.excludePrior("protected_park_proximity", value=(None, 1000) )
        suitable_areas.excludePrior("protected_habitat_proximity", value=(None, 0) )
        suitable_areas.excludePrior("protected_bird_proximity", value=(None, 0) )
        suitable_areas.excludePrior("protected_biosphere_proximity", value=(None, 0) )
        suitable_areas.excludePrior("protected_landscape_proximity", value=(None, 0) )
        suitable_areas.excludePrior("protected_natural_monument_proximity", value=(None, 0) )

        # Generate turbine locations
        suitable_areas.distributeItems(separation=turbine_separation_m)
        turbine_coordinates = suitable_areas.itemCoords

        # CRS check
        if suitable_areas.srs.IsSame(gk.srs.EPSG3035) == 1: # this should always be the case in our analysis
            # EPSG3035's unit is meters, and sub-meter precision is not necessary.
            turbine_coordinates = turbine_coordinates.astype(int)

        # Save turbine placement
        (
            pd
            .DataFrame(turbine_coordinates)
            .rename(columns={0: "x_m", 1: "y_m"})
            .to_csv(output_path_of_turbine_locations_csv, index=True, header=True)
        )

    except:
        # Sometimes the input geodataframe is too small for the analysis, 
        # which will raise the error "geokit.core.extent.GeoKitExtentError: Unit size is larger than extent width"
        # For them, an empty csv file is created
        (
        pd
        .DataFrame(columns=['x_m', 'y_m'])
        .to_csv(output_path_of_turbine_locations_csv, index=True, header=True)
        )


if __name__ == "__main__":
    onshore_wind_turbine_locations(
        input_path_of_the_examined_area = snakemake.input.municipal_geojson,
        directory_path_of_prior_datasets = snakemake.params.prior_datasets_path,
        turbine_separation_m = snakemake.params.onshore_turbine_distance,
        output_path_of_turbine_locations_csv = snakemake.output[0]
    )
