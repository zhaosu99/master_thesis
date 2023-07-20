import glaes as gl


def open_field_solar_PV_areas(input_path_of_the_examined_area, directory_path_of_prior_datasets, onshore_wind_turbine_path, output_path_of_suitable_areas_tif):

    #This function generates suitable areas for open field solar PV plants

    gl.Priors.loadDirectory(directory_path_of_prior_datasets)

    suitable_areas = gl.ExclusionCalculator(input_path_of_the_examined_area, srs=3035, pixelSize=100, limitOne=False)
    
    # Geographical exclusions
    suitable_areas.excludePrior("elevation_threshold", value=(1750, None) ) # alpine forest line assumed at 1750 m
    suitable_areas.excludePrior("slope_threshold", value=(15, None) )
    suitable_areas.excludePrior("dni_threshold", value=(None, 1) )
    
    # Topographical factors exclusions
    suitable_areas.excludePrior("river_proximity", value=(None, 300) )
    suitable_areas.excludePrior("lake_proximity", value=(None, 300) )
    suitable_areas.excludePrior("ocean_proximity", value=(None, 1000) )
    suitable_areas.excludePrior("wetland_proximity", value=(None, 300) )
    suitable_areas.excludePrior("woodland_proximity", value=(None, 300) )

    # Agricultural factors exclusions
    suitable_areas.excludePrior("agriculture_permanent_crop_proximity", value=(None, 300) )
    suitable_areas.excludePrior("agriculture_heterogeneous_proximity", value=(None, 300) )
    suitable_areas.excludePrior("agriculture_arable_proximity", value=(None, 300) )
    
    # Human activities factors exclusions
    suitable_areas.excludePrior("settlement_proximity", value=(None, 1000) )
    suitable_areas.excludePrior("touristic_proximity", value=(None, 2000) )
    suitable_areas.excludePrior("industrial_proximity", value=(None, 300) )
    suitable_areas.excludePrior("mining_proximity", value=(None, 300) )
    suitable_areas.excludePrior("camping_proximity", value=(None, 300) )
    suitable_areas.excludePrior("leisure_proximity", value=(None, 300) )
    suitable_areas.excludePrior("railway_proximity", value=(None, 300) )
    suitable_areas.excludePrior("roads_main_proximity", value=(None, 300) )
    suitable_areas.excludePrior("roads_secondary_proximity", value=(None, 300) )
    suitable_areas.excludePrior("airport_proximity", value=(None, 5000) )
    suitable_areas.excludePrior("airfield_proximity", value=(None, 5000) )
    suitable_areas.excludePrior("power_line_proximity", value=(None, 300) )
      
    # Protected areas exclusions
    suitable_areas.excludePrior("protected_park_proximity", value=(None, 1000) )
    suitable_areas.excludePrior("protected_habitat_proximity", value=0 )
    suitable_areas.excludePrior("protected_bird_proximity", value=0 )
    suitable_areas.excludePrior("protected_biosphere_proximity", value=0 )
    suitable_areas.excludePrior("protected_landscape_proximity", value=0 )
    suitable_areas.excludePrior("protected_natural_monument_proximity", value=0 )   

    # Onshore wind sites exclusion
    suitable_areas.excludeVectorType(onshore_wind_turbine_path, buffer=100)

    suitable_areas.save(output_path_of_suitable_areas_tif, overwrite=True)

            
if __name__ == "__main__":
    open_field_solar_PV_areas(
        input_path_of_the_examined_area = snakemake.input.region, 
        directory_path_of_prior_datasets = snakemake.params.prior_datasets_path,
        onshore_wind_turbine_path = snakemake.input.onshore_turbine_locations_geojson, 
        output_path_of_suitable_areas_tif = snakemake.output[0]
        )

