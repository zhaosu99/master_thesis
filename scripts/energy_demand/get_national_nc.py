import xarray as xr
import geopandas as gpd


def select_national_nc(input_gdf, input_nc, output_nc):

    def get_geobounds(input_gdf):

        min_x, min_y, max_x, max_y = input_gdf.copy().to_crs(4326).total_bounds
        def get_smaller_quarters(x):
            return x // 0.25 * 0.25
        def get_larger_quarters(x):
            return (x // 0.25 + 1) * 0.25
        min_lon = get_smaller_quarters(min_x)
        min_lat = get_smaller_quarters(min_y)
        max_lon = get_larger_quarters(max_x)
        max_lat = get_larger_quarters(max_y)

        return (min_lon, min_lat, max_lon, max_lat)

    min_lon, min_lat, max_lon, max_lat = get_geobounds(input_gdf)
    input_dataset = xr.open_dataset(input_nc)
    # Beware the sequences of latitudes and longitudes are different
    output_dataset = input_dataset.copy().sel(latitude=slice(max_lat, min_lat), longitude=slice(min_lon, max_lon))

    output_dataset.to_netcdf(output_nc)


def batch_process(input_geojson, input_nc_list, output_nc_list):

    input_gdf = gpd.read_file(input_geojson)
    for n in range(len(input_nc_list)):
        select_national_nc(input_gdf, input_nc_list[n], output_nc_list[n])


if __name__ == "__main__":
    batch_process(
        input_geojson = snakemake.input.geoboundaries_geojson, 
        input_nc_list = snakemake.input.raw_nc_list, 
        output_nc_list = snakemake.output.output_nc_list
        )
