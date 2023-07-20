import pandas as pd
import numpy as np
import xarray as xr
from pathlib import Path


def calculate_aggregated_hourly_residential_heat_demand(input_nc_list, sfh_ratios_csv, output_nc):

    # This function calculates the weighted average of two ndarrays
    def get_sfh_ratio(country_code, sfh_ratios_csv):
        sfh_ratios = pd.read_csv(sfh_ratios_csv).loc[lambda x: x.code == country_code, 'sfh_ratio']
        return float(sfh_ratios)

    def get_aggregated_ndarray(ndarray_1, ndarray_2, ratio_1):
        return ndarray_1 * ratio_1 + ndarray_2 * ( 1 - ratio_1)

    sfh_xarray = xr.open_dataset(input_nc_list[0])
    mfh_xarray = xr.open_dataset(input_nc_list[1])

    var_name1 = list(sfh_xarray.data_vars)[0]
    sfh_ndarray = sfh_xarray.data_vars[var_name1].values
    var_name2 = list(mfh_xarray.data_vars)[0]
    mfh_ndarray = mfh_xarray.data_vars[var_name2].values

    country_code = Path(input_nc_list[0]).stem[: 3]
    sfh_ratio = get_sfh_ratio(country_code, sfh_ratios_csv)

    result_ndarray = []
    length = len(sfh_ndarray)
    for n in range(length):
        result_ndarray.append(
            get_aggregated_ndarray(sfh_ndarray[n], mfh_ndarray[n], sfh_ratio)
            )
    
    xr.DataArray(np.array(result_ndarray), coords=sfh_xarray.coords, name='hourly_factor').to_netcdf(output_nc)


if __name__ == "__main__":
    calculate_aggregated_hourly_residential_heat_demand(
        input_nc_list = snakemake.input.water, 
        sfh_ratios_csv = snakemake.input.sfh_ratios,
        output_nc = snakemake.output[0]
        )
    calculate_aggregated_hourly_residential_heat_demand(
        input_nc_list = snakemake.input.space, 
        sfh_ratios_csv = snakemake.input.sfh_ratios,
        output_nc = snakemake.output[1]
        )