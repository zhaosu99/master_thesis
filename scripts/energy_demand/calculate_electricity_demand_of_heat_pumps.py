import numpy as np
import xarray as xr


def calculate_electricity_demand_of_heat_pumps(heat_demand_nc, cop_nc, output_nc):

    # This function calculates the electricity demand ratio of input heat pump type
    heat_demand_xarray = xr.open_dataset(heat_demand_nc)
    cop_xarray = xr.open_dataset(cop_nc)

    var_name1 = list(heat_demand_xarray.data_vars)[0]
    heat_demand_ndarray = heat_demand_xarray.data_vars[var_name1].values
    var_name2 = list(cop_xarray.data_vars)[0]
    cop_ndarray = cop_xarray.data_vars[var_name2].values

    # times 1000000 to increase accuracy, which will be divided in the next step
    sum_2darray = heat_demand_ndarray.sum(axis=0)
    for n in range(len(heat_demand_ndarray)):
        heat_demand_ndarray[n] = heat_demand_ndarray[n] * 1000000 / sum_2darray

    # In the wanted_ndarray, each value represents the electricity/heat value
    # If the ratio = 0.3, then it means: if the annual heat demand in this place is 1 unit, 
    # in this year, it will need 0.3 unit of electricity to cover such demand
    wanted_ndarray = heat_demand_ndarray / cop_ndarray

    xr.DataArray(wanted_ndarray, coords=cop_xarray.coords, name='hourly_factor').to_netcdf(output_nc)


def synthesize(space_heat_demand_nc, water_heat_demand_nc, space_heat_pump_cop_nc, water_heat_pump_cop_nc, ouput_nc_list):

    # Space heat
    calculate_electricity_demand_of_heat_pumps(
        heat_demand_nc = space_heat_demand_nc, 
        cop_nc = space_heat_pump_cop_nc, 
        output_nc = ouput_nc_list[0]
        )

    # Water heat
    calculate_electricity_demand_of_heat_pumps(
        heat_demand_nc = water_heat_demand_nc, 
        cop_nc = water_heat_pump_cop_nc, 
        output_nc = ouput_nc_list[1]
        )


if __name__ == "__main__":
    synthesize(
        space_heat_demand_nc = snakemake.input.residential_space_heat_demand, 
        water_heat_demand_nc = snakemake.input.residential_water_heat_demand, 
        space_heat_pump_cop_nc = snakemake.input.space_heat_pump_cop, 
        water_heat_pump_cop_nc = snakemake.input.water_heat_pump_cop, 
        ouput_nc_list = snakemake.output.residental)

    synthesize(
        space_heat_demand_nc = snakemake.input.commercial_space_heat_demand, 
        water_heat_demand_nc = snakemake.input.commercial_water_heat_demand, 
        space_heat_pump_cop_nc = snakemake.input.space_heat_pump_cop, 
        water_heat_pump_cop_nc = snakemake.input.water_heat_pump_cop, 
        ouput_nc_list = snakemake.output.commercial)