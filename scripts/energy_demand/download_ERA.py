import cdsapi
import datetime


def weather(data_package, variable, dates, product_type, file):

    c = cdsapi.Client()

    params = {
        'format': 'netcdf',
        'variable': variable,
        "year": dates["year"],
        "month": dates["month"],
        "time": dates["time"],
        'product_type': product_type,
        'area': [75, -15, 30, 40]
    }

    if (variable == '2m_temperature') | (variable == 'soil_temperature_level_1'):
        params["day"] = ["%.2d" % day for day in range(1, 32)]
    c.retrieve(data_package, params, file)


def wind(output_path):

    # Select all months from 2000 to 2018 by the date of the first day of the month
    data_package = 'reanalysis-era5-single-levels-monthly-means'
    variable = "10m_wind_speed"
    product_type = 'monthly_averaged_reanalysis'
    dates = {
        'year': [str(year) for year in range(2000, 2019)],
        'month': ["%.2d" % month for month in range(1, 13)],
        'time': [datetime.time(i).strftime('%H:%M') for i in range(1)]
        }

    # Call the general weather download function with wind specific parameters
    weather(data_package, variable, dates, product_type, output_path)


def temperatures(output_path_list):

    for n in range(19):

        #Select period
        data_package = 'reanalysis-era5-single-levels'
        variable = '2m_temperature'
        product_type = 'reanalysis'
        dates = {
                'year': ["%.2d" % y for y in range(2000, 2019)][n],
                'month': ["%.2d" % month for month in range(1, 13)],
                'day':  ["%.2d" % day for day in range(1, 32)],
                'time': [datetime.time(i).strftime('%H:%M') for i in range(24)]
                }

        # Call the general weather download function with temperature specific parameters
        weather(data_package, variable, dates, product_type, output_path_list[n])
        

def all(wind_output_path, temp_output_path_list):
    wind(wind_output_path)
    temperatures(temp_output_path_list)


if __name__ == "__main__":
    all(
        wind_output_path = snakemake.output.wind, 
        temp_output_path_list = snakemake.output.temp
        )