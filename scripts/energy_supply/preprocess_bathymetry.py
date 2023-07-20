import rasterio as rio
import numpy as np
from rasterio.warp import calculate_default_transform, reproject, Resampling


def reproject_raster(source_tif, dst_crs, destination_tif):

    # This function reprojects the studied raster to new coordinate reference system
    with rio.open(source_tif) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        with rio.open(destination_tif, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rio.band(src, i),
                    destination=rio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest)


def twice_reproject(source_tif, temperal_tif, destination_tif):

    # The original raster doesn't have specific crs in its metadata, but can be identified as EPSG: 4326
    # Thus, here it is reprojected twice, finally to EPSG: 3035
    reproject_raster(source_tif, 'EPSG:4326', temperal_tif)
    reproject_raster(temperal_tif, 'EPSG:3035', destination_tif)


if __name__ == "__main__":
    twice_reproject(
        source_tif = snakemake.input[0], 
        temperal_tif = snakemake.output[0], 
        destination_tif = snakemake.output[1]
        )