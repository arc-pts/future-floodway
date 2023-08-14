from osgeo import gdal
import numpy as np
from os import PathLike

def cov_per_cell(
    multiband_raster: PathLike,
    output_raster: PathLike,
    consider_nodata: bool = False
) -> PathLike:
    if ".tif" in multiband_raster: multiband_raster=multiband_raster.replace(".tif", ".TIF")
    if ".tif" in output_raster: output_raster=output_raster.replace(".tif", ".TIF")
    # Open the multiband raster dataset
    dataset = gdal.Open(multiband_raster, gdal.GA_ReadOnly)
    if dataset is None:
        print("Failed to open the raster dataset.")
        return
    
    # Get the number of bands and dimensions of the raster
    num_bands = dataset.RasterCount
    width = dataset.RasterXSize
    height = dataset.RasterYSize
    
    # Create an empty array to store the pixel values of all bands
    raster_values = np.zeros((num_bands, height, width), dtype=np.float32)
    
    # Read the pixel values of all bands into the array
    for band_index in range(num_bands):
        band = dataset.GetRasterBand(band_index + 1)
        raster_values[band_index, :, :] = band.ReadAsArray()
        nodata_value = band.GetNoDataValue()
        if consider_nodata:
            raster_values[raster_values == nodata_value] = 0
        else:
            raster_values[raster_values == nodata_value] = np.nan

    # Calculate the mean and standard deviation for each cell across all bands
    mean = np.nanmean(raster_values, axis = 0)
    std_dev = np.nanstd(raster_values, axis = 0)

    # Calculate the cov for each cell in the raster array
    cov_array = np.divide(std_dev,mean)
    
    # Create a new single-band raster dataset for the covs
    driver = gdal.GetDriverByName('GTiff')
    cov_dataset = driver.Create(output_raster, width, height, 1, gdal.GDT_Float32)

    # Write the covs to the single-band raster dataset
    cov_band = cov_dataset.GetRasterBand(1)
    cov_band.WriteArray(cov_array[:, :], 0, 0)

    # Set the geotransform and projection of the output raster same as the input
    cov_dataset.SetGeoTransform(dataset.GetGeoTransform())
    cov_dataset.SetProjection(dataset.GetProjection())
    
    # Close the datasets
    dataset = None
    cov_dataset = None
    
    print("Coefficient of variation calculated and saved to:", output_raster)
    return output_raster