from osgeo import gdal
import numpy as np

# Calculates the z-score per raster cell for a user-defined band in a multi-band raster. Returns 
#   a single band raster of calculated z-scores.
# Inputs: 
#   raster_path - the path to the multi-band raster
#   raster_band - the band that should have the z-scores calculated for
#   output_path - the single-band, z-score raster
#   nodata_value - the value representing missing data, default = -9999

def zscore_per_cell(raster_path, raster_band, output_path, nodata_value = -9999):
    # Open the multiband raster dataset
    dataset = gdal.Open(raster_path, gdal.GA_ReadOnly)
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
    
    # Set cells with nodata_value to NaN
    if nodata_value:
        raster_values[raster_values == nodata_value] = np.nan

    # Calculate the mean and standard deviation for each cell across all bands
    mean = np.nanmean(raster_values, axis = 0)
    std_dev = np.nanstd(raster_values, axis = 0)

    # Calculate the z-scores for each cell in the raster array
    zscore_array = (raster_values - mean) / std_dev
    
    # Create a new single-band raster dataset for the z-scores
    driver = gdal.GetDriverByName('GTiff')
    zscore_dataset = driver.Create(output_path, width, height, 1, gdal.GDT_Float32)

    # Write the z-scores to the single-band raster dataset
    raster_band_index = raster_band - 1
    zscore_band = zscore_dataset.GetRasterBand(1)
    zscore_band.WriteArray(zscore_array[raster_band_index, :, :], 0, 0)

    # Set the geotransform and projection of the output raster same as the input
    zscore_dataset.SetGeoTransform(dataset.GetGeoTransform())
    zscore_dataset.SetProjection(dataset.GetProjection())
    
    # Close the datasets
    dataset = None
    zscore_dataset = None
    
    print("Z-scores calculated and saved to:", output_path)