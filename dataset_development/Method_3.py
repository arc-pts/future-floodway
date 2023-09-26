from osgeo import gdal
import numpy as np

# Inputs
raster_path = 'data/0_source/composite_AUS.tif'
output_path = 'data/2_output/Method_3_Scores.tif'

# FUNCTION: N/A
#
# DESCRIPTION: Calculates the 'Method 3' score for a  multi-band raster. An 
#   equal weight is applied to each band. Returns a single band raster of 
#   calculated scores. Any value not in the defined ranges is assigned a 
#   value of zero.
#
# INPUTS: 
#   raster_path - the path to the multi-band raster
#   output_path - the single-band, z-score raster
#
# CODE:
#
# Open the multiband raster dataset
dataset = gdal.Open(raster_path, gdal.GA_ReadOnly)
if dataset is None:
    print("Failed to open the raster dataset.")

# Get the number of bands and dimensions of the raster
num_bands = dataset.RasterCount
width = dataset.RasterXSize
height = dataset.RasterYSize
    
# Create an empty array to store the cell values of all bands
raster_values = np.zeros((num_bands, height, width), dtype=np.float32)

for band_index in range(num_bands):
    band = dataset.GetRasterBand(band_index + 1)
    array = band.ReadAsArray()

    # Define conditions (based on score trendline)
    condition1 = array == 1
    condition2 = array == 2
    condition3 = array == 3
    condition4 = array == 4
    condition5 = array == 5
    condition6 = array == 6

    # Define calculations for each condition (based on score trendline)
    calculation1 = 1
    calculation2 = 2
    calculation3 = 4
    calculation4 = 6
    calculation5 = 8
    calculation6 = 10

    # Apply conditional calculations
    result_array = np.where(condition1, calculation1,
                    np.where(condition2, calculation2, 
                     np.where(condition3, calculation3, 
                      np.where(condition4, calculation4,
                       np.where(condition5, calculation5, 
                        np.where(condition6, calculation6, 0))))))

    raster_values[band_index, :, :] = result_array

# Calculate the mean across all bands as the final score
final_score = np.nanmean(raster_values, axis = 0)

# Create a new single-band raster dataset for the output
driver = gdal.GetDriverByName('GTiff')
score_dataset = driver.Create(output_path, width, height, 1, gdal.GDT_Float32)

# Write the scores to the single-band raster dataset
score_band = score_dataset.GetRasterBand(1)
score_band.WriteArray(final_score, 0, 0)

# Set the geotransform and projection of the output raster same as the input
score_dataset.SetGeoTransform(dataset.GetGeoTransform())
score_dataset.SetProjection(dataset.GetProjection())

# Close the datasets
dataset = None
score_dataset = None

print("Scores calculated and saved to:", output_path)