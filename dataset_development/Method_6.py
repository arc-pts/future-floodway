from osgeo import gdal
import numpy as np

AC_grid_path = 'data/0_source/ESP_AC.tif'
DV_grid_path = 'data/0_source/ESP_DV.tif'
#DV_output_path = 'data/1_interim/Method_1_DV.tif'
#AC_output_path = 'data/1_interim/Method_1_AC.tif'
output_path = 'data/2_output/Method_6_Scores_ESP.tif'

# FUNCTION: N/A
#
# DESCRIPTION: Calculates the 'Method 1' score for a  multi-band DxV raster
#   and a single band Annual Chance raster. Prior to running this script, the 
#   rasters were clipped to identical boundaries and cell-size, ensuring that 
#   they align into a numpy array. A single band raster of calculated scores
#   is returned. Any value not in the defined ranges is assigned a value of 
#   zero. Equal weight was given to the two rasters. The option exists to 
#   export the interim rasters, but is commented out for now.
#
# INPUTS: 
#   AC_grid_path - the path to the single band Annual Chance raster
#   DV_grid_path - the path to the multi-band DxV raster
#   output_path - the single-band score raster to be created
#
# CODE:
#
# ---------------------
# Flood Frequency Score
# ---------------------
# Processes the annual chance grid

 # Open the percent annual chance raster dataset
AC_input_dataset = gdal.Open(AC_grid_path, gdal.GA_ReadOnly)
if AC_input_dataset is None:
    print("Failed to open the raster dataset.")

# Get the raster dimensions
AC_width = AC_input_dataset.RasterXSize
AC_height = AC_input_dataset.RasterYSize

# Read data
AC_band = AC_input_dataset.GetRasterBand(1)
AC_array = AC_band.ReadAsArray()

# Define condition
AC_condition1 = (AC_array >= 0.2) & (AC_array <= 48)
AC_condition2 = AC_array > 48

# Define calculation for condition
AC_calculation1 = 1.597626 * np.log(AC_array/100) + 11.172
AC_calculation2 = 10

# Apply conditional calculation
AC_scores = np.where(AC_condition1, AC_calculation1, 
             np.where(AC_condition2, AC_calculation2, 0))


# ------------
# Hazard Score
# ------------
# Processes the DV grid

# Open the multiband DV grid dataset
DV_input_dataset = gdal.Open(DV_grid_path, gdal.GA_ReadOnly)
if DV_input_dataset is None:
    print("Failed to open the raster dataset.")

# Get the number of bands and dimensions of the raster
DV_num_bands = DV_input_dataset.RasterCount
DV_width = DV_input_dataset.RasterXSize
DV_height = DV_input_dataset.RasterYSize
    
# Create an empty array to store the cell values of all bands
DV_raster_values = np.zeros((DV_num_bands, DV_height, DV_width), dtype=np.float32)

for band_index in range(DV_num_bands):
    DV_band = DV_input_dataset.GetRasterBand(band_index + 1)
    DV_array = DV_band.ReadAsArray()

    # Define conditions (based on score trendline)
    DV_condition1 = (DV_array > 0) & (DV_array <= 3.2)
    DV_condition2 = (DV_array > 3.2) & (DV_array <= 43.04)
    DV_condition3 = DV_array > 43.04

    # Define calculations for each condition (based on score trendline)
    DV_calculation1 = DV_array * 0.7813
    DV_calculation2 = 2.89249 * np.log(DV_array) - 0.8825
    DV_calculation3 = 10

    # Apply conditional calculations
    DV_result_array = np.where(DV_condition1, DV_calculation1,
                    np.where(DV_condition2, DV_calculation2, 
                     np.where(DV_condition3, DV_calculation3, 0)))

    DV_raster_values[band_index, :, :] = DV_result_array

# Calculate the mean across all bands as the final score
DV_scores = np.nanmean(DV_raster_values, axis = 0)


# --------------
# Combine Scores
# --------------
# Check that the DxV and Annual Chance grids align
if AC_width == DV_width and AC_height == DV_height:

    # Create an empty array to store the individual classification results
    raster_values = np.zeros((2, AC_height, AC_width), dtype=np.float32)
    raster_values[0, :, :] = AC_scores
    raster_values[1, :, :] = DV_scores

    condition1 = raster_values[0, :, :] == 10

    # Apply conditional calculations
    final_score = np.where(condition1, 10, np.nanmean(raster_values, axis = 0))

    # Create a new single-band raster dataset for the output
    driver = gdal.GetDriverByName('GTiff')
    score_dataset = driver.Create(output_path, DV_width, DV_height, 1, gdal.GDT_Float32)

    # Write the scores to the single-band raster dataset
    score_band = score_dataset.GetRasterBand(1)
    score_band.WriteArray(final_score, 0, 0)

    # Set the geotransform and projection of the output raster same as the input
    score_dataset.SetGeoTransform(DV_input_dataset.GetGeoTransform())
    score_dataset.SetProjection(DV_input_dataset.GetProjection())

    # Close the datasets
    dataset = None
    score_dataset = None
        
    print("Scores calculated and saved to:", output_path)
else:
    print("Raster shapes to not match. Unable to calculate score.")