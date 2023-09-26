# A generic helper function to calculate given percentile values per chunk
# Inputs: 
#   Single-band raster
#   Polygon shapefile containing raster chunk polygons

# Output:
#   A group of percentile rasters

import arcpy
from os import PathLike


def ChunkPercentilesToRaster(value_raster: PathLike, chunk_polygon: PathLike, output_folder: PathLike) -> PathLike:
    # Environment Settings
    desc = arcpy.Describe(value_raster)
    spatialReference = desc.spatialReference
    cellSize = desc.meanCellHeight
    extent = desc.extent

    if output_folder:
        workspace = output_folder
    else:
        workspace = desc.path

    arcpy.AddMessage("Output will be saved in {0}".format(workspace))

    arcpy.env.snapRaster = value_raster
    arcpy.env.outputCoordinateSystem = spatialReference
    arcpy.env.extent = extent
    arcpy.env.cellSize = cellSize
    arcpy.env.workspace = workspace
    arcpy.env.overwriteOutput = True

    percentiles = [10, 20, 30, 40, 50, 60, 70, 80, 90]

    for percent in percentiles:
        arcpy.AddMessage("Calculating {0} percentile".format(percent))
        outName = r"\Floodway_" + str(percent) + "Pct.tif"
        outZonalStats = arcpy.sa.ZonalStatistics(chunk_polygon, "FID", value_raster, "PERCENTILE",
                                "DATA", "CURRENT_SLICE", percent)
        # Get area higher than the percentile
        outMinus = arcpy.sa.Minus(value_raster, outZonalStats)
        # This represents areas greater than the percentile. Actually 100 - percentile
        outSetNull = arcpy.sa.SetNull(outMinus, percent, "VALUE < 0")
        outSetNull.save(workspace + outName)

    arcpy.AddMessage("Output saved to {0}".format(output_folder))
    return output_folder


def addToMap(layer):
    arcpy.AddMessage("Adding to Map")
    m = arcpy.mp.ArcGISProject("CURRENT").listMaps('*')[0]
    m.addDataFromPath(arcpy.Describe(layer).catalogPath)


if __name__ == "__main__":

    value_raster = arcpy.GetParameterAsText(0)
    chunk_polygon = arcpy.GetParameterAsText(1)
    output_folder = arcpy.GetParameterAsText(2)

    wse_chunk_polygon = ChunkPercentilesToRaster(value_raster, chunk_polygon, output_folder)
    # addToMap(wse_chunk_polygon)




