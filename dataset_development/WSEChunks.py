"""
Script to create polygon chunks based on WSE intervals
Inputs: 
  Single-band WSE raster (Horizontal render)
  Interval (1ft, 2ft, etc)
  Scoped streamlines: Used to clip out polygons

Output:
  Polygon shapefile containing interval chunks. Holes filled in (except at polygon boundaries)
"""
import arcpy
import math

def WSEChunks(WSE_raster, interval=2):
    # Environment Settings
    desc = arcpy.Describe(WSE_raster)
    spatialReference = desc.spatialReference
    cellSize = desc.meanCellHeight
    extent = desc.extent
    workspace = desc.path

    arcpy.env.snapRaster = WSE_raster
    arcpy.env.outputCoordinateSystem = spatialReference
    arcpy.env.extent = extent
    arcpy.env.cellSize = cellSize
    arcpy.env.workspace = workspace
    arcpy.env.overwriteOutput = True

    minWSE = math.floor(float(arcpy.GetRasterProperties_management(WSE_raster, "MINIMUM")[0]))
    maxWSE = math.ceil(float(arcpy.GetRasterProperties_management(WSE_raster, "MAXIMUM")[0]))

    # Reclassify. Set range to lower bound
    arcpy.AddMessage("Reclassifying raster")
    reclassify_range = []
    workingWSE = minWSE

    while workingWSE < maxWSE:
        reclassify_range.append([workingWSE, workingWSE+interval, int(workingWSE)])
        workingWSE += interval

    # arcpy.AddMessage(reclassify_range)

    WSE_Reclassify = arcpy.sa.Reclassify(WSE_raster, "Value", arcpy.sa.RemapRange(reclassify_range))
    # WSE_Reclassify.save(workspace + r"\WSE_Reclassify.tif")

    # Raster to Polygon and fill in holes
    arcpy.AddMessage("Raster to Polygon")
    WSE_chunks = arcpy.conversion.RasterToPolygon(WSE_Reclassify, workspace + r"\WSE_chunks.shp")
    WSE_chunks_eliminate = arcpy.EliminatePolygonPart_management(WSE_chunks, workspace + r"\WSE_chunks_filled.shp", "PERCENT", part_area_percent=90)

    # Clip polygons to scoped stream

    # # Delete Intermediate Data
    arcpy.Delete_management(WSE_chunks)

    # Add to map
    arcpy.AddMessage("Adding to Map")
    m = arcpy.mp.ArcGISProject("CURRENT").listMaps('*')[0]
    m.addDataFromPath(arcpy.Describe(WSE_chunks_eliminate).catalogPath)

    return WSE_chunks_eliminate


if __name__ == "__main__":

    WSE_raster = arcpy.GetParameterAsText(0)
    interval = arcpy.GetParameter(1)

    WSEChunks(WSE_raster, interval)