from os import PathLike, path
import arcpy
from arcpy.sa import *
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")


def flow_derivative(
    base_d_raster: PathLike,
    base_v_raster: PathLike,
    base_dv_raster: PathLike, 
    plus_d_raster: PathLike, 
    plus_v_raster: PathLike, 
    plus_dv_raster: PathLike,
    out_directory: PathLike
) -> None:
    base_d_raster = Raster(base_d_raster); plus_d_raster = Raster(plus_d_raster)
    base_v_raster = Raster(base_v_raster); plus_v_raster = Raster(plus_v_raster)
    base_dv_raster = Raster(base_dv_raster); plus_dv_raster = Raster(plus_dv_raster)
    flow_diff = Con(base_dv_raster, (plus_dv_raster-base_dv_raster)/base_dv_raster)
    flow_diff.save(path.join(out_directory, "flow_diff.tif"))
    depth_diff = Con(base_d_raster, (plus_d_raster-base_d_raster)/base_d_raster)
    depth_diff.save(path.join(out_directory, "depth_diff.tif"))
    vel_diff = Con(base_v_raster, (plus_v_raster-base_v_raster)/base_v_raster)
    vel_diff.save(path.join(out_directory, "vel_diff.tif"))
    fd_rast = Con(base_dv_raster, (flow_diff/depth_diff))
    fd_rast.save(path.join(out_directory, "flow_derivative.tif"))