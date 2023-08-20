from os import path, PathLike
import numpy as np
from scipy.spatial import KDTree
from chunk_development.cell import cell
from chunk_development.face import face
import arcpy
from arcpy.sa import *
from arcpy.ia import *
arcpy.env.overwriteOutput = True
arcpy.env.addOutputsToMap = 0
arcpy.env.workspace = arcpy.env.scratchGDB
scratch_folder = arcpy.env.scratchFolder


class BadProjectionError(Exception):
    pass


def faces_to_polylines(face_dict: dict, projection_file: PathLike) -> dict:
    """
    face_dict: {face_id : face object}
    """
    spatial_reference=arcpy.Describe(projection_file).spatialReference
    face_polylines = dict()
    # cnt = 0
    # total = len(face_dict.keys())
    for face_id, face_obj in face_dict.items():
        # cnt+=1; print(f"    creating face polyline {cnt} of {total}")
        face_polylines[face_id] = arcpy.Polyline(arcpy.Array([arcpy.Point(*pnt) for pnt in face_obj.coordinates]), spatial_reference)
    return face_polylines


def cells_to_polygons(cell_dict: dict, projection_file: PathLike) -> dict:
    spatial_reference=arcpy.Describe(projection_file).spatialReference
    cell_polygons = dict()
    # cnt = 0
    # total = len(cell_dict.keys())
    for cell_id, cell_obj in cell_dict.items():
        # cnt+=1; print(f"    creating cell polygon {cnt} of {total}")
        cell_polygons[cell_id] = arcpy.Polygon(arcpy.Array([arcpy.Point(*cp) for cp in cell_obj.cleaned_coordinates]), spatial_reference)
    return cell_polygons


def cells_to_mask_shape(
    cell_ids: set,
    cell_polygons: dict,
    out_directory: PathLike = None,
    out_name: str = "mask"
) -> arcpy.Polygon:
    """
    cell_ids: {cell_ids}
    cell_polygons: {cell_id : cell polygon object}
    """
    if not out_directory: out_directory = arcpy.env.workspace
    mask = arcpy.Dissolve_management(list(cell_polygons[cell_id] for cell_id in cell_ids), path.join(out_directory, out_name))
    return arcpy.Describe(mask).catalogPath


def cell_points_to_shape(
    cells: dict[int:cell], 
    projection_file: PathLike, 
    out_directory: PathLike = None, 
    out_name: str = "cells_pnts"
) -> PathLike:
    """
    cells: {cell_id : cell object}
    """
    if not out_directory: out_directory = arcpy.env.workspace
    cell_pnts_shp = arcpy.CreateFeatureclass_management(
        out_directory, 
        out_name, 
        "point", 
        spatial_reference=arcpy.Describe(projection_file).spatialReference
    )
    arcpy.AddField_management(cell_pnts_shp, "cell_id", "LONG")
    with arcpy.da.InsertCursor(cell_pnts_shp, ["SHAPE@", "cell_id"]) as cur:
        for cell_id, cell_obj in cells.items():
            cur.insertRow([arcpy.PointGeometry(arcpy.Point(*cell_obj.center_coordinates)), cell_id])
    return arcpy.Describe(cell_pnts_shp).catalogPath


def chunks_to_polygons(
    chunk_dict: dict[int:set[int]],
    cell_polygons: dict[int:arcpy.Polygon],
    write_shapefile: bool = False,
    out_directory: PathLike = None
) -> dict:
    """
    chunk_dict: {cell_id : {cell_ids}}
    cell_polygons: {cell_id : cell polygon object}
    """
    chunk_polygons = dict()
    # cnt = 0
    # total = len(chunk_dict.keys())
    chunk_features = list()
    for chunk_id, cell_ids in chunk_dict.items():
        # cnt+=1; print(f"    creating chunk polygon {cnt} of {total}")
        chunk_poly = arcpy.Dissolve_management(list(cell_polygons[cell_id] for cell_id in cell_ids), f"memory/chunk_poly_{chunk_id}")
        chunk_features.append(arcpy.Describe(chunk_poly).catalogPath)
        chunk_poly = next(arcpy.da.SearchCursor(chunk_poly, ["shape@"]))[0]
        chunk_polygons[chunk_id] = chunk_poly
    if write_shapefile and out_directory:
        arcpy.Merge_management(chunk_features, path.join(out_directory, "chunks"))
    return chunk_polygons # {chunk id : chunk polygon object}


def add_burn_id(layer: PathLike, field_name: str = "burn_id") -> None:
    try: arcpy.DeleteField_management(layer, field_name)
    except: pass
    arcpy.AddField_management(layer, field_name, "LONG")
    arcpy.CalculateField_management(layer, field_name, f"!{arcpy.Describe(layer).oidFieldName}!+1")


def copy_shape(
    source: PathLike, 
    out_directory: PathLike, 
    out_name: str = "out_shape"
) -> PathLike:
    dest = path.join(out_directory, out_name)
    arcpy.CopyFeatures_management(source, dest)
    return dest


def check_projection(scoped_stream_network: PathLike) -> None:
    sr = arcpy.Describe(scoped_stream_network).spatialReference
    if not sr.PCSCode or sr.PCSCode == 0:
        raise BadProjectionError("input stream network is not georeferenced in a state plane coordinate system. please reproject to state plane and try again.")


def streams_to_points_array(
    stream_network_layer: PathLike,
    mesh_perimeter_shape: PathLike,
    point_interval_ft: int = 100
) -> np.ndarray:
    start_pnts = arcpy.GeneratePointsAlongLines_management(
        stream_network_layer, 
        "start_pnts", 
        Distance=f"{point_interval_ft} Feet", 
        Include_End_Points="END_POINTS"
    )
    start_pnts_clipped = arcpy.SelectLayerByLocation_management(start_pnts, "INTERSECT", mesh_perimeter_shape)
    return np.array(tuple(row[0] for row in arcpy.da.SearchCursor(start_pnts_clipped, "Shape@XY")))


def perimeter_to_shapefile(
    mesh_perimeter_coordinates: np.ndarray, # np.ndarray[np.ndarray[float]]
    projection_file: PathLike, 
    out_directory: PathLike = arcpy.env.workspace, 
    out_name: str = "mesh_perimeter"
) -> PathLike:
    mesh_perimeter_shp = arcpy.CreateFeatureclass_management(
        out_directory, 
        out_name, 
        "POLYGON", 
        spatial_reference=arcpy.Describe(projection_file).spatialReference
    )
    with arcpy.da.InsertCursor(mesh_perimeter_shp, ["SHAPE@"]) as cur:
        cur.insertRow([arcpy.Polygon(arcpy.Array([arcpy.Point(*point) for point in mesh_perimeter_coordinates]))])
    return arcpy.Describe(mesh_perimeter_shp).catalogPath


def merge_shapes(
    shapes: list[PathLike], 
    projection_file: PathLike, 
    out_directory: PathLike, 
    out_name: str = "chunks"
) -> PathLike:
    out = path.join(out_directory, out_name)
    arcpy.Merge_management(shapes, out)
    return out


def point_spatial_index(input_points: np.ndarray = None, xy_array: np.ndarray = None):
    """
    input_points - this is an array of point xy coord arrays array(array(x,y), array(x,y), ...) that need to be associated with their nearest neighbor
    xy_array - an array of point xy coord arrays array(array(x,y), array(x,y), ...) of potential nearest neighbors to the input_points
    """
    mytree = KDTree(xy_array)
    nearest_neighbor_dists, nearest_neighbor_indexes = mytree.query(input_points)
    return nearest_neighbor_indexes


def scratch_me() -> None:
    scratch_db = arcpy.env.scratchGDB
    dir, filename = path.split(scratch_db)
    try:
        if arcpy.Exists(scratch_db):
            try: arcpy.Delete_management(scratch_db)
            except Exception as e: arcpy.AddMessage(f"could not rewrite scratch database: {e}")
        arcpy.CreateFileGDB_management(dir, filename)
    except: pass