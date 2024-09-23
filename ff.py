from os import PathLike, makedirs
import numpy as np
import datetime
from chunk_development.ras_mesh_data import fetch_mesh_data
from chunk_development.spatial_utils import merge_shapes, copy_shape
from chunk_development.WSEChunks import WSEChunks
from chunk_development.spatial_utils import *
from dataset_development.flow_derivative_per_cell import flow_derivative
from dataset_development.zscore_per_cell import zscore_per_cell
from dataset_development.cov_per_cell import cov_per_cell


def main() -> None:

    datasets.get_percentile_DV2_per_chunk_rasters(
        ras_plan_hdf_file = r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_inputs\test_hdf_file\lbr1.p04.hdf",
        scoped_stream_network = r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_inputs\streams_layer\streams_scope.shp",
        DV2_raster = r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_inputs\RAS\100yr\D _ V^2 (Max).Terrain.hydroDEM.tif",
        percentiles = [10,20,30,40,50,60,70,80,90],
        out_directory = r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_outputs\datasets\percentile_DV2_per_chunk"
    )

    return


class chunks(object):

    def __init__(self) -> None:
        pass

    @staticmethod
    def from_backwater_cells(
        ras_plan_hdf_file: PathLike, 
        scoped_stream_network: PathLike,
        out_directory: PathLike, 
        out_name: str = "backwater_chunks"
    ) -> PathLike:
        mesh_areas = fetch_mesh_data(ras_plan_hdf_file)
        chunk_shapes = list()
        for me in mesh_areas:
            me.trace_backwater_chunks(scoped_stream_network)
            chunk_shapes.append(
                me.backwater_chunks_to_shapefile(
                    scoped_stream_network,
                    r"C:\temp",
                    out_name+f"_{me.name}"
                )
            )
        if len(chunk_shapes) > 1:
            chunk_shape = merge_shapes(
                chunk_shapes,
                scoped_stream_network,
                out_directory,
                out_name
            )
        else:
            chunk_shape = copy_shape(
                chunk_shapes[0],
                out_directory,
                out_name
            )
        return chunk_shape

    @staticmethod
    def from_neighbor_cells(
        ras_plan_hdf_file: PathLike, 
        scoped_stream_network: PathLike,
        out_directory: PathLike, 
        out_name: str = "neighbor_chunks"
    ) -> PathLike:
        mesh_areas = fetch_mesh_data(ras_plan_hdf_file)
        chunk_shapes = list()
        for me in mesh_areas:
            me.trace_neighbor_chunks(scoped_stream_network)
            chunk_shapes.append(
                me.neighbor_chunks_to_shapefile(
                    scoped_stream_network,
                    r"C:\temp",
                    out_name+f"_{me.name}"
                )
            )
        if len(chunk_shapes) > 1:
            chunk_shape = merge_shapes(
                chunk_shapes,
                scoped_stream_network,
                out_directory,
                out_name
            )
        else:
            chunk_shape = copy_shape(
                chunk_shapes[0],
                out_directory,
                out_name
            )
        return chunk_shape

    @staticmethod
    def from_wse_raster(
        WSE_raster: PathLike, 
        scoped_stream_network: PathLike,
        output_folder: PathLike, 
        out_name: str, 
        interval: int = 2
    ) -> PathLike:
        chunk_shape = WSEChunks(
            WSE_raster, 
            scoped_stream_network,
            output_folder, 
            out_name, 
            interval
        )
        return chunk_shape


class datasets(object):

    def __init__(self) -> None:
        pass

    @staticmethod
    def get_percentile_DV2_per_chunk_rasters(
        ras_plan_hdf_file: PathLike,
        scoped_stream_network: PathLike,
        DV2_raster: PathLike,
        percentiles: list[int] = [10,20,30,40,50,60,70,80,90],
        out_directory: PathLike = None,
    ) -> None:
        intermediate_folder = path.join(out_directory, "_intermediate_data")
        if not path.exists(intermediate_folder):
            makedirs(intermediate_folder)

        me = fetch_mesh_data(ras_plan_hdf_file)[0]
        cell_pnts_shp = me.cell_points_to_shapefile(scoped_stream_network, intermediate_folder)
        me.trace_floodplains(scoped_stream_network)
        fp_mask = me.floodplain_cells_to_mask_shape(scoped_stream_network, intermediate_folder)
        me.trace_neighbor_chunks(scoped_stream_network)
        chunk_polys = me.neighbor_chunks_to_polygons(scoped_stream_network, True, intermediate_folder)

        start = datetime.datetime.now(); print(f"started getting raster percentiles for chunks at {start}...")
        rast = arcpy.Raster(DV2_raster)
        chunk_percentile_vals = dict() # {chunk_id : [percentile values]}
        # cnt = 0
        # total = len(chunk_polys.keys())
        for chunk_id, chunk_poly in chunk_polys.items():
            # cnt+=1; print(f"    getting raster percentiles for chunk {cnt} of {total}")
            chunk_rast = ExtractByMask(rast, chunk_poly)
            chunk_percentile_vals[chunk_id] = np.nanpercentile(chunk_rast.read(), percentiles)
        print(f"completed raster percentiles for chunks in {datetime.datetime.now() - start}.")

        start = datetime.datetime.now(); print(f"started getting cell percentiles at {start}...")
        cell_percentile_vals = dict() # {cell_id : [percentile values]}
        # cnt = 0
        # total = len(me.cell_residence_neighbor_chunks.keys())
        for cell_id, chunk_ids in me.cell_residence_neighbor_chunks.items():
            # cnt+=1; print(f"    getting cell percentiles for cell {cnt} of {total}")
            cell_percentile_vals[cell_id] = np.mean(np.stack([chunk_percentile_vals[chunk_id] for chunk_id in chunk_ids], axis=1), axis=1)
        print(f"completed cell percentiles in {datetime.datetime.now() - start}.")

        start = datetime.datetime.now(); print(f"started creating tin points at {start}...")
        p_fields = list(f"per_{int(100-p)}" for p in percentiles)
        [arcpy.AddField_management(cell_pnts_shp, pf, "FLOAT") for pf in p_fields]
        query = f"cell_id IN ({str(list(me.cell_residence_neighbor_chunks.keys())).strip('[]')})"
        tin_pnts = arcpy.SelectLayerByAttribute_management(cell_pnts_shp, where_clause=query)
        with arcpy.da.UpdateCursor(tin_pnts, ["cell_id"]+p_fields) as uc:
            for row in uc:
                row[1:] = cell_percentile_vals[row[0]]
                uc.updateRow(row)
        arcpy.CopyFeatures_management(tin_pnts, path.join(intermediate_folder, "tin_pnts"))
        print(f"completed tin points in {datetime.datetime.now() - start}.")

        start = datetime.datetime.now(); print(f"started creating output rasters at {start}...")
        cell_size = arcpy.Describe(rast).meanCellWidth
        arcpy.env.snapRaster = rast
        # cnt = 0
        # total = len(p_fields)
        for field_name in p_fields:
            # cnt+=1; print(f"    creating output raster {cnt} of {total}")
            # create stats raster
            stats_rast_path = path.join(intermediate_folder, f"stats_rast_{field_name}.tif")
            tin = arcpy.CreateTin_3d("stats_tin", arcpy.Describe(scoped_stream_network).spatialReference, [[tin_pnts, field_name, "Mass_Points", ""], [fp_mask, "", "Hard_Clip", ""]])
            stats_rast = arcpy.TinRaster_3d(tin, stats_rast_path, sample_distance= f"CELLSIZE {cell_size}")
            # create out raster
            out_rast_path = path.join(out_directory, f"{field_name}.tif")
            out_rast = Con(Raster(DV2_raster)>Raster(stats_rast), 1)
            out_rast.save(out_rast_path)
        print(f"completed output rasters in {datetime.datetime.now() - start}.")

    @staticmethod
    def get_flow_derivative(
        base_d_raster: PathLike,
        base_v_raster: PathLike,
        base_dv_raster: PathLike, 
        plus_d_raster: PathLike, 
        plus_v_raster: PathLike, 
        plus_dv_raster: PathLike,
        out_directory: PathLike
    ) -> PathLike:
        out_raster = flow_derivative(
            base_d_raster,
            base_v_raster,
            base_dv_raster, 
            plus_d_raster,
            plus_v_raster,
            plus_dv_raster,
            out_directory
        )
        return out_raster

    @staticmethod
    def get_z_score_per_cell(
        multiband_raster: PathLike,
        target_raster_band: int,
        output_raster: PathLike,
        consider_nodata: bool = False
    ) -> PathLike:
        out_raster = zscore_per_cell(
            multiband_raster,
            target_raster_band,
            output_raster,
            consider_nodata
        )
        return out_raster

    @staticmethod
    def get_cov_per_cell(
        multiband_raster: PathLike,
        output_raster: PathLike,
        consider_nodata: bool = False
    ) -> PathLike:
        out_raster = cov_per_cell(
            multiband_raster,
            output_raster,
            consider_nodata
        )
        return out_raster



if __name__ == "__main__":
    main()