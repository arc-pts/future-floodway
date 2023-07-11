import h5py
import numpy as np
import os
from datetime import datetime
from cell import cell
from face import face
from spatial_utils import *
import sys
MAX_RECURSION_DEPTH = 1000
sys.setrecursionlimit(MAX_RECURSION_DEPTH)


class MissingStreamNetworkError(Exception):
    pass


class mesh_area(object):
    def __init__(self, name: str):
        self.name = name
        self.coordinates: np.ndarray[np.ndarray[float]] = None
        self.faces: dict[int:face] = dict() # {face_id : face object}
        self.cells: dict[int:cell] = dict() # {cell_id : cell object}
        self.trace_sources: list[int] = list() # [cell_ids]
        self.floodplain_cell_ids: set[int] = set() # {cell_ids}
        self.neighbor_chunks: dict[int:set[int]] = dict() # {cell_id : {cell_ids}}
        self.backwater_chunks: dict[int:set[int]] = dict() # {cell_id : {cell_ids}}
        self.recursion_depth: int = 0


    def get_trace_sources(self, stream_network_layer: os.PathLike) -> None:
        """
        self.trace_sources - a list of cell_ids to be used in tracing floodplain cells
        """
        start = datetime.now(); print(f"getting trace sources at {start}...")
        if not stream_network_layer:
            raise MissingStreamNetworkError("Missing a stream network layer.")
        
        # check proj
        check_projection(stream_network_layer)

        # create mesh perimeter shape
        mesh_perimeter_shp = perimeter_to_shapefile(self.coordinates, stream_network_layer)

        # get source points for tracing
        input_points = streams_to_points_array(
            stream_network_layer,
            mesh_perimeter_shp
        )

        xy_array_coords = np.array(tuple(cell_obj.center_coordinates for cell_obj in self.cells.values()), dtype=object)
        nn_indexes = point_spatial_index(input_points, xy_array_coords)
        self.trace_sources = list(list(self.cells.keys())[i] for i in nn_indexes if not np.isnan(self.cells[i].min_elev))
        print(f"completed trace sources in {datetime.now() - start}.")


    def trace_floodplains(self, stream_network_layer: os.PathLike = None) -> None:
        start = datetime.now(); print(f"started tracing floodplains at {start}...")
        if not self.trace_sources:
            self.get_trace_sources(stream_network_layer)
        global visited; visited = set()
        global trace_benchmark; trace_benchmark = dict()
        for cell_id in self.trace_sources:
            try:
                self._trace_floodplain_cells(cell_id); 
            except RecursionError:
                print(f"Max recursion depth reached for Cell ID {cell_id} which prevented flow tracing from that source.")
        print(f"completed tracing floodplains in {datetime.now() - start}.")


    def _trace_floodplain_cells(self, cell_id: cell, initial_face_id: int = None, backwater_elev: float = None) -> bool:
        trace_benchmark.setdefault(cell_id, -99.9)
        self.recursion_depth += 1
        if self.recursion_depth > MAX_RECURSION_DEPTH - 100: self.recursion_depth -= 1; return True
        cell_object = self.cells[cell_id]

        if (not np.isnan(cell_object.min_elev)) and (not cell_id in visited) or (not backwater_elev and trace_benchmark[cell_id] < cell_object.max_wse) or (backwater_elev and trace_benchmark[cell_id] < backwater_elev):
            visited.add(cell_id)
            self.floodplain_cell_ids.add(cell_id)
            trace_benchmark[cell_id] = backwater_elev or cell_object.max_wse

            for face_id in cell_object.face_ids:
                stack_limit_reached = False
                face_object = self.faces[face_id]
                if not initial_face_id or face_object.id != initial_face_id:
                    opposite_cell_object = self.cells[cell_object.opposite_cell_ids[face_id]]
                    if (trace_benchmark[cell_id] >= opposite_cell_object.max_wse) and (trace_benchmark[cell_id] > face_object.min_elev) and not backwater_elev: # check for downstream cells here if not plotting backwater
                        stack_limit_reached = self._trace_floodplain_cells(opposite_cell_object.id, face_object.id)
                    elif (trace_benchmark[cell_id] < opposite_cell_object.max_wse) and (trace_benchmark[cell_id] > opposite_cell_object.min_elev) and (trace_benchmark[cell_id] > face_object.min_elev): # handle backwater cells here
                        stack_limit_reached = self._trace_floodplain_cells(opposite_cell_object.id, face_object.id, trace_benchmark[cell_id])

                    if stack_limit_reached and not backwater_elev:
                        self.trace_sources.append(cell_id); self.recursion_depth -= 1; return False
                    elif stack_limit_reached and backwater_elev:
                        self.recursion_depth -= 1; return True

        self.recursion_depth -= 1; return False


    def trace_backwater_chunks(self, stream_network_layer: os.PathLike = None) -> None:
        start = datetime.now(); print(f"started tracing backwater chunks at {start}...")
        if not self.trace_sources:
            self.get_trace_sources(stream_network_layer)
        global trace_benchmark; trace_benchmark = dict()
        global chunks; chunks = dict()
        for cell_id in self.trace_sources:
            try:
                global ceiling; ceiling = -99.9
                global ceiling_cell_id; ceiling_cell_id = None
                global floor; floor = self.cells[cell_id].max_wse
                global root_cell_id; root_cell_id = cell_id
                global visited; visited = set()
                self._seek_backwater_ceiling(cell_id)
                visited = set()
                self._trace_cells_within_range(ceiling_cell_id)
            except RecursionError:
                print(f"Max recursion depth reached for Cell ID {cell_id} which prevented chunk tracing from that source.")
        self.backwater_chunks = chunks.copy()
        print(f"completed tracing backwater chunks in {datetime.now() - start}.")


    def _seek_backwater_ceiling(self, cell_id: cell) -> None:
        cell_object = self.cells[cell_id]
        visited.add(cell_id)
        global ceiling; ceiling = max(cell_object.max_wse, ceiling)
        if cell_object.max_wse == ceiling:
            global ceiling_cell_id; ceiling_cell_id = cell_object.id
        for face_id in cell_object.face_ids:
            face_object = self.faces[face_id]
            opposite_cell_object = self.cells[cell_object.opposite_cell_ids[face_id]]
            if (not np.isnan(opposite_cell_object.min_elev)) and (opposite_cell_object.id not in visited) and (opposite_cell_object.max_wse >= floor) and (opposite_cell_object.min_elev < floor) and (face_object.min_elev < floor):
                self._seek_backwater_ceiling(opposite_cell_object.id)


    def trace_neighbor_chunks(self, stream_network_layer: os.PathLike = None) -> None:
        start = datetime.now(); print(f"started tracing neighbor chunks at {start}...")
        if not self.trace_sources:
            self.get_trace_sources(stream_network_layer)
        global trace_benchmark; trace_benchmark = dict()
        global chunks; chunks = dict()
        for cell_id in self.trace_sources:
            try:
                global ceiling, floor
                self._seek_neighbor_range(cell_id)
                global root_cell_id; root_cell_id = cell_id
                self._trace_cells_within_range(cell_id)
            except ValueError:
                print(f"no neighbors found for Cell ID {cell_id} which prevented chunk tracing from that source.")
            except RecursionError:
                print(f"Max recursion depth reached for Cell ID {cell_id} which prevented chunk tracing from that source.")
        self.neighbor_chunks = chunks.copy()
        print(f"completed tracing neighbor chunks in {datetime.now() - start}.")


    def _seek_neighbor_range(self, cell_id: cell) -> None:
        cell_object = self.cells[cell_id]
        global ceiling, floor; ceiling = floor = cell_object.max_wse
        for face_id in cell_object.face_ids:
            face_object = self.faces[face_id]
            opposite_cell_object = self.cells[cell_object.opposite_cell_ids[face_id]]
            if (not np.isnan(opposite_cell_object.min_elev)) and (cell_object.max_wse > opposite_cell_object.min_elev) and (cell_object.max_wse > face_object.min_elev):
                ceiling = max(ceiling, opposite_cell_object.max_wse)
                floor = min(floor, opposite_cell_object.max_wse)


    def _trace_cells_within_range(self, cell_id: cell, initial_face_id: int = None, backwater_elev: float = None) -> None:
        chunks.setdefault(root_cell_id, set())
        trace_benchmark.setdefault(cell_id, -99.9)
        cell_object = self.cells[cell_id]
        if (not np.isnan(cell_object.min_elev)) and (cell_id not in chunks[root_cell_id]) or (not backwater_elev and trace_benchmark[cell_id] < cell_object.max_wse) or (backwater_elev and trace_benchmark[cell_id] < backwater_elev):
            chunks[root_cell_id].add(cell_id)
            trace_benchmark[cell_id] = backwater_elev or cell_object.max_wse
            for face_id in cell_object.face_ids:
                face_object = self.faces[face_id]
                if not initial_face_id or face_object.id != initial_face_id:
                    opposite_cell_object = self.cells[cell_object.opposite_cell_ids[face_id]]
                    if (ceiling >= opposite_cell_object.max_wse >= floor):
                        if (trace_benchmark[cell_id] >= opposite_cell_object.max_wse) and (trace_benchmark[cell_id] > face_object.min_elev) and not backwater_elev: # check for downstream cells here if not plotting backwater
                            self._trace_cells_within_range(opposite_cell_object.id, face_object.id)
                        elif (trace_benchmark[cell_id] < opposite_cell_object.max_wse) and (trace_benchmark[cell_id] > opposite_cell_object.min_elev) and (trace_benchmark[cell_id] > face_object.min_elev): # handle backwater cells here
                            self._trace_cells_within_range(opposite_cell_object.id, face_object.id, trace_benchmark[cell_id])


    def cells_to_shapefile(self, cell_ids: set[int], projection_file: os.PathLike, out_directory: os.PathLike, out_name: str = "cells") -> os.PathLike:
        start = datetime.now(); print(f"creating cell shapefile at {start}...")
        out_shape = cells_to_shape(
            dict((cell_id, list([self.faces[face_id] for face_id in self.cells[cell_id].face_ids])) for cell_id in cell_ids),
            projection_file,
            out_directory,
            out_name
        )
        print(f"cell shapefile complete in {datetime.now() - start}.")
        return out_shape
    

    def backwater_chunks_to_shapefile(self, projection_file: os.PathLike, out_directory: os.PathLike, out_name: str = "backwater_chunks") -> os.PathLike:
        start = datetime.now(); print(f"creating chunk shapefile at {start}...")
        out_shape = chunks_to_shape(
            self.backwater_chunks,
            self.cells,
            self.faces,
            projection_file,
            out_directory,
            out_name
        )
        print(f"chunk shapefile complete in {datetime.now() - start}.")
        return out_shape


    def neighbor_chunks_to_shapefile(self, projection_file: os.PathLike, out_directory: os.PathLike, out_name: str = "neighbor_chunks") -> os.PathLike:
        start = datetime.now(); print(f"creating chunk shapefile at {start}...")
        out_shape = chunks_to_shape(
            self.neighbor_chunks,
            self.cells,
            self.faces,
            projection_file,
            out_directory,
            out_name
        )
        print(f"chunk shapefile complete in {datetime.now() - start}.")
        return out_shape


def fetch_mesh_data(ras_plan_hdf_file: os.PathLike) -> list[mesh_area]:
    start = datetime.now(); print(f"getting RAS mesh area data at {start}...")
    mesh_areas = list()
    with h5py.File(ras_plan_hdf_file, 'r') as f:
        mesh_area_names = [str(raw_name).strip("b'") for raw_name in f["/Geometry/2D Flow Areas/Attributes"]["Name"][()]]
        mesh_id = -1
        for name in mesh_area_names:
            mesh_id += 1
            facepoints_coordinates = f[f"/Geometry/2D Flow Areas/{name}/FacePoints Coordinate"][()]
            facepoints_index = f[f"/Geometry/2D Flow Areas/{name}/Faces FacePoint Indexes"][()]
            faces_perimeter_info = f[f"/Geometry/2D Flow Areas/{name}/Faces Perimeter Info"][()]
            faces_perimeter_values = f[f"/Geometry/2D Flow Areas/{name}/Faces Perimeter Values"][()]
            face_min_elev = f[f"/Geometry/2D Flow Areas/{name}/Faces Minimum Elevation"][()]
            face_cells_index = f[f"/Geometry/2D Flow Areas/{name}/Faces Cell Indexes"][()]
            cell_pnt_start, cell_pnt_cnt = f[f"/Geometry/2D Flow Areas/Cell Info"][()][mesh_id]
            internal_cell_point_coordinates = f[f"/Geometry/2D Flow Areas/Cell Points"][()][cell_pnt_start:cell_pnt_start+cell_pnt_cnt]
            all_cell_point_coordinates = f[f"/Geometry/2D Flow Areas/{name}/Cells Center Coordinate"][()]
            cell_max_wsel = f[f"/Results/Unsteady/Output/Output Blocks/Base Output/Summary Output/2D Flow Areas/{name}/Maximum Water Surface"][()]
            cell_min_elev = f[f"/Geometry/2D Flow Areas/{name}/Cells Minimum Elevation"][()]
            cell_face_info = f[f"/Geometry/2D Flow Areas/{name}/Cells Face and Orientation Info"][()]
            cell_face_values = f[f"/Geometry/2D Flow Areas/{name}/Cells Face and Orientation Values"][()]

            # get faces
            face_id = -1
            face_dict: dict[int:face] = dict() # [face_id : face object]
            for item in facepoints_index:
                face_id += 1
                face_obj = face(face_id)
                starting_row, count = faces_perimeter_info[face_id]
                face_obj.coordinates.append(tuple(facepoints_coordinates[item[0]]))
                if count > 0:
                    for row in range(starting_row, starting_row + count):
                        face_obj.coordinates.append(tuple(faces_perimeter_values[row]))
                face_obj.coordinates.append(tuple(facepoints_coordinates[item[1]]))
                face_obj.min_elev = face_min_elev[face_id]
                face_obj.opposite_cell_ids = {face_cells_index[face_id][0] : face_cells_index[face_id][1], face_cells_index[face_id][1] : face_cells_index[face_id][0]}
                face_dict[face_id] = face_obj
            
            # get cells
            cell_id = -1
            cell_dict: dict[int:cell] = dict() # [cell_id : cell object]
            for row in all_cell_point_coordinates:
                cell_id += 1
                cell_obj = cell(cell_id)
                cell_obj.center_coordinates = tuple(row)
                cell_obj.max_wse = cell_max_wsel[0][cell_id]
                cell_obj.min_elev = cell_min_elev[cell_id]
                starting_row, count = cell_face_info[cell_id]
                for face_id,_ in cell_face_values[starting_row:starting_row+count]:
                    cell_obj.face_ids.add(face_id)
                    cell_obj.opposite_cell_ids[face_id] = face_dict[face_id].opposite_cell_ids[cell_id]
                    cell_obj.connecting_face_ids[face_dict[face_id].opposite_cell_ids[cell_id]] = face_id
                cell_dict[cell_id] = cell_obj

            # get mesh areas
            mesh_obj = mesh_area(name)
            mesh_obj.coordinates = f[f"/Geometry/2D Flow Areas/{name}/Perimeter"][()]
            mesh_obj.faces = face_dict.copy()
            mesh_obj.cells = cell_dict.copy()
            mesh_areas.append(mesh_obj)

    print(f"mesh area data complete in {datetime.now() - start}.")
    return mesh_areas


if __name__ == "__main__":
    pass