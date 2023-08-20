import h5py
import numpy as np
import os
from datetime import datetime
from chunk_development.cell import cell
from chunk_development.face import face
from chunk_development.spatial_utils import *
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
        self.neighbor_chunks: dict[int:set[int]] = dict() # {root cell_id : {cell_ids}}
        self.cell_residence_neighbor_chunks: dict[int:set[int]] = dict() # {cell_id : {chunk_ids}}
        self.backwater_chunks: dict[int:set[int]] = dict() # {root cell_id : {cell_ids}}
        self.face_polylines: dict = None # {face_id : polyline object}
        self.cell_polygons: dict = None # {cell_id : polygon object}
        self.neighbor_chunk_polygons: dict = None # {chunk_id : polygon object}
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


    def trace_neighbor_chunks(self, stream_network_layer: os.PathLike = None) -> None:
        start = datetime.now(); print(f"started tracing neighbor chunks at {start}...")
        if not self.trace_sources:
            self.get_trace_sources(stream_network_layer)
        global trace_benchmark; trace_benchmark = dict()
        global chunks; chunks = dict()
        global cell_residence_chunks; cell_residence_chunks = dict()
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
        self.cell_residence_neighbor_chunks = cell_residence_chunks.copy()
        print(f"completed tracing neighbor chunks in {datetime.now() - start}.")


    def create_face_polylines(self, projection_file: os.PathLike) -> dict:
        start = datetime.now(); print(f"creating face polylines at {start}...")
        self.face_polylines = faces_to_polylines(self.faces, projection_file)
        print(f"face polylines complete in {datetime.now() - start}.")
        return self.face_polylines


    def create_cell_polygons(self, projection_file: os.PathLike) -> dict:
        start = datetime.now(); print(f"creating cell polygons at {start}...")
        self.cell_polygons = cells_to_polygons(self.cells, projection_file)
        print(f"cell polygons complete in {datetime.now() - start}.")
        return self.cell_polygons


    def floodplain_cells_to_mask_shape(self, stream_network_layer: os.PathLike = None, out_directory: os.PathLike = None, out_name: str = "floodplain_mask") -> os.PathLike:
        start = datetime.now(); print(f"creating mask polygon at {start}...")
        if not self.floodplain_cell_ids:
            self.trace_floodplains(stream_network_layer)
        if not self.cell_polygons:
            self.create_cell_polygons(stream_network_layer)
        mask_shape = cells_to_mask_shape(
            self.floodplain_cell_ids,
            self.cell_polygons,
            out_directory,
            out_name
        )
        print(f"mask polygon complete in {datetime.now() - start}.")
        return mask_shape


    def cell_points_to_shapefile(self, projection_file: os.PathLike, out_directory: os.PathLike, out_name: str = "cell_pnts") -> os.PathLike:
        start = datetime.now(); print(f"creating cell point shapefile at {start}...")
        out_shape = cell_points_to_shape(
            self.cells,
            projection_file,
            out_directory,
            out_name
        )
        print(f"cell point shapefile complete in {datetime.now() - start}.")
        return out_shape


    def neighbor_chunks_to_polygons(self, stream_network_layer: os.PathLike = None) -> dict:
        start = datetime.now(); print(f"creating chunk polygons at {start}...")
        if not self.neighbor_chunks:
            self.trace_neighbor_chunks(stream_network_layer)
        if not self.cell_polygons:
            self.create_cell_polygons(stream_network_layer)
        self.neighbor_chunk_polygons = chunks_to_polygons(self.neighbor_chunks, self.cell_polygons)
        print(f"chunk polygons complete in {datetime.now() - start}.")
        return self.neighbor_chunk_polygons


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
            cell_residence_chunks.setdefault(cell_id, set()).add(root_cell_id)
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

            # get cell cleaned coordinates
            face_id = -1
            starting_row = 0
            count = 0
            for item in facepoints_index:
                face_id += 1
                face_coords = []
                face_coords.append(list(facepoints_coordinates[item[0]]))
                starting_row = faces_perimeter_info[face_id][0]
                count = faces_perimeter_info[face_id][1]
                if count > 0:
                    for row in range(starting_row, starting_row + count):
                        face_coords.append(list(faces_perimeter_values[row]))
                face_coords.append(list(facepoints_coordinates[item[1]]) )
                for cell_id in [int(i) for i in face_cells_index[face_id]]:
                    cell_obj = cell_dict.setdefault(cell_id, cell(cell_id))
                    cell_obj.face_coordinates.append(face_coords)

            global cleaned_coords; global raw_coords

            def clean_coords(end_coord: list = None) -> None:
                if not end_coord:
                    target = raw_coords.pop(0)
                else:
                    i = -1
                    for face in raw_coords:
                        i+=1
                        if end_coord == face[0]:
                            target = raw_coords.pop(i)
                        elif end_coord == face[-1]:
                            target = raw_coords.pop(i)
                            target.reverse()
                [cleaned_coords.append(coord) for coord in target]
                if raw_coords: clean_coords(target[-1])

            for cell_id, cell_obj in cell_dict.items():
                cleaned_coords = []
                raw_coords = cell_obj.face_coordinates.copy()
                clean_coords()
                cell_dict[cell_id].cleaned_coordinates = cleaned_coords

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