from os import PathLike
import pandas as pd
from rashdf import RasGeomHdf
import numpy as np
import geopandas as gpd
from typing import Union

def get_conveyance_and_mannings_curves(
    geom_hdf: PathLike,
    to_dataframe: bool = False
) -> Union[pd.DataFrame, dict]:
    with RasGeomHdf(geom_hdf) as ghdf:
        mesh_names = ghdf.mesh_area_names()
        face_curves = {
            "mesh_name": [], 
            "face_id": [], 
            "elevation": [], 
            "mannings_n": [], 
            "area": [],
            "wetted_perimeter": [],
            "hydraulic_radius": [],
            "conveyance": []
        } if to_dataframe else {}
        for name in mesh_names:
            FACES_AREA_ELEV_INFO = ghdf[rf"/Geometry/2D Flow Areas/{name}/Faces Area Elevation Info"][()]
            FACES_AREA_ELEV_VALUES = ghdf[rf"/Geometry/2D Flow Areas/{name}/Faces Area Elevation Values"][()]
            if to_dataframe:
                face_curves["mesh_name"] += [name] * FACES_AREA_ELEV_VALUES.shape[0]
                for face_id, (start, cnt) in enumerate(FACES_AREA_ELEV_INFO):
                    face_curves["face_id"] += [face_id] * cnt
                    elev, area, wet_perim, mann_n = FACES_AREA_ELEV_VALUES[start:start+cnt].T
                    face_curves["elevation"] += list(elev)
                    face_curves["mannings_n"] += list(mann_n)
                    face_curves["area"] += list(area)
                    face_curves["wetted_perimeter"] += list(wet_perim)
                    face_curves["hydraulic_radius"] += list(np.nan_to_num(area / wet_perim))
                    face_curves["conveyance"] += list(np.nan_to_num(1.486 / mann_n * (area * ((area / wet_perim) ** (2 / 3)))))
            else:
                face_curves[name] = {}
                for face_id, (start, cnt) in enumerate(FACES_AREA_ELEV_INFO):
                    elev, area, wet_perim, mann_n = FACES_AREA_ELEV_VALUES[start:start+cnt].T
                    face_curves[name][face_id] = dict(
                        elevation = elev,
                        mannings_n = mann_n,
                        area = area,
                        wetted_perimeter = wet_perim,
                        hydraulic_radius = np.nan_to_num(area / wet_perim),
                        conveyance = np.nan_to_num(1.486 / mann_n * (area * ((area / wet_perim) ** (2 / 3))))
                    )
        return pd.DataFrame(face_curves) if to_dataframe else face_curves
    
def evaluate_conveyance_reduction(
    initial_geom_hdf: PathLike,
    updated_geom_hdf: PathLike
) -> gpd.GeoDataFrame:
    with RasGeomHdf(initial_geom_hdf) as ghdf1, RasGeomHdf(updated_geom_hdf) as ghdf2:
        mesh_faces = ghdf1.mesh_cell_faces()
        assert mesh_faces.to_json() == ghdf2.mesh_cell_faces().to_json()
    df1 = get_conveyance_and_mannings_curves(initial_geom_hdf)
    df2 = get_conveyance_and_mannings_curves(updated_geom_hdf)
    mesh_faces["delta_q_perc"] = mesh_faces.apply(
        lambda row: round(
            (
                np.median(df2[row["mesh_name"]][row["face_id"]]["conveyance"]) /
                np.median(df1[row["mesh_name"]][row["face_id"]]["conveyance"])
            ) - 1,
            2
        ) * 100,
        axis = 1
    )
    mesh_faces["delta_n_perc"] = mesh_faces.apply(
        lambda row: round(
            (
                df2[row["mesh_name"]][row["face_id"]]["mannings_n"][0] /
                df1[row["mesh_name"]][row["face_id"]]["mannings_n"][0]
            ) - 1,
            2
        ) * 100,
        axis = 1
    )
    return mesh_faces

if __name__ == "__main__":
    gdf = evaluate_conveyance_reduction(
        r"C:\Users\USJB713989\Downloads\OneDrive_2024-09-25\Briar Creek\Base Geometry\Briar_Creek_WS.g02.hdf",
        r"C:\Users\USJB713989\Downloads\OneDrive_2024-09-25\Briar Creek\H1 to H5 Nval 10%  Reduction\Briar_Creek_WS.g03.hdf"
    )
    print(gdf)
    gdf.to_file(r"C:\Users\USJB713989\Downloads\OneDrive_2024-09-25\Briar Creek\conv.shp")