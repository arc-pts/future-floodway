from os import PathLike
import pandas as pd
from rashdf import RasGeomHdf
import numpy as np
import geopandas as gpd
from typing import Union
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

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
        # assert mesh_faces.to_json() == ghdf2.mesh_cell_faces().to_json()
    d1 = get_conveyance_and_mannings_curves(initial_geom_hdf)
    d2 = get_conveyance_and_mannings_curves(updated_geom_hdf)
    mesh_faces["med_elev"] = mesh_faces.apply(
        lambda row: round(
            np.mean(
                [
                    np.median(d2[row["mesh_name"]][row["face_id"]]["elevation"]),
                    np.median(d1[row["mesh_name"]][row["face_id"]]["elevation"])
                ]
            ),
            2
        ),
        axis = 1
    )
    mesh_faces["delta_q_perc"] = mesh_faces.apply(
        lambda row: round(
            (
                np.interp(
                    x=row["med_elev"],
                    xp=d2[row["mesh_name"]][row["face_id"]]["elevation"],
                    fp=d2[row["mesh_name"]][row["face_id"]]["conveyance"]
                ) /
                np.interp(
                    x=row["med_elev"],
                    xp=d1[row["mesh_name"]][row["face_id"]]["elevation"],
                    fp=d1[row["mesh_name"]][row["face_id"]]["conveyance"]
                )
            ) - 1,
            2
        ) * 100,
        axis = 1
    )
    mesh_faces["delta_n_perc"] = mesh_faces.apply(
        lambda row: round(
            (
                d2[row["mesh_name"]][row["face_id"]]["mannings_n"][0] /
                d1[row["mesh_name"]][row["face_id"]]["mannings_n"][0]
            ) - 1,
            2
        ) * 100,
        axis = 1
    )
    return mesh_faces

def plot_curves(
    base_q_vals: np.ndarray,
    base_n_vals: np.ndarray,
    base_z_vals: np.ndarray,
    updated_q_vals: np.ndarray,
    updated_n_vals: np.ndarray,
    updated_z_vals: np.ndarray
) -> Figure:
    figure, axis = plt.subplots(1, 2)
    figure.set_size_inches(14, 4)

    axis[0].set_title("conveyance")
    axis[0].set_xlabel("conveyance")
    axis[0].set_ylabel("elevation")
    axis[0].ticklabel_format(useOffset=False, style="plain")
    axis[0].plot(
        base_q_vals, 
        base_z_vals,
        label="base"
    )
    axis[0].plot(
        updated_q_vals, 
        updated_z_vals,
        label="updated"
    )
    axis[0].legend()

    axis[1].set_title("mannings n")
    axis[1].set_xlabel("mannings n")
    axis[1].set_ylabel("elevation")
    axis[1].ticklabel_format(useOffset=False, style="plain")
    axis[1].plot(
        base_n_vals, 
        base_z_vals,
        label="base"
    )
    axis[1].plot(
        updated_n_vals, 
        updated_z_vals,
        label="updated"
    )
    axis[1].legend()
    return figure

def get_mesh_names(geom_hdf: PathLike) -> list[str]:
    with RasGeomHdf(geom_hdf) as ghdf:
        return ghdf.mesh_area_names()

def get_face_ids(geom_hdf: PathLike, mesh_name: str) -> list[int]:
    with RasGeomHdf(geom_hdf) as ghdf:
        return list(range(ghdf[rf"/Geometry/2D Flow Areas/{mesh_name}/Faces Low Elevation Centroid"].size))
