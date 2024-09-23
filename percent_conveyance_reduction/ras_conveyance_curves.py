from os import PathLike
import pandas as pd
from rashdf import RasGeomHdf
import numpy as np
import geopandas as gpd

def get_conveyance_and_mannings_curves(geom_hdf: PathLike) -> pd.DataFrame:
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
        }
        for name in mesh_names:
            FACES_AREA_ELEV_INFO = ghdf[rf"/Geometry/2D Flow Areas/{name}/Faces Area Elevation Info"][()]
            FACES_AREA_ELEV_VALUES = ghdf[rf"/Geometry/2D Flow Areas/{name}/Faces Area Elevation Values"][()]
            face_curves["mesh_name"] += [name] * FACES_AREA_ELEV_VALUES.shape[0]
            for face_id, (start, cnt) in enumerate(FACES_AREA_ELEV_INFO):
                face_curves["face_id"] += [face_id] * cnt
                elev, area, wet_perim, mann_n = FACES_AREA_ELEV_VALUES[start:start+cnt].T
                face_curves["elevation"] += list(elev)
                face_curves["mannings_n"] += list(mann_n)
                face_curves["area"] += list(area)
                face_curves["wetted_perimeter"] += list(wet_perim)
                face_curves["hydraulic_radius"] += list(area / wet_perim)
                face_curves["conveyance"] += list(1.486 / mann_n * (area * ((area / wet_perim) ** (2 / 3))))
        return pd.DataFrame(face_curves).fillna(0)
    
def evaluate_conveyance_reduction(
    initial_geom_hdf: PathLike,
    updated_geom_hdf: PathLike
) -> gpd.GeoDataFrame:
    with RasGeomHdf(initial_geom_hdf) as ghdf1, RasGeomHdf(updated_geom_hdf) as ghdf2:
        mesh_faces = ghdf1.mesh_cell_faces()
        assert mesh_faces.to_json() == ghdf2.mesh_cell_faces().to_json()
    df1 = get_conveyance_and_mannings_curves(initial_geom_hdf)
    df2 = get_conveyance_and_mannings_curves(updated_geom_hdf)
    assert df1["elevation"].to_list() == df2["elevation"].to_list()
    mesh_faces["delta_conveyance"] = mesh_faces["face_id"].apply(
        lambda face_id: 1 - (
            np.trapz(
                y=df2.loc[df2["face_id"] == face_id]["conveyance"].to_numpy(),
                x=df2.loc[df2["face_id"] == face_id]["elevation"].to_numpy()
            ) / np.trapz(
                y=df1.loc[df1["face_id"] == face_id]["conveyance"].to_numpy(),
                x=df1.loc[df1["face_id"] == face_id]["elevation"].to_numpy()
            )
        )
    )
    mesh_faces["delta_mann_n"] = mesh_faces["face_id"].apply(
        lambda face_id: 1 - (
            np.trapz(
                y=df2.loc[df2["face_id"] == face_id]["mannings_n"].to_numpy(),
                x=df2.loc[df2["face_id"] == face_id]["elevation"].to_numpy()
            ) / np.trapz(
                y=df1.loc[df1["face_id"] == face_id]["mannings_n"].to_numpy(),
                x=df1.loc[df1["face_id"] == face_id]["elevation"].to_numpy()
            )
        )
    )
    return mesh_faces

if __name__ == "__main__":
    pass