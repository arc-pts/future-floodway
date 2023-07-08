class face(object):
    def __init__(self, face_id: int = None) -> None:
        self.id: int = face_id
        self.coordinates: list[tuple[float]] = list()
        self.min_elev: float = None
        self.opposite_cell_ids: dict[int:int] = dict() # {cell_id : opposite cell_id}