class cell(object):
    def __init__(self, cell_id: int) -> None:
        self.id: int = cell_id
        self.center_coordinates: tuple[float] = None
        self.max_wse: float = None
        self.min_elev: float = None
        self.face_ids: set[int] = set()
        self.opposite_cell_ids: dict[int:int] = dict() # {face_id : opposite cell_id}
        self.connecting_face_ids: dict[int:int] = dict() # {opposite cell_id : face_id}
        self.chunk_fellows: set[int] = set([cell_id]) # {cell_id, cell_id, ...}