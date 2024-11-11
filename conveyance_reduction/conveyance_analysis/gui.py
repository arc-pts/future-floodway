from tkinter import Tk, Button, Label, Text, ttk, StringVar
from tkinter.filedialog import askopenfilename
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import ctypes
from os import PathLike
from sys import exit
from conveyance_analysis.ras_conveyance_curves import (
    get_conveyance_and_mannings_curves,
    plot_curves,
    get_mesh_names,
    get_face_ids
)

class DesktopPathError(Exception):
    pass

class Gui(Tk):
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)   

        self._ras_curves = None  
        self._mesh_names = None  
        self._face_ids = None   

        # root config
        self.title("Generate Breaklines")
        self.bind_all("<Button-1>", lambda event: event.widget.focus_set())
        # self.geometry("600x280")
        self.attributes('-alpha',0.90)
        self.configure(background='gray30')

        # hdf1 path
        self.hdf1_label = Label(
            self, 
            text="Base Model Plan HDF File:", 
            bg='gray30', 
            fg="cornflower blue", 
            font=12
        )
        self.hdf1_label.pack()
        self.hdf1_path_txt = Text(self, height=2, width=100, padx=5, pady=5, wrap="word")
        # self.hdf1_path_txt.insert("1.0", r"C:\Users\USJB713989\Michael Baker International\PTS3 Innovations - FW HZ SO3 - FW HZ SO3\Data\Flood Hazard Zones\Briar Creek\Base Geometry\Briar_Creek_WS.g02.hdf") # default path
        self.hdf1_path_txt.pack(padx=5, pady=5)
        self.hdf1_browse = Button(self, text="Browse", command=self.browse_hdf1_path)
        self.hdf1_browse.pack(pady=5)

        # hdf2 path
        self.hdf2_label = Label(
            self, 
            text="Updated Model Plan HDF File:", 
            bg='gray30', 
            fg="cornflower blue", 
            font=12
        )
        self.hdf2_label.pack()
        self.hdf2_path_txt = Text(self, height=2, width=100, padx=5, pady=5, wrap="word")
        # self.hdf2_path_txt.insert("1.0", r"C:\Users\USJB713989\Michael Baker International\PTS3 Innovations - FW HZ SO3 - FW HZ SO3\Data\Flood Hazard Zones\Briar Creek\H1 to H5 Nval 10% Increase\Briar_Creek_WS.g03.hdf") # default path
        self.hdf2_path_txt.pack(padx=5, pady=5)
        self.hdf2_browse = Button(self, text="Browse", command=self.browse_hdf2_path)
        self.hdf2_browse.pack(pady=5)

        # execute button
        self.plot_button = Button(
            self, 
            text="Plot Curves", 
            command=self.initiate_plot
        )
        self.plot_button.pack(fill="both", padx=5, pady=5)

        # close button
        self.close_button = Button(self, text="Close", command=exit)
        self.close_button.pack(side= "bottom", fill="both", padx=5, pady=5)

        # main loop
        self.mainloop()

    @property
    def ras_curves(self):
        if self._ras_curves is None:
            self._ras_curves = {
                0: get_conveyance_and_mannings_curves(
                    self.hdf1_path_txt.get("1.0","end").strip().strip('"')
                ),
                1: get_conveyance_and_mannings_curves(
                    self.hdf2_path_txt.get("1.0","end").strip().strip('"')
                ),
            }
        return self._ras_curves

    @property
    def mesh_names(self) -> None:
        if self._mesh_names is None:
            if self._ras_curves is not None:
                self._mesh_names = list(self.ras_curves[0].keys())
            else:
                self._mesh_names = get_mesh_names(
                    self.hdf1_path_txt.get("1.0","end").strip().strip('"')
                )
        return self._mesh_names

    @property
    def face_ids(self) -> None:
        if self._face_ids is None:
            if self._ras_curves is not None:
                self._face_ids = list(self.ras_curves[0][
                    self.mesh_name_str.get().strip().strip('"')
                ].keys())
            else:
                self._face_ids = get_face_ids(
                    self.hdf1_path_txt.get("1.0","end").strip().strip('"'),
                    self.mesh_name_str.get().strip().strip('"')
                )
        return self._face_ids

    def get_hdf_file_path(self) -> PathLike:
        Tk().withdraw() # keep the root window from appearing
        return askopenfilename(
            initialdir = "/",
            title = "Select RAS plan HDF file.",
            filetypes = [("HDF files","*.hdf")]
        )

    def show_error(self, error: str) -> None:
        MessageBox = ctypes.windll.user32.MessageBoxW
        MessageBox(None, error, 'ERROR', 0x40000)
        exit()

    def browse_hdf1_path(self):
        self.hdf1_path_txt.delete("1.0","end")
        self.hdf1_path_txt.insert("1.0", self.get_hdf_file_path())

    def browse_hdf2_path(self):
        self.hdf2_path_txt.delete("1.0","end")
        self.hdf2_path_txt.insert("1.0", self.get_hdf_file_path())

    def mesh_name_changed(self, *args, **kwargs):
        self._face_ids = None
        self.update_plot()

    def initiate_plot(self) -> None:
        try:

            # self.geometry("1200x800")

            # mesh name
            Label(
                self, 
                text="Mesh Name", 
                bg='gray30', 
                fg="cornflower blue", 
                font=12
            ).pack()
            self.mesh_name_str = StringVar()
            self.mesh_name_cb = ttk.Combobox(self, textvariable=self.mesh_name_str)
            self.mesh_name_cb['values'] = self.mesh_names
            self.mesh_name_cb.current(0)
            self.mesh_name_str.trace_add("write", self.mesh_name_changed)
            self.mesh_name_cb.pack()

            # face id
            Label(
                self, 
                text="Cell Face ID", 
                bg='gray30', 
                fg="cornflower blue", 
                font=12
            ).pack()
            self.face_id_txt = Text(self, height=1, width=12, padx=5, pady=5)
            self.face_id_txt.insert("1.0", "0")
            self.face_id_txt.pack()
            self.face_id_txt.bind("<FocusOut>", self.update_plot)

            self.update_plot()

        except Exception as e:
            self.show_error(str(e))

    def update_plot(self, *args, **kwargs) -> None:
        try:
            self.hdf1_label.pack_forget()
            self.hdf1_path_txt.pack_forget()
            self.hdf1_browse.pack_forget()
            self.hdf2_label.pack_forget()
            self.hdf2_path_txt.pack_forget()
            self.hdf2_browse.pack_forget()
            self.plot_button.pack_forget()
            self.cw.pack_forget()
        except:
            pass
        mesh_name=self.mesh_name_str.get().strip().strip('"')
        face_id=int(self.face_id_txt.get("1.0","end").strip().strip('"'))
        try:
            fig = plot_curves(
                self.ras_curves[0][mesh_name][face_id]["conveyance"],
                self.ras_curves[0][mesh_name][face_id]["mannings_n"],
                self.ras_curves[0][mesh_name][face_id]["elevation"],
                self.ras_curves[1][mesh_name][face_id]["conveyance"],
                self.ras_curves[1][mesh_name][face_id]["mannings_n"],
                self.ras_curves[1][mesh_name][face_id]["elevation"]
            )
        except KeyError as ke:
            fig = plot_curves(error_string=f"There is no cell face with ID = {ke}.")

        canvas = FigureCanvasTkAgg(fig, self)   
        canvas.draw() 
        self.cw = canvas.get_tk_widget()
        self.cw.pack(padx=5, pady=5)     
