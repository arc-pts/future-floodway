from tkinter import Tk, Button, Label, Text
from tkinter.filedialog import askopenfilename, askdirectory
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,  
    NavigationToolbar2Tk
)
import ctypes
from os import PathLike
from conveyance_reduction.ras_conveyance_curves import (
    get_conveyance_and_mannings_curves,
    plot_curves
)

class DesktopPathError(Exception):
    pass

class Gui(Tk):
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)   

        # ras curve data
        self._ras_curves = None     

        # root config
        self.title("Generate Breaklines")
        self.geometry("600x400")
        self.attributes('-alpha',0.90)
        self.configure(background='gray30')

        # hdf1 path
        Label(
            self, 
            text="Base Model Plan HDF File:", 
            bg='gray30', 
            fg="cornflower blue", 
            font=12
        ).pack()
        self.hdf1_path_txt = Text(self, height=2, width=53, padx=5, pady=5)
        self.hdf1_path_txt.insert("1.0", r"C:\Users\USJB713989\Michael Baker International\PTS3 Innovations - FW HZ SO3 - FW HZ SO3\Data\Flood Hazard Zones\Briar Creek\Base Geometry\Briar_Creek_WS.g02.hdf") # default path
        self.hdf1_path_txt.pack()
        Button(self, text="Browse", command=self.browse_hdf1_path).pack(pady=5)

        # hdf2 path
        Label(
            self, 
            text="Updated Model Plan HDF File:", 
            bg='gray30', 
            fg="cornflower blue", 
            font=12
        ).pack()
        self.hdf2_path_txt = Text(self, height=2, width=53, padx=5, pady=5)
        self.hdf2_path_txt.insert("1.0", r"C:\Users\USJB713989\Michael Baker International\PTS3 Innovations - FW HZ SO3 - FW HZ SO3\Data\Flood Hazard Zones\Briar Creek\H1 to H5 Nval 10% Increase\Briar_Creek_WS.g03.hdf") # default path
        self.hdf2_path_txt.pack()
        Button(self, text="Browse", command=self.browse_hdf2_path).pack(pady=5)

        # close button
        close_button = Button(self, text="Close", command=self.destroy)
        close_button.pack(side="bottom", fill="both", padx=5, pady=5)

        # execute button
        plot_button = Button(
            self, 
            text="Plot Curves", 
            command=self.initiate_plot
        )
        plot_button.pack(side="bottom", fill="both", padx=5, pady=5)

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

    def get_hdf_file_path(self) -> PathLike:
        Tk().withdraw() # keep the root window from appearing
        return askopenfilename(
            initialdir = "/",
            title = "Select RAS plan HDF file.",
            filetypes = [("HDF files","*.hdf")]
        )

    def get_output_dir_path(self) -> PathLike:
        Tk().withdraw() # keep the root window from appearing
        return askdirectory(
            initialdir = "/",
            title = "Select output directory."
        )

    def show_complete(self) -> None:
        MessageBox = ctypes.windll.user32.MessageBoxW
        MessageBox(None, 'Complete!', ' ', 0x40000)

    def show_error(self, error: str) -> None:
        MessageBox = ctypes.windll.user32.MessageBoxW
        MessageBox(None, error, 'ERROR', 0x40000)
        self.destroy()

    def browse_out_dir(self):
        self.out_dir.delete("1.0","end")
        self.out_dir.insert("1.0", self.get_output_dir_path())

    def browse_hdf1_path(self):
        self.hdf1_path_txt.delete("1.0","end")
        self.hdf1_path_txt.insert("1.0", self.get_dem_file_path())

    def browse_hdf2_path(self):
        self.hdf2_path_txt.delete("1.0","end")
        self.hdf2_path_txt.insert("1.0", self.get_dem_file_path())

    def initiate_plot(self) -> None:
        try:

            # mesh name
            Label(
                self, 
                text="Mesh Name", 
                bg='gray30', 
                fg="cornflower blue", 
                font=12
            ).pack()
            self.mesh_name_txt = Text(self, height=2, width=53, padx=5, pady=5)
            self.mesh_name_txt.insert("1.0", list(self.ras_curves[0].keys())[0])
            self.mesh_name_txt.pack()
            self.mesh_name_txt.bind("<<Modified>>", self.update_plot)

            # face id
            Label(
                self, 
                text="Cell Face ID", 
                bg='gray30', 
                fg="cornflower blue", 
                font=12
            ).pack()
            self.face_id_txt = Text(self, height=2, width=53, padx=5, pady=5)
            self.face_id_txt.insert("1.0", "0")
            self.face_id_txt.pack()
            self.face_id_txt.bind("<<Modified>>", self.update_plot)

            self.update_plot()

        except Exception as e:
            self.show_error(str(e))

    def update_plot(self, *args, **kwargs) -> None:
        try:
            self.cw.pack_forget()
            # self.tbw.pack_forget()
        except:
            pass
        mesh_name=self.mesh_name_txt.get("1.0","end").strip().strip('"')
        face_id=int(self.face_id_txt.get("1.0","end").strip().strip('"'))
        fig = plot_curves(
            self.ras_curves[0][mesh_name][face_id]["conveyance"],
            self.ras_curves[0][mesh_name][face_id]["mannings_n"],
            self.ras_curves[0][mesh_name][face_id]["elevation"],
            self.ras_curves[1][mesh_name][face_id]["conveyance"],
            self.ras_curves[1][mesh_name][face_id]["mannings_n"],
            self.ras_curves[1][mesh_name][face_id]["elevation"]
        )

        canvas = FigureCanvasTkAgg(fig, self)   
        canvas.draw() 
        self.cw = canvas.get_tk_widget()
        self.cw.pack()     
        # toolbar = NavigationToolbar2Tk(canvas, self) 
        # toolbar.update()    
        # self.tbw = canvas.get_tk_widget()
        # self.tbw.pack() 

        self.mesh_name_txt.edit_modified(False)
        self.face_id_txt.edit_modified(False)