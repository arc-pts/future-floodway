# conda env setup
Create a conda environment with the required dependencies:
```
conda create -n geo python=3.9
conda activate geo
conda install geopandas
conda install geospatial -c conda-forge
```

# build
```
(geo) > pyinstaller cli.py --name conveyance_analysis --onefile --noconsole
```
