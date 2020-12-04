# Night-time_lights_Electrification
Project funded by the World Bank Group and ESMAP to estimate electricity access rate with night-time lights and other satellite-based datasets in developing countries.

The code provided downloads night-time lights data (NTL, VIIRS from NOAA) from Google Earth Engine and correct the inter-annual shift in noisiness using population data.

## Anaconda installation

Download Anaconda Individual Edition for 64 or 32 bits at:
https://www.anaconda.com/products/individual

Install Anaconda and make sure to check the box to add Anaconda3 to the PATH environment variable during the installation process.

## Setting up an operational geospatial Python environment

Launch a command prompt (type 'cmd' in the Windows Start menu and then type enter). Copy and execute:

	conda create --name geo37 python=3.7

This command will create a new Python environment, called geo37, using Python 3.7 (and not the latest version 3.8).
So far, the code is not working properly on Python 3.8 because of some packages dependencies issues.

Activate the new environment by executing:

	conda activate geo37

Then run the following lines to set up the environment:

	conda update --all
	conda install jupyter
	conda config --prepend channels conda-forge
	pip install earthengine-api --upgrade
	pip install geetools

From https://www.lfd.uci.edu/~gohlke/pythonlibs/, download (chose 32 or 64 bits files):

	rasterio‑1.1.8‑cp37‑cp37m‑win_amd64.whl
	GDAL‑3.1.4‑cp37‑cp37m‑win_amd64.whl
	Shapely‑1.7.1‑cp37‑cp37m‑win_amd64.whl
	Fiona‑1.8.18‑cp37‑cp37m‑win_amd64.whl
	geopandas‑0.8.1‑py3‑none‑any.whl
	pyproj‑3.0.0.post1‑cp37‑cp37m‑win_amd64.whl
	Rtree‑0.9.4‑cp37‑cp37m‑win_amd64.whl

Run:

	cd %userprofile%\Downloads
	pip install GDAL‑3.1.4‑cp37‑cp37m‑win_amd64.whl
	pip install rasterio‑1.1.8‑cp37‑cp37m‑win_amd64.whl
	pip install Shapely‑1.7.1‑cp37‑cp37m‑win_amd64.whl
	pip install Fiona‑1.8.18‑cp37‑cp37m‑win_amd64.whl
	pip install geopandas‑0.8.1‑py3‑none‑any.whl
	pip install pyproj‑3.0.0.post1‑cp37‑cp37m‑win_amd64.whl
	pip install Rtree‑0.9.4‑cp37‑cp37m‑win_amd64.whl

	conda install git -y
	pip install git+https://github.com/wpgp/wpgpDownloadPy
	pip install hdx-python-api
	conda install matplotlib

Close the command prompt. The environment is ready to be used.

## Using Night-time_lights_Electrification

Download Night-time_lights_Electrification-main.zip from Github and unzip it a dedicated folder. 
This will be the working directory and downloaded files and outputs will appear there.

Then run Jupyter Notebook by opening a command prompt (cmd) and typing:

	conda activate geo37
	jupyter notebook

This should open a new tab in your default browser.
You should be able to navigate to the Night-time_lights_Electrification-main folder through the Jupyter Notebook interface.
If not, that means you probably put the Night-time_lights_Electrification-main folder somewhere at the root of your hard drive disk
and not in subdirectory of your Documents folder.
You can either change the location of your Night-time_lights_Electrification-main folder, 
or close the Command prompt window, relaunch it and type the following commands:

	cd C:\
	conda activate geo37
	jupyter notebook

You should now be able to navigate to the jupyter notebook Night-time_lights_Electrification-main folder.

The folder contains 3 Jupyter Notebook files (.ipnyb) and one Python file for secondary functions (.py).

Open 1-Google_Earth_Engine_processing by clicking on it.

You can navigate through the code and execute selected cells by pressing shift+enter.

You should be able to run the first lines of code in Jupyter Notebook. This should open a new tab to connect to your Google account.
Once connected, copy the code which is displayed on the page, and paste it in the Jupyter Notebook, in the verification code box.
If you never connected your Google account to Google Earth Engine, this should raise an error.
Go on https://earthengine.google.com/signup/ to sign-up.
Once signed-up, rerun the lines and repaste the verification code.

Running entirely 1-Google_Earth_Engine_processing should launch the processing of the raw NTL images,
which will be downloaded directly on your computer in the right folder, or on your Google Drive account if they are too big.
You can chose the countries (through their 3 letters ISO code) and the years which will be processed.
Folders are automatically created in the working document.
Geometries of certain countries are already available. If you want to create another specific area, enter the GPS coordinates as demonstrated.
Resolution (in meters/pixel) can be changed, as well as aggregation method to combine monthly images.

2-Population_data_reprojection downloads and reprojects population datasets (HRSL and WorldPop) that will be used later to correct the noisiness and level shifts of NTL images in time.

3-NTL_inter-annual_correction corrects NTL images by withdrawing the average radiance value of uninhabited areas in the surroundings.
A visualization tool is suggested at the end.



