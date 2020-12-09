# Night-time lights to assess electrification rates in developing countries
Project funded by the World Bank Group and ESMAP to estimate electricity access rate with night-time lights and other satellite-based datasets in developing countries.

The code provided:
1) Downloads raw night-time lights data (VIIRS satellite operated by NOAA) from Google Earth Engine
2) Downloads population data (WorldPop from CIESIN and HRSL from Facebook) and reproject them to match the resolution of night-time lights data
3) Corrects the inter-annual shift in noisiness of night-time light data using an original algorithm based on uninhabited areas

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

Then run the following lines to set up the environment and type 'y' for Yes if requested:

	conda update --all	
	conda install jupyter	
	pip install earthengine-api --upgrade	
	pip install geetools

From https://www.lfd.uci.edu/~gohlke/pythonlibs/, download the following files (chose 32 or 64 bits):

	rasterio‑1.1.8‑cp37‑cp37m‑win_amd64.whl
	GDAL‑3.1.4‑cp37‑cp37m‑win_amd64.whl
	Shapely‑1.7.1‑cp37‑cp37m‑win_amd64.whl
	Fiona‑1.8.18‑cp37‑cp37m‑win_amd64.whl
	geopandas‑0.8.1‑py3‑none‑any.whl
	pyproj‑3.0.0.post1‑cp37‑cp37m‑win_amd64.whl
	Rtree‑0.9.4‑cp37‑cp37m‑win_amd64.whl

Execute:

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

Close the command prompt. The environment is now ready to be used.

## Using Night-time_lights_Electrification code

Download Night-time_lights_Electrification-main.zip from Github and unzip it a dedicated folder. 
This will be the working directory and downloaded files and outputs will appear there if possible.

Then run Jupyter Notebook by opening a command prompt (cmd) and typing:

	conda activate geo37
	jupyter notebook

This should open a new tab in your default browser.
You should be able to navigate to the Night-time_lights_Electrification-main folder through the Jupyter Notebook interface.
If not, that means you probably put the Night-time_lights_Electrification-main folder somewhere at the root of your hard drive disk
and not in subfolder of your Documents.
You can either change the location of your Night-time_lights_Electrification-main folder, 
or close the Command prompt window, relaunch it and type the following commands instead:

	cd C:\
	conda activate geo37
	jupyter notebook

The folder contains:
1) A Jupyter Notebook file, Night-time_light_correction.ipnyb, containing the main program
2) Two Python files, raster_manipulation.py and files_creation.py, for secondary functions
3) Geometry files for Belgium as an example (gadm36_BEL_0)

### Night-time_light_correction

Open Night-time_light_correction.ipynb by clicking on it through the Jupyter interface.

You can navigate through the code and execute selected cells with code by pressing 'shift + enter'.

The first cell should open a new tab to connect to your Google account.
If you never connected your Google account to Google Earth Engine, this should raise an error.
Go on https://earthengine.google.com/signup/ to sign-up and rerun the code.
Sign-in and copy the code displayed on the page, and paste it in the Jupyter Notebook, in the verification code box as explained.

Running the first part of the program should launch the processing of the raw NTL images,
which will be downloaded directly on your computer in the right folder (the working directory where you extracted Night-time_lights_Electrification-main.zip) if possible.
You can chose the countries (through their 3 letters ISO code) and the years which will be processed.
Folders will be automatically created in the working document.
Geometries of certain countries are already available. If you want to create another specific area, enter the GPS coordinates as demonstrated in the code or provide the path to an adequate file (shapefile, geopackage, json...).
Resolution (in meters/pixel) can be changed, as well as aggregation method to combine monthly images.
Large countries might require an additional step as their files are too big to be downloaded directly (>30 Mo).
If so, the files will be transfered to your Google Drive account. You will find them under the folder called accordingly to the 3 letters ISO code.
Download the files and extract them manually it the proper folder on your computer.

The second part downloads and reprojects population datasets (HRSL and WorldPop) that will be used later to correct the noisiness and level shifts of NTL images in time. By default HRSL data are downloaded from Humanitarian Data Exchange (https://data.humdata.org/) and WorldPop files from its organization website (https://www.worldpop.org/).

As WorldPop files' size is substantial, only the data of 2017 are downloaded and serve as the default population data to correct night-time lights data alongside HRSL.

Night-time lights data are clipped accordingly to the global administrative areas found on https://gadm.org/ for complete countries.

Third part corrects NTL images by withdrawing the average radiance value of uninhabited areas in the surroundings.
A visualization tool is suggested at the end to compare raw data and corrected ones.
