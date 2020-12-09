from pathlib import Path
import rasterio
from rasterio.mask import mask
from rasterio.warp import reproject, Resampling
import geopandas as gpd
import numpy as np
import json

def add_value_to_raster(raster, px, py, value):
    if (px < raster.shape[0]
    and px >= 0
    and py < raster.shape[1]
    and py >= 0):
        if np.isnan(raster[px,py]):
            raster[px,py]=0
        raster[px,py] += value

def save_raster(path, raster, affine, crs=None, nodata=0):
    """Save a raster to the specified file.
    From Gridfinder and Christopher Arderne's Github
    Parameters
    ----------
    file : str
        Output file path
    raster : numpy.array
        2D numpy array containing raster values
    affine: affine.Affine
        Affine transformation for the raster
    crs: str, proj.Proj, optional (default EPSG4326)
        CRS for the raster
    """

    path = Path(path)
    if not path.parents[0].exists():
        path.parents[0].mkdir(parents=True, exist_ok=True)

    if not crs:
        crs = "+proj=latlong"

    filtered_out = rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=raster.shape[0],
        width=raster.shape[1],
        count=1,
        dtype=raster.dtype,
        crs=crs,
        transform=affine,
        nodata=nodata,
    )
    filtered_out.write(raster, 1)
    filtered_out.close()
    
def clip_raster(raster, boundary, boundary_layer=None):
    """Clip the raster according to the given administrative boundary.
    From Gridfinder and Christopher Arderne's Github
    Parameters
    ----------
    raster : string, pathlib.Path or rasterio.io.DataSetReader
        Location of or already opened raster.
    boundary : string, pathlib.Path or geopandas.GeoDataFrame
        The polygon by which to clip the raster.
    boundary_layer : string, optional
        For multi-layer files (like GeoPackage), specify the layer to be used.
    Returns
    -------
    tuple
        Three elements:
            clipped : numpy.ndarray
                Contents of clipped raster.
            affine : affine.Affine()
                Information for mapping pixel coordinates
                to a coordinate system.
            crs : dict
                Dict of the form {'init': 'epsg:4326'} defining the coordinate
                reference system of the raster.
    """

    if isinstance(raster, Path):
        raster = str(raster)
    if isinstance(raster, str):
        raster = rasterio.open(raster)

    if isinstance(boundary, Path):
        boundary = str(boundary)
    if isinstance(boundary, str):
        if ".gpkg" in boundary:
            driver = "GPKG"
        else:
            driver = None  # default to shapefile
            boundary_layer = ""  # because shapefiles have no layers

        boundary = gpd.read_file(boundary, layer=boundary_layer, driver=driver)

    if not (boundary.crs == raster.crs or boundary.crs == raster.crs.data):
        boundary = boundary.to_crs(crs=raster.crs)
    coords = [json.loads(boundary.to_json())["features"][0]["geometry"]]

    # mask/clip the raster using rasterio.mask
    clipped, affine = mask(dataset=raster, shapes=coords, crop=True)

    if len(clipped.shape) >= 3:
        clipped = clipped[0]

    return clipped, affine, raster.crs

def exact_raster_reprojection(pop, ntl):
    """Reproject and downsample a (population) raster exactly according to the resolution of another source (night-time light image).
    Parameters
    ----------
    pop : (Population) raster to upsample, opened with Rasterio
    ntl : (Night-time light) raster with the right resolution, opened with Rasterio
    
    Returns as a first output a numpy.array containing the data from pop_rd downsampled and reprojected to match ntl raster's resolution.
    Second output is a numpy.array containing the percentage of the area covered with buildings. This output is meaningful only if the population
    data is based on footprint recognition (typically HRSL).
    Trying to read big tif files might lead to a memory error.
    """
    
    try:
        pop_rd = pop.read(1)
    except MemoryError:
        print("Unable to open file due to a lack of available memory.")
        print("Trying to open file row by row (float32).")
        pop_rd = np.zeros((pop.shape[0],pop.shape[1]), dtype=np.float32)
        from rasterio.windows import Window
        for row in range(pop.shape[0]):
            w = Window(0, row, pop.shape[1]+1, 1)
            r = pop.read(1, window=w)
            pop_rd[row,:] = r
    
    #Creating the final array that will contain the data from the source
    raster_pop_ntl = np.full_like(ntl.read(1), np.nan)
    rester_built_area = np.full_like(ntl.read(1), np.nan)
    
    # Calculating the area of one HRSL cell    
    hrsl_area_build = np.abs(pop.transform[0]*pop.transform[4])

    #Calculating the transformation ratios
    percent_x = pop.transform[0]/ntl.transform[0]
    percent_y = pop.transform[4]/ntl.transform[4]

    percent_0_x = (pop.xy(0, 0, offset='ul')[0] - ntl.xy(0,0, offset='ul')[0])/ntl.transform[0]
    percent_0_y = (pop.xy(0, 0, offset='ul')[1] - ntl.xy(0,0, offset='ul')[1])/ntl.transform[0]

    for row_pop in range(pop.shape[0]):
        for col_pop in range(pop.shape[1]):
            if pop_rd[row_pop][col_pop] > 0:
                row_ntl = row_pop*percent_y-percent_0_y
                col_ntl = col_pop*percent_y+percent_0_x
                row_ind = np.int(np.floor(row_ntl))
                col_ind = np.int(np.floor(col_ntl))
                
                #If pop cell entirely within ntl cell
                
                if (col_ntl-np.floor(col_ntl)<1-percent_x
                and row_ntl-np.floor(row_ntl)<1-percent_y):
                    if (row_ind < raster_pop_ntl.shape[0]
                    and row_ind >= 0
                    and col_ind < raster_pop_ntl.shape[1]
                    and col_ind >= 0):
                        add_value_to_raster(raster_pop_ntl, row_ind, col_ind, pop_rd[row_pop][col_pop])
                        add_value_to_raster(rester_built_area, row_ind, col_ind, hrsl_area_build)
                        
                #If pop cell astride two ntl cells next to each other on the x-axis (left-right)

                elif (col_ntl-np.floor(col_ntl)>1-percent_x
                and row_ntl-np.floor(row_ntl)<1-percent_y):
                    
                    #Allocating pop to the left cell
                    
                    if (row_ind < raster_pop_ntl.shape[0]
                    and row_ind >= 0
                    and col_ind < raster_pop_ntl.shape[1]
                    and col_ind >= 0):
                        add_value_to_raster(raster_pop_ntl, row_ind, col_ind, (1-(col_ntl-np.floor(col_ntl)))/percent_x*pop_rd[row_pop][col_pop])
                        add_value_to_raster(rester_built_area, row_ind, col_ind, (1-(col_ntl-np.floor(col_ntl)))/percent_x*hrsl_area_build)
                        
                    #Allocating pop to the right cell
                        
                    if (row_ind < raster_pop_ntl.shape[0]
                    and row_ind >= 0
                    and col_ind+1 < raster_pop_ntl.shape[1]
                    and col_ind+1 >= 0):
                        add_value_to_raster(raster_pop_ntl, row_ind, col_ind+1, (1-(1-(col_ntl-np.floor(col_ntl)))/percent_x)*pop_rd[row_pop][col_pop])
                        add_value_to_raster(rester_built_area, row_ind, col_ind+1, (1-(1-(col_ntl-np.floor(col_ntl)))/percent_x)*hrsl_area_build)
                        
                #If pop cell astride two ntl cells next to each other on the y-axis (bottom-top)
                
                elif (col_ntl-np.floor(col_ntl)<1-percent_x
                and row_ntl-np.floor(row_ntl)>1-percent_y):
                    
                    #Allocating pop to the upper cell
                    
                    if (row_ind < raster_pop_ntl.shape[0]
                    and row_ind >= 0
                    and col_ind < raster_pop_ntl.shape[1]
                    and col_ind >= 0):
                        add_value_to_raster(raster_pop_ntl, row_ind, col_ind, (1-(row_ntl-np.floor(row_ntl)))/percent_y*pop_rd[row_pop][col_pop])
                        add_value_to_raster(rester_built_area, row_ind, col_ind, (1-(row_ntl-np.floor(row_ntl)))/percent_y*hrsl_area_build)
        
                    #Allocating pop to the lower cell
                    
                    if (row_ind+1 < raster_pop_ntl.shape[0]
                    and row_ind+1 >= 0
                    and col_ind < raster_pop_ntl.shape[1]
                    and col_ind >= 0):
                        add_value_to_raster(raster_pop_ntl, row_ind+1, col_ind, (1-(1-(row_ntl-np.floor(row_ntl)))/percent_y)*pop_rd[row_pop][col_pop])
                        add_value_to_raster(rester_built_area, row_ind+1, col_ind, (1-(1-(row_ntl-np.floor(row_ntl)))/percent_y)*hrsl_area_build)
                        
                #If pop cell astride four ntl cells
                
                elif (col_ntl-np.floor(col_ntl)>1-percent_x
                and row_ntl-np.floor(row_ntl)>1-percent_y):
                    
                    #Allocating pop to the upper left cell
                    
                    if (row_ind < raster_pop_ntl.shape[0]
                    and row_ind >= 0
                    and col_ind < raster_pop_ntl.shape[1]
                    and col_ind >= 0):
                        add_value_to_raster(raster_pop_ntl, row_ind, col_ind, (1-(col_ntl-np.floor(col_ntl)))/percent_x*(1-(row_ntl-np.floor(row_ntl)))/percent_y*pop_rd[row_pop][col_pop])
                        add_value_to_raster(rester_built_area, row_ind, col_ind, (1-(col_ntl-np.floor(col_ntl)))/percent_x*(1-(row_ntl-np.floor(row_ntl)))/percent_y*hrsl_area_build)

                    #Allocating pop to the upper right cell
                        
                    if (row_ind < raster_pop_ntl.shape[0]
                    and row_ind >= 0
                    and col_ind+1 < raster_pop_ntl.shape[1]
                    and col_ind+1 >= 0):
                        add_value_to_raster(raster_pop_ntl, row_ind, col_ind+1, (1-(1-(col_ntl-np.floor(col_ntl)))/percent_x)*(1-(row_ntl-np.floor(row_ntl)))/percent_y*pop_rd[row_pop][col_pop])
                        add_value_to_raster(rester_built_area, row_ind, col_ind+1, (1-(1-(col_ntl-np.floor(col_ntl)))/percent_x)*(1-(row_ntl-np.floor(row_ntl)))/percent_y*hrsl_area_build)

                    #Allocating pop to the lower left cell
                        
                    if (row_ind+1 < raster_pop_ntl.shape[0]
                    and row_ind+1 >= 0
                    and col_ind < raster_pop_ntl.shape[1]
                    and col_ind >= 0):
                        add_value_to_raster(raster_pop_ntl, row_ind+1, col_ind, (1-(col_ntl-np.floor(col_ntl)))/percent_x*(1-(1-(row_ntl-np.floor(row_ntl)))/percent_y)*pop_rd[row_pop][col_pop])
                        add_value_to_raster(rester_built_area, row_ind+1, col_ind, (1-(col_ntl-np.floor(col_ntl)))/percent_x*(1-(1-(row_ntl-np.floor(row_ntl)))/percent_y)*hrsl_area_build)

                    #Allocating pop to the lower right cell
                    
                    if (row_ind+1 < raster_pop_ntl.shape[0]
                    and row_ind+1 >= 0
                    and col_ind+1 < raster_pop_ntl.shape[1]
                    and col_ind+1 >= 0):
                        add_value_to_raster(raster_pop_ntl, row_ind+1, col_ind+1, (1-(1-(col_ntl-np.floor(col_ntl)))/percent_x)*(1-(1-(row_ntl-np.floor(row_ntl)))/percent_y)*pop_rd[row_pop][col_pop])
                        add_value_to_raster(rester_built_area, row_ind+1, col_ind+1, (1-(1-(col_ntl-np.floor(col_ntl)))/percent_x)*(1-(1-(row_ntl-np.floor(row_ntl)))/percent_y)*hrsl_area_build)

        #Printing percentage of processed cells
                        
        if row_pop/100 == np.floor(row_pop/100):
            print(int(row_pop/pop.shape[0]*100),'%      ',end="\r")
    
    #Dividing by the unitary area of HRSL cells
    
    rester_built_area = rester_built_area/np.abs(ntl.transform[0]*ntl.transform[4])
    raster_pop_ntl = np.nansum(pop_rd)/np.nansum(raster_pop_ntl)*raster_pop_ntl
    
    print('100%      ',end="\r")
    
    return raster_pop_ntl, rester_built_area

def fast_raster_reprojection(pop, ntl):
    """Reproject and downsample a (population) raster fastly but inexactly (interpolation) according to the resolution of another source (night-time light image).
    Parameters
    ----------
    pop : (Population) raster to upsample, opened with Rasterio
    ntl : (Night-time light) raster with the right resolution, opened with Rasterio
    
    Returns as a numpy.array containing the data from pop downsampled and reprojected to match ntl raster's resolution.
    This output is meaningful only if pop and ntl have a similar resolution (1 to 5), pop is continuous its data sufficiently smooth.
    Otherwise, interpolation might produce non-negligible errors.
    Trying to read big tif files might lead to a memory error.
    """
    
    try:
        pop_rd = pop.read(1)
    except MemoryError:
        print("Unable to open file due to a lack of available memory.")
        print("Trying to open file row by row (float32).")
        pop_rd = np.zeros((pop.shape[0],pop.shape[1]), dtype=np.float32)
        from rasterio.windows import Window
        for row in range(pop.shape[0]):
            w = Window(0, row, pop.shape[1]+1, 1)
            r = pop.read(1, window=w)
            pop_rd[row,:] = r
    
    pop_rd = (pop_rd>=0)*pop_rd
    
    src_transform = pop.transform
    src_crs = {'init': 'EPSG:4326'}

    dst_shape = ntl.shape
    dst_transform = ntl.transform
    dst_crs = {'init': 'EPSG:4326'}
    raster_pop_ntl = np.zeros(dst_shape)

    reproject(
        pop_rd,
        raster_pop_ntl,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=dst_transform,
        dst_crs=dst_crs,
        resampling=Resampling.nearest)

    raster_pop_ntl = raster_pop_ntl*ntl.transform[0]**2/pop.transform[0]**2
    raster_pop_ntl = (raster_pop_ntl>=0)*raster_pop_ntl
    raster_pop_ntl = np.nansum(pop_rd)/np.nansum(raster_pop_ntl)*raster_pop_ntl

    return raster_pop_ntl
