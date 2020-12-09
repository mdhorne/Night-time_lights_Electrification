import ee
import os
import geetools
from pathlib import Path
import shutil
import shapely

def create_geometry_from_file(files):
    geometry_dict = {}
    for file in files:
        region = file[0]
        country = region[0:3]
        aoi = gpd.read_file(Path(file[1]))
        geometry_buffer = aoi.geometry[0]
        while 'MultiPolygon' in str(type(geometry_buffer)):
            geometry_buffer = geometry_buffer.buffer(0.5)
        
        geometry = [[coord[0],coord[1]] for coord in geometry_buffer.exterior.coords]
        geometry_dict[region] = ee.Geometry.Polygon([geometry])
        
        if Path("./"+country+'/misc_data/AOI_'+region+'.gpkg') in [f for f in Path("./"+country+'/misc_data').iterdir()]:
            os.remove(Path("./"+country+'/misc_data/AOI_'+region+'.gpkg'))
        aoi.to_file(Path("./"+country+'/misc_data/AOI_'+region+'.gpkg'), layer='gadm36_'+region+'_0', driver="GPKG")
        
        #for geo_file in Path(file[1]).parent.iterdir():
        #    if geo_file.stem == Path(file[1]).stem:
        #        shutil.move(geo_file, Path('./'+file[0]+'/misc_data/'+geo_file.stem))
        print("Geometry added for", region, '\n')
    return geometry_dict

def create_geometry_from_coordinates(coords):
    geometry_dict = {}
    for coord in coords:
        region = coord[0]
        country = region[0:3]
        geometry_dict[region] = ee.Geometry.Polygon(coord[1])
        poly = shapely.geometry.Polygon(coord[1])
        geodf = {'Name': [region], 'geometry': [poly]}
        geodf = gpd.GeoDataFrame(geodf, crs="EPSG:4326")
        if Path("./"+country+'/misc_data/AOI_'+region+'.gpkg') in [file for file in Path("./"+country+'/misc_data').iterdir()]:
            os.remove(Path("./"+country+'/misc_data/AOI_'+region+'.gpkg'))
        geodf.to_file(Path("./"+country+'/misc_data/AOI_'+region+'.gpkg'), layer='gadm36_'+region+'_0', driver="GPKG")  
        print("Geometry added for", region, '\n')
    return geometry_dict

def create_geometry_from_ISO(countries):
    geometry_dict = {}
    for country in countries:
        if Path('./'+country+'/misc_data/AOI_'+country+'.gpkg') not in [file for file in Path('./'+country+'/misc_data').iterdir()]:
            print('Downloading gadm36_'+country+'_gpkg.zip', end='\r')
            url = 'https://biogeo.ucdavis.edu/data/gadm3.6/gpkg/gadm36_'+country+'_gpkg.zip'
            r = requests.get(url, allow_redirects=True)
            open('./'+country+'/misc_data/gadm36_'+country+'_gpkg.zip', 'wb').write(r.content)
            zipfile.ZipFile(Path('./'+country+'/misc_data/gadm36_'+country+'_gpkg.zip'), 'r').extractall(Path('./'+country+'/misc_data'))
            os.remove(Path('./'+country+'/misc_data/gadm36_'+country+'_gpkg.zip'))
            os.remove(Path('./'+country+'/misc_data/license.txt'))
            os.rename(Path('./'+country+'/misc_data/gadm36_'+country+'.gpkg'),Path('./'+country+'/misc_data/AOI_'+country+'.gpkg'))
            print('gadm36_'+country+'.gpkg download completed')
        else:
            print("Loading AOI_"+country+".gpkg")
        aoi_in = Path("./"+country+"/misc_data/AOI_"+country+".gpkg")
        aoi = gpd.read_file(aoi_in, layer='gadm36_'+country+'_0')
        
        geometry_buffer = aoi.geometry[0]
        while 'MultiPolygon' in str(type(geometry_buffer)):
            geometry_buffer = geometry_buffer.buffer(0.5)
        
        geometry = [[coord[0],coord[1]] for coord in geometry_buffer.exterior.coords]
        geometry_dict[country] = ee.Geometry.Polygon([geometry])
        print("Geometry added for", country, '\n')
    return geometry_dict
    

def create_directories(region_names):
    print('Following countries will be processed:')
    print(region_names)

    country_names = set([region[0:3] for region in region_names])
    
    # Create directories for each country
    for country in country_names:
        if Path('./'+country) not in Path('./').iterdir():
            os.mkdir('./'+country)
        if Path('./'+country+'/misc_data') not in Path('./'+country).iterdir():
            os.mkdir('./'+country+'/misc_data')
        if Path('./'+country+'/worldpop') not in Path('./'+country).iterdir():
            os.mkdir('./'+country+'/worldpop')
        if Path('./'+country+'/NTL') not in Path('./'+country).iterdir():
            os.mkdir('./'+country+'/NTL')
    print('')
            
def add_geometry(name, coordinates, geometry_dict):
    geometry_dict[name] = ee.Geometry.Polygon(
        [
        [[coord[0], coord[1]] for coord in coordinates]
        ])
    return geometry_dict
            
def launch_gee_processes(geometry_dict, years, ntl_data='VCMSLCFG', mask_cld_value=0, reducer='mean', resolution=480, monthly=False):

    # Create a mask for NTL data. mask_cld_value sets the required number of clear sky observations.
    def mask_NTL(image):
        rad = image.select("avg_rad") # Average radiance value
        cld = image.select("cf_cvg") # Number of clear (cloud-free) observations per pixel
        cld = cld.expression('CLD > '+str(mask_cld_value), {'CLD': cld.select('cf_cvg')})
        return rad.mask(cld)
    
    regions_drive = []
    # Process the region's NTL data through Google Earth Engine
    for region in geometry_dict.keys():
        geometry = geometry_dict[region]
        country = region[0:3]

        for year in years:

            name = 'NTL_'+str(year)+'_'+region
            path = '.\\'+country+'\\NTL'

            if Path('./'+country+'/NTL/'+name+'.tif') not in [file for file in Path('./'+country+'/NTL').iterdir()] and region not in regions_drive:

                # For NTL data, chose between VCMCFG data and VCMSLCFG data. VCMSLCFG are corrected from stray light.
                # Monthly images are aggregated by their mean value, try median()/max()/percentile([70])... at line 21.

                NTL_image = ee.ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/"+ntl_data).filter(ee.Filter.calendarRange(year, year, 'year')).map(mask_NTL)
                
                if reducer == 'mean':
                    NTL_image = NTL_image.reduce(ee.Reducer.mean())
                elif reducer == None:
                    NTL_image = NTL_image.reduce(ee.Reducer.mean())
                elif reducer == 'max':
                    NTL_image = NTL_image.reduce(ee.Reducer.max())
                elif reducer == 'min':
                    NTL_image = NTL_image.reduce(ee.Reducer.min())
                elif reducer == 'median':
                    NTL_image = NTL_image.reduce(ee.Reducer.median())
                elif 'percentile' in reducer:
                    NTL_image = NTL_image.reduce(ee.Reducer.percentile([int(reducer[12:-2])]))
                
                NTL_image = NTL_image.clip(geometry).set('year', year)
                
                print('Downloading '+name+'\t\t',end='\r')

                try:
                    task = geetools.batch.image.toLocal(NTL_image,
                                                        name = name, # Path of the output
                                                        path = path, # Name of the output
                                                        scale = resolution) # Resolution in meters per pixel

                # NTL images are downloaded directly in the proper folder
                # Temporary download files might appear in the working directory            

                    print('Download completed '+name+'\t\t')
                    os.remove(Path('./'+country+'/NTL/'+name+'.zip'))
                    NTL_file = [file for file in Path('./'+country+'/NTL/'+name+'').iterdir()][0]
                    shutil.move(NTL_file, Path('./'+country+'/NTL/'+name+'.tif'))
                    os.rmdir(Path('./'+country+'/NTL/'+name))

                except:
                    print(name+' too heavy for direct download, file will be sent to Google Drive and available in a few minutes.')
                    regions_drive.append(country)

            elif country not in regions_drive:
                print(name+'.tif already in folder')

    if len(regions_drive)>0:
        print('')
        print('Sending heavy files on Google Drive. They will be available in a few minutes.')

        for region in regions_drive:
            geometry = geometry_dict[region]
            country = region[0:3]
            
            for year in years:
                ee_year = ee.List.sequence(year, year)
                                
                if reducer == 'mean':
                    NTL_proc = ee.ImageCollection.fromImages(ee_year.map(lambda y : ee.ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/"+ntl_data).filter(ee.Filter.calendarRange(y, y, 'year')).map(mask_NTL).reduce(ee.Reducer.mean()).clip(geometry).set('year', year)))
                elif reducer == None:
                    NTL_proc = ee.ImageCollection.fromImages(ee_year.map(lambda y : ee.ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/"+ntl_data).filter(ee.Filter.calendarRange(y, y, 'year')).map(mask_NTL).reduce(ee.Reducer.mean()).clip(geometry).set('year', year)))
                elif reducer == 'max':
                    NTL_proc = ee.ImageCollection.fromImages(ee_year.map(lambda y : ee.ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/"+ntl_data).filter(ee.Filter.calendarRange(y, y, 'year')).map(mask_NTL).reduce(ee.Reducer.max()).clip(geometry).set('year', year)))
                elif reducer == 'min':
                    NTL_proc = ee.ImageCollection.fromImages(ee_year.map(lambda y : ee.ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/"+ntl_data).filter(ee.Filter.calendarRange(y, y, 'year')).map(mask_NTL).reduce(ee.Reducer.min()).clip(geometry).set('year', year)))
                elif reducer == 'median':
                    NTL_proc = ee.ImageCollection.fromImages(ee_year.map(lambda y : ee.ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/"+ntl_data).filter(ee.Filter.calendarRange(y, y, 'year')).map(mask_NTL).reduce(ee.Reducer.median()).clip(geometry).set('year', year)))
                elif 'percentile' in reducer:
                    NTL_proc = ee.ImageCollection.fromImages(ee_year.map(lambda y : ee.ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/"+ntl_data).filter(ee.Filter.calendarRange(y, y, 'year')).map(mask_NTL).reduce(ee.Reducer.percentile([int(reducer[12:-2])])).clip(geometry).set('year', year)))                                                            

                task = geetools.batch.Export.imagecollection.toDrive(collection=NTL_proc,
                                                                     region=geometry,
                                                                     folder=country,
                                                                     namePattern='NTL_'+str(year)+"_"+region, # Name of the output
                                                                     scale=resolution, # Resolution in meters per pixel
                                                                     crs='EPSG:4326') # Spatial reference system
                #print(task[0].status()['description'],task[0].status()['state'], 'file should be available soon on Google Drive')

        print('')
        print('Some files were too heavy for direct download:',regions_drive)
        print('Please download the images on Google Drive and put them into the appropriate folder:')
        print('\t./COUNTRY-ISO-CODE/NTL')

    print('')
    print('Done')

####

from pathlib import Path
import shutil
import requests
import zipfile
import rasterio
import os
import numpy as np
from rasterio.warp import reproject, Resampling
import geopandas as gpd
from raster_manipulation import save_raster, clip_raster, exact_raster_reprojection, fast_raster_reprojection


# Ignore warnings from numpy (division by 0 or nan)
_ = np.seterr(divide='ignore', invalid='ignore')





def pop_data_download(region_names, wp_year=2017):
    
    from hdx.utilities.easy_logging import setup_logging
    setup_logging()
    from hdx.hdx_configuration import Configuration
    Configuration.create(hdx_site='prod', user_agent='Read-only user', hdx_read_only=True)
    from hdx.data.dataset import Dataset
    
    import wpgpDownload
    from wpgpDownload.utils.convenience_functions import download_country_covariates as download_worldpop
    from wpgpDownload.utils.convenience_functions import refresh_csv
    refresh_csv()

    hdx_datasets = Dataset.search_in_hdx('hrsl', rows=500)
    hdx_resources = Dataset.get_all_resources(hdx_datasets)
    
    print('')

    country_names = set([region[0:3] for region in region_names])

    for country in country_names:
        print(country)

        for res in hdx_resources:
            if 'population_'+country.lower() in res['name'] and '.zip' in res['name'] and 'csv' not in res['name']:
                print('Downloading HRSL',res['name'], end='\r')
                url, path = res.download()
                print('HRSL',res['name'],'download completed       ')
                shutil.move(Path(path),Path('./'+country+'/misc_data/population_'+country.lower()+'.zip'))
                zipfile.ZipFile(Path('./'+country+'/misc_data/population_'+country.lower()+'.zip'), 'r').extractall(Path('./'+country+'/misc_data'))
                for file in Path('./'+country+'/misc_data').iterdir():
                    if 'population_'+country.lower() in file.name and file.suffix != '.tif':
                        os.remove(file)
        
        if type(wp_year) == list:
            years = wp_year
        elif type(wp_year) == int: 
            years = [wp_year]

        #NTL_files = [file for file in Path("./"+country+"/NTL").iterdir() if "NTL" in file.name]
        #
        #years = []
        #for NTL_file in NTL_files:
        #    years.append(NTL_file.name[4:8])
        #years = [year for year in set(years)]
        #years.sort()

        for year in years:
            print('Downloading WorldPop '+country+' '+str(year)+'\t\t',end='\r')
            download_worldpop(ISO=country,out_folder='.\\'+country+'\\worldpop',prod_name='ppp_'+str(year))
            print('WorldPop '+country+' '+str(year)+' download completed\t\t')
        
        print("")
        
    print('Done')
    
def clipping_rasters(region_names):
    for region in region_names:
        country = region[0:3]
        aoi_name = "AOI_"+region+".gpkg"
        aoi_in = Path("./"+country+"/misc_data") / aoi_name
        aoi = gpd.read_file(aoi_in, layer='gadm36_'+region+'_0')
        for NTL_file in Path("./"+country+"/NTL").iterdir():
            if region == NTL_file.stem[9:]:
                clipped, aff, crs = clip_raster(NTL_file, aoi)
                save_raster(NTL_file, clipped, aff)
                print(NTL_file.name, 'clipped')
    
def reproject_rasters(region_names):
    
    for region in region_names:
        country=region[0:3]
        print('')
        print(region)
        
        path_in = Path("./"+region)
        
        aoi_name = "AOI_"+region+".gpkg"
        aoi_in = Path("./"+country+"/misc_data") / aoi_name
        aoi = gpd.read_file(aoi_in, layer='gadm36_'+region+'_0')
        
        # Opening first NTL raster to get the targeted resolution
        
        NTL_file = [file for file in Path("./"+country+"/NTL").iterdir() if region == file.stem[9:]][0]
        
        ntl = rasterio.open(NTL_file)
        
        # Reprojecting HRSL raster

        misc_data_files = [file.name for file in Path("./"+country+"/misc_data").iterdir()]
        if ("HRSL_"+region+"_NTL-resolution.tif" not in misc_data_files):
            if len([file for file in misc_data_files if 'population_'+country.lower() in file])==1:
                hrsl_name = [file for file in misc_data_files if 'population_'+country.lower() in file][0]
                hrsl_in = Path("./"+country+"/misc_data") / hrsl_name
                hrsl = rasterio.open(hrsl_in)
                
                print('Reprojecting ', hrsl_name)
                hrsl_ntl_reso, _ = exact_raster_reprojection(hrsl, ntl)
                
                hrsl_out = "HRSL_"+region+"_NTL-resolution.tif"
                path_out = Path("./"+country+"/misc_data") / hrsl_out
                save_raster(path_out, hrsl_ntl_reso, ntl.transform)
                clipped, aff, _ = clip_raster(path_out, aoi)
                save_raster(path_out, clipped, aff)
                
            else:
                print('Unable to find one HRSL file')
                    
        # Reprojecting WorldPop rasters
        
        worldpop_files = [file for file in Path("./"+country+"/worldpop").iterdir() if 'ppp' in file.name]
        if len(worldpop_files) > 0:
            for file in worldpop_files:
                if "worldpop_"+file.stem[-4:]+'_'+region+"_NTL-resolution.tif" not in [file.name for file in Path("./"+country+"/worldpop").iterdir()]:
                    worldpop_in = Path("./"+country+"/worldpop") / file.name
                    worldpop = rasterio.open(worldpop_in)

                    print('Reprojecting ',file.name)
                    #worldpop_ntl_reso, _ = exact_raster_reprojection(worldpop, ntl) # For exact but slow reprojection
                    worldpop_ntl_reso = fast_raster_reprojection(worldpop, ntl)

                    ntl_out = "worldpop_"+file.stem[-4:]+'_'+region+"_NTL-resolution.tif"
                    path_out = Path("./"+country+"/worldpop") / ntl_out
                    save_raster(path_out, worldpop_ntl_reso, ntl.transform)
                    clipped, aff, _ = clip_raster(path_out, aoi)
                    save_raster(path_out, clipped, aff)
        else:
            print('Unable to find a WorldPop file')
            
    print('')
    print('Done')
    
def interannual_correction(region_names, population_threshold=20):

    for region in region_names:
        country = region[0:3]
        
        NTL_files = [file for file in Path("./"+country+"/NTL").iterdir()
                    if "local-correction" not in file.name and 'NTL' in file.name
                    and region == file.stem[9:]]
        
        misc_data_files = [file.name for file in Path("./"+country+"/misc_data").iterdir()]
        
        # Opening area of interest (AOI)
        
        aoi_name = "AOI_"+region+".gpkg"
        aoi_in = Path("./"+country+"/misc_data") / aoi_name
        aoi = gpd.read_file(aoi_in, layer='gadm36_'+region+'_0')
        
        for NTL_file in NTL_files:
            
            # Opening the population raster to identify unhabited areas
            
            worldpop_name = 'worldpop_2017_'+region+'_NTL-resolution.tif' # WorldPop data from the exact same year is not necessary. 2017 is used here.
            worldpop_in = Path("./"+country+"/worldpop") / worldpop_name
            worldpop = rasterio.open(worldpop_in).read(1)
            worldpop = (worldpop>population_threshold)*worldpop
            
            if "HRSL_"+region+"_NTL-resolution.tif" in misc_data_files:
                raster_pop_name = 'HRSL_'+region+'_NTL-resolution.tif'
                raster_pop_in = Path("./"+country+"/misc_data") / raster_pop_name
                raster_pop = rasterio.open(raster_pop_in).read(1)
            else:
                raster_pop = worldpop
                
            # Opening NTL raster
            
            ntl = rasterio.open(NTL_file)
            ntl_rd = ntl.read(1)
            
            # Replacing 0 values, i.e where no data is available, by np.nan
            
            ntl_nodata = 1.0*(rasterio.open(NTL_file).read_masks()==255)[0]
            ntl_nodata[ntl_nodata == 0.] = np.nan
            
            # Creating 2 masks, one for populated areas and another for unpopulated areas
            
            mask_ntl_pop = np.full_like(ntl_rd, np.nan)
            mask_ntl_pop[raster_pop>0] = 1
            mask_ntl_unpop = np.full_like(ntl_rd, np.nan)
            mask_ntl_unpop[raster_pop<=0] = 1
            mask_ntl_unpop[np.isnan(raster_pop)] = 1
            mask_ntl_unpop = mask_ntl_unpop*ntl_nodata

            ntl_rd_pop = ntl_rd*mask_ntl_pop
            ntl_rd_unpop = ntl_rd*mask_ntl_unpop
            ntl_corr = np.full_like(ntl_rd, np.nan)
            
            # Beginning the correction process
            
            # Each NTL value is corrected based on the following principle:
            # Look for unpopulated cells in the surroundings of the cell to be corrected
            # Surroundings is determined by a (2 x diameter) x 2 x diameter) square around the cell to be corrected
            # >250 unpopulated cells should be found in the surroundings, otherwise, diameter increases until >250 unpopulated cells are found
            # The lowest and highest deciles of the >250 NTL values are removed
            # Withdraw of the mean value of radiance of the remaining unpopulated cells from the cell being corrected
            # 
            # In short terms, NTL values are corrected by withdrawing the average value of radiance of unpopulated cells closeby
            # Unpopulated cells are corrected by considering they should be dark (radiance is naught).
            
            diameter = 15

            for row in range(ntl_corr.shape[0]):
                for col in range(ntl_corr.shape[1]):
                    if ntl_nodata[row, col] == 1:
                        diameter_temp = diameter
                        unpop_surroundings = ntl_rd_unpop[max(0,row-diameter_temp):min(ntl_rd_pop.shape[0],row+diameter_temp+1), max(0,col-diameter_temp):min(ntl_rd_pop.shape[1],col+diameter_temp+1)]
                        while np.count_nonzero(~np.isnan(unpop_surroundings))<250:
                            diameter_temp += 1
                            unpop_surroundings = ntl_rd_unpop[max(0,row-diameter_temp):min(ntl_rd_pop.shape[0],row+diameter_temp+1), max(0,col-diameter_temp):min(ntl_rd_pop.shape[1],col+diameter_temp+1)]

                        list_unpop_surroundings = unpop_surroundings[np.logical_not(np.isnan(unpop_surroundings))]
                        list_unpop_surroundings.sort()
                        list_unpop_surroundings = list_unpop_surroundings[int(len(list_unpop_surroundings)/10):-int(len(list_unpop_surroundings)/10)]
                        ntl_corr[row,col] = ntl_rd[row,col]-list_unpop_surroundings.mean()

                if row/100 == np.floor(row/100):
                    print('Correcting',NTL_file.stem, int(row/ntl_rd_pop.shape[0]*100),'%      ',end="\r")
            
            # If radiance is found negative, radiance is set to 0
            
            ntl_corr[ntl_corr<0] = 0
            
            # Radiance of unpopulated cells is set to 0
            
            ntl_corr[raster_pop<=5] = 0
            
            # Saving corrected raster
            
            ntl_out = NTL_file.stem+'_local-correction'+NTL_file.suffix
            path_out = Path("./"+country+"/NTL/Corrected") / ntl_out
            save_raster(path_out, ntl_corr, ntl.transform)
            clipped, aff, crs = clip_raster(path_out, aoi)
            save_raster(path_out, clipped, aff)

    print('')
    print('')
    print('Done')
    
def correction_plot(region_names):

    import matplotlib.pyplot as plt
    for region in region_names:
        country = region[0:3]

        NTL_files = [file for file in Path("./"+country+"/NTL/Corrected").iterdir()
                if "local-correction" in file.name]
        
        misc_data_files = [file.name for file in Path("./"+country+"/misc_data").iterdir()]
        
        years = []
        for NTL_file in NTL_files:
            years.append(NTL_file.name[4:8])
        years = [year for year in set(years)]
        years.sort()

        results_all = []
        results_unpop = []
        results_pop = []
        
        results_all_corr = []
        results_unpop_corr = []
        results_pop_corr = []
        
        for year in years:
        
            # Opening the population raster to identify unhabited areas
            
            population_threshold = 20 # Population threshold for WorldPop rasters to distinguish inhabited/uninhabited areas
            worldpop_name = 'worldpop_2017_'+region+'_NTL-resolution.tif' # WorldPop data from the exact same year is not necessary. 2017 is used here.
            worldpop_in = Path("./"+country+"/worldpop") / worldpop_name
            worldpop = rasterio.open(worldpop_in).read(1)
            worldpop = (worldpop>population_threshold)*worldpop

            if "HRSL_"+region+"_NTL-resolution.tif" in misc_data_files:
                raster_pop_name = 'HRSL_'+region+'_NTL-resolution.tif'
                raster_pop_in = Path("./"+country+"/misc_data") / raster_pop_name
                raster_pop = rasterio.open(raster_pop_in).read(1)
            else:
                raster_pop = worldpop

            # Opening NTL rasters
            
            ntl_name = "NTL_"+year+"_"+region+".tif"
            ntl_in = Path("./"+country+'/NTL') / ntl_name
            ntl = rasterio.open(ntl_in)
            ntl_rd = ntl.read(1)
            
            ntl_corr_name = "NTL_"+year+"_"+region+"_local-correction.tif"
            ntl_corr_in = Path("./"+country+'/NTL/Corrected') / ntl_corr_name
            ntl_corr = rasterio.open(ntl_corr_in)
            ntl_corr_rd = ntl_corr.read(1)

            # Replacing 0 values, i.e where no data is available, by np.nan

            ntl_nodata = 1.0*(rasterio.open(ntl_in).read_masks()==255)[0]
            ntl_nodata[ntl_nodata == 0.] = np.nan

            # Creating 2 masks, one for populated areas and another for unpopulated areas

            mask_ntl_pop = np.full_like(ntl_rd, np.nan)
            mask_ntl_pop[raster_pop>5] = 1
            mask_ntl_pop = mask_ntl_pop*ntl_nodata
            mask_ntl_unpop = np.full_like(ntl_rd, np.nan)
            mask_ntl_unpop[raster_pop<=5] = 1
            mask_ntl_unpop[np.isnan(raster_pop)] = 1
            mask_ntl_unpop = mask_ntl_unpop*ntl_nodata
            
            ntl_rd = ntl_rd*ntl_nodata
            ntl_rd_pop = ntl_rd*mask_ntl_pop
            ntl_rd_unpop = ntl_rd*mask_ntl_unpop
        
            results_all.append(np.nanmean(ntl_rd))
            results_unpop.append(np.nanmean(ntl_rd_unpop))
            results_pop.append(np.nanmean(ntl_rd_pop))

            ntl_corr_rd_pop = ntl_corr_rd*mask_ntl_pop
            ntl_corr_rd_unpop = ntl_corr_rd*mask_ntl_unpop
        
            results_all_corr.append(np.nanmean(ntl_corr_rd))
            results_unpop_corr.append(np.nanmean(ntl_corr_rd_unpop))
            results_pop_corr.append(np.nanmean(ntl_corr_rd_pop))

        plt.figure(figsize=(15,5))
        plt.suptitle(region)
        plt.subplot(121)
        plt.title('Uncorrected')
        plt.plot([int(year) for year in years], results_pop, '^-')
        plt.plot([int(year) for year in years], results_all, 'o-')
        plt.plot([int(year) for year in years], results_unpop, 'v-')
        plt.legend(['Populated areas','All areas','Unpopulated areas'])
        plt.ylabel('Mean of radiance')

        plt.subplot(122)
        plt.title('Corrected')
        plt.plot([int(year) for year in years], results_pop_corr, '^-')
        plt.plot([int(year) for year in years], results_all_corr, 'o-')
        plt.plot([int(year) for year in years], results_unpop_corr, 'v-')
        plt.legend(['Populated areas','All areas','Unpopulated areas'])
        plt.ylabel('Mean of radiance')
        
        plt.show()