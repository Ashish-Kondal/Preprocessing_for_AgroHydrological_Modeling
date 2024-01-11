#=================================================================
# Converting ASCII Files of NMME data to NETCDF and Merging GridMET
# Author: Ashish Kondal
# Date: September 06, 2022
# Description: This scripts merge GridMET and NMME data to make continuous time series of 4-5 years.
#              And then save the file into NetCDF.
#=================================================================

#================================================================
#   IMPORTING LIBRARIES/PACKAGES
#================================================================

from __future__ import print_function
import sys
import os, string
import datetime as dt
import netCDF4 as nc
import pandas as pd
import numpy as np
import calendar
import itertools as itt
from tqdm import tqdm
import time

#================================================================
#   DEFINE VARIABLES AND INPUT FILES
#================================================================

input_directory = 'path_to/VIC_CropSyst_Forcings/Subsetting_SingleGrids_24thRes/';      # Path to input directory (output directory of "Subsetting_NMME_GridMET.py")
fullvariable_name = ['precipitation_amount','max_air_temperature','min_air_temperature','wind_speed','specific_humidity','surface_downwelling_shortwave_flux_in_air','max_relative_humidity','min_relative_humidity'];
basins = ['okanogon','wallawalla','yakima','pnw']

init_month = sys.argv[1];                                       # Retrieving NMME's initialization month from sbatch script
init_month = int(init_month);
f'initialization Month = {init_month}'

varnum = sys.argv[2]                                            # Retrieving "variables" number (goes from 1-8) from sbatch script
varnum = int(varnum)-1
f'Variable Number = {varnum}'

basin_id = sys.argv[3]                                          # Retrieving "basins" number (goes from 1-4; see "basins") from batch script
basin_id = int(basin_id)-1
f'basin = {basins[basin_id]}'

gridnum = sys.argv[4]                                          # Retrieving grid block number (goes from 1-160) from batch script
gridnum = int(gridnum)-1
#f'basin = {basins[gridnum]}' # Goes from 0-159


# Latitudinal and Longitudinal Extent of Study Area
lat_blocks = np.arange(39.5,50.5,1)
lon_blocks = np.arange(-125,-108,1)

# Selecting coordinating based on current "gridnum"
cord_ind = pd.DataFrame(itt.product(range(len(lat_blocks)-1), range(len(lon_blocks)-1)))
ii = cord_ind.iloc[gridnum,0]
jj = cord_ind.iloc[gridnum,1]

# Headers of input files    
nmme_header = ['Lat','Lon','Initialization_Year','Initialization_Month','Forecast_Year','Forecast_Month','Forecast_Day','Precipitation_mm', 'Max_Temperature_K', 'Min_Temperature_K', 'Wind_Speed_mps','Specific_Humidity_Kg_per_Kg','ShortWave_Radiation_W_per_m2','Max_Relative_Humidity','Min_Relative_Humidity'];
gridmet_header = ['Lat','Lon','Year','DOY','precipitation_amount','max_air_temperature','min_air_temperature','wind_speed','specific_humidity','surface_downwelling_shortwave_flux_in_air','max_relative_humidity','min_relative_humidity'];
req_cols_nmme = ['Lat','Lon','Initialization_Year','Initialization_Month','Forecast_Year','Forecast_Month','Forecast_Day',nmme_header[varnum+7]]
req_cols_gmet = ['Lat','Lon','Year','DOY',gridmet_header[varnum+4]]
print(gridmet_header[varnum+4])


#===============================================================================================
#  function to convert ASCII files into NETCDF after appending gridMET data to NMME forecast
#===============================================================================================

def asc2nc(outpath,init_month,basin_id,gmet_data,varnum,ii,jj):

    if ((init_month == 6) & (varnum <= 4)):
        nmme_filename = input_directory + 'NMME/InitMonth_0' + str(init_month) +"/" + fullvariable_name[varnum] + '/NMME_Data_{0:.5f}_{1:.5f}.txt'.format(lat_blocks[ii],lon_blocks[jj])
    else:
        nmme_filename = input_directory + 'NMME/InitMonth_0' + str(init_month) +"/" + fullvariable_name[varnum] + '/NMME_Data_{0:.5f}_{1:.5f}.txt'.format(lat_blocks[ii],lon_blocks[jj])
    
    print(nmme_filename)
    nmme_data = pd.read_csv(nmme_filename, header=0, index_col=None, usecols = req_cols_nmme)    
    filepath = outpath +'/InitMonth_0' + str(init_month) + '/' + fullvariable_name[varnum]
    grid_filepath = outpath +'/InitMonth_0' + str(init_month) + '/' + fullvariable_name[varnum] + '/' + 'Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj])
    
    # Check if output directory exists or not. If yes, delete previously generated files to avoid overwritting/appending. If not, then create output directory
    try:
        if os.path.isdir(outpath + '/InitMonth_0' + str(init_month)) == False:
            os.mkdir(outpath + '/InitMonth_0' + str(init_month))
        if os.path.exists(filepath) == False:    
            os.mkdir(filepath)
        if os.path.exists(grid_filepath) == False:    
            os.mkdir(grid_filepath)              
        else:
            try:
                for f in os.listdir(grid_filepath):
                    os.remove(os.path.join(grid_filepath,f))
            except IOError:
                raise IOError("Cannot Delete Previously Generated Files")
    except:
        pass
        
    year = np.arange(1982,2011,1)                 # Year Range of the data

    for y in year:                                 # Running loop for each year
        # calculating number of days and years of data needs to append before NMME data, it depends upon initialization month
        inidate = dt.date(y-3,1,1)                  
        enddate = dt.date(y,12,31)
        delta = enddate - inidate
        days = delta.days + 1
        print(y)
        
        # Hardcoding the starting and ending DOY for appending gridMET data
        if init_month <= 4:
            if calendar.isleap(y)==True:
                front_gmet_str = {"1":31, "2":60, "3":91, "4":121, "5":152, "6":182, "7":213, "8":244, "9":274, "10":305, "11":335, "12":366}
                back_gmet_str = {"1":275, "2":306, "3":336, "4":1, "5":32, "6":61, "7":92, "8":122, "9":153, "10":183, "11":214, "12":245}
            else:
                front_gmet_str = {"1":31, "2":59, "3":90, "4":120, "5":151, "6":181, "7":212, "8":243, "9":273, "10":304, "11":334, "12":365}
                back_gmet_str = {"1":274, "2":305, "3":335, "4":1, "5":32,"6":60,"7":91,"8":121,"9":152,"10":182,"11":213,"12":244}
        else:
            if calendar.isleap(y+1)==True:
                front_gmet_str =  {"1":31, "2":60, "3":91, "4":121, "5":152, "6":182, "7":213, "8":244, "9":274, "10":305, "11":335, "12":366}
                back_gmet_str = {"1":275, "2":306, "3":336, "4":1, "5":32, "6":61, "7":92, "8":122, "9":153, "10":183, "11":214, "12":245}
            elif calendar.isleap(y)==True:
                front_gmet_str =  {"1":31, "2":60, "3":91, "4":121, "5":152, "6":182, "7":213, "8":244, "9":274, "10":305, "11":335, "12":366}
                back_gmet_str = {"1":274, "2":305, "3":335, "4":1, "5":32,"6":60,"7":91,"8":121,"9":152,"10":182,"11":213,"12":244}
            else:
                front_gmet_str = {"1":31, "2":59, "3":90, "4":120, "5":151, "6":181, "7":212, "8":243, "9":273, "10":304, "11":334, "12":365}
                back_gmet_str = {"1":274, "2":305, "3":335, "4":1, "5":32,"6":60,"7":91,"8":121,"9":152,"10":182,"11":213,"12":244}

        #===============================================================================================
        #    Based on initialization month, subsetting gridMET data that will be appended to NMME data
        #===============================================================================================         
        tic = time.perf_counter()
        if init_month <= 3:
            gmet_back = gmet_data[(gmet_data.Year == y)]
            gmet_front = gmet_data[(gmet_data.Year == y-3) | (gmet_data.Year == y-2) | (gmet_data.Year == y-1) | ((gmet_data.Year == y) & (gmet_data.DOY <= front_gmet_str[str(init_month)]))]
            gmet_back_yr = gmet_back[(gmet_back.DOY >= back_gmet_str[str(init_month)])]
            nmme_yr = nmme_data[(nmme_data.Initialization_Year == y)]
            gmet_back_yr = gmet_back_yr.dropna()
            gmet_front = gmet_front.dropna()
            if varnum == 4:
                nmme_yr = nmme_yr.fillna(0)
            else:
                nmme_yr = nmme_yr.dropna()
        elif init_month == 4:
            gmet_front = gmet_data[(gmet_data.Year == y-3) | (gmet_data.Year == y-2) | (gmet_data.Year == y-1) | ((gmet_data.Year == y) & (gmet_data.DOY <= front_gmet_str[str(init_month)]))]
            nmme_yr = nmme_data[(nmme_data.Initialization_Year == y)]
            gmet_front = gmet_front.dropna()
            if varnum == 4:
                nmme_yr = nmme_yr.fillna(0)
            else:
                nmme_yr = nmme_yr.dropna()
        elif init_month > 4:
            gmet_front = gmet_data[(gmet_data.Year == y-3) | (gmet_data.Year == y-2) | (gmet_data.Year == y-1) | ((gmet_data.Year == y) & (gmet_data.DOY <= front_gmet_str[str(init_month)]))]
            gmet_back = gmet_data[(gmet_data.Year == y+1)]
            gmet_back_yr = gmet_back[(gmet_back.DOY >= back_gmet_str[str(init_month)])]     
            nmme_yr = nmme_data[(nmme_data.Initialization_Year == y)]
            gmet_back_yr = gmet_back_yr.dropna()
            gmet_front = gmet_front.dropna()
            if varnum == 4:
                nmme_yr = nmme_yr.fillna(0)
            else:
                nmme_yr = nmme_yr.dropna()
            enddate = dt.date(y+1,12,31)
            delta = enddate - inidate
            days = delta.days+1
        else:
            pass

        toc = time.perf_counter()
        print(f"Extracting 4-5 Years data took {toc - tic:0.4f} seconds")   
        print(days)
                
      
        lat = list(set(nmme_yr['Lat']))
        lon = list(set(nmme_yr['Lon']))

        lat.sort()
        lon.sort()

        all_data = np.zeros([days,len(lat),len(lon)], dtype=float)
        all_data[:,:,:] = np.nan
        
        # Checking if NMME and GridMET has same coordinates or not.
        try:
           lat_nmme = list(set(nmme_yr['Lat']))
           lon_nmme = list(set(nmme_yr['Lon']))
           lat_nmme.sort()
           lon_nmme.sort()
           isittrue = np.allclose(lat,lat_nmme) & np.allclose(lon,lon_nmme)
           f"GridMET and NMME Latitude are same {isittrue}"
        except Exception as e:
            print("NMME and GridMET Grids are not same")

        #===============================================================================================
        #           Combining NMME and GridMET data
        #=============================================================================================== 
    
        tic = time.perf_counter()
        for i, j in itt.product(range(len(lat)), range(len(lon))):
            nmme_yr_grid = nmme_yr[np.isclose(nmme_yr.Lat,lat[i]) & np.isclose(nmme_yr.Lon,lon[j])]
            nmme_yr_grid = nmme_yr_grid.drop_duplicates(subset=['Forecast_Year','Forecast_Month','Forecast_Day'], keep='first')                          
            nmme_var = nmme_header[varnum+7]
            nmme_final = nmme_yr_grid[nmme_var]
            del nmme_yr_grid
            if len(nmme_final) == 0:
                pass
            else:
                if init_month == 4:
                    gmet_yr_grid = gmet_front[np.isclose(gmet_front.Lat,lat[i]) & np.isclose(gmet_front.Lon,lon[j])]
                    output_file = grid_filepath + '/' + 'Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj]) + '_Merged_VIC_Forcings_init_month_' + str(init_month) + '_' + fullvariable_name[varnum] + '_' + str(y-3) + '_' + str(y) + '.nc' 
                    gmet_var = gridmet_header[varnum+4]
                    gridMET_final = gmet_yr_grid[gmet_var]
                    if (len(gridMET_final)> 0):
                        all_data[:,i,j] = pd.concat([gridMET_final,nmme_final],axis=0)
                    else:
                        all_data[:,i,j] = np.nan
                        print("GridMET has NaN values",lat[i],lon[j],'\n')
                   
                    del gmet_yr_grid
                    del gridMET_final
                    del nmme_final
                elif init_month < 4:
                    gmet_yr_grid = gmet_front[np.isclose(gmet_front.Lat,lat[i]) & np.isclose(gmet_front.Lon,lon[j])]
                    gmet_back_yr_grid = gmet_back_yr[np.isclose(gmet_back_yr.Lat,lat[i]) & np.isclose(gmet_back_yr.Lon,lon[j])]
                    output_file = grid_filepath + '/'  + 'Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj]) + '_Merged_VIC_Forcings_init_month_' + str(init_month) + '_' + fullvariable_name[varnum] + '_' + str(y-3) + '_' + str(y) + '.nc' 
                    gmet_var = gridmet_header[varnum+4]
                    gridMET_final = gmet_yr_grid[gmet_var]
                    gridMET_back_Final = gmet_back_yr_grid[gmet_var]
                    if ((len(gridMET_final)> 0) & (len(gridMET_back_Final)> 0)):
                        all_data[:,i,j] = pd.concat([gridMET_final,nmme_final,gridMET_back_Final],axis=0)
                    else:
                        all_data[:,i,j] = np.nan
                        print("GridMET has NaN values",lat[i],lon[j],'\n')
                    
                    del gmet_yr_grid
                    del gridMET_final
                    del gridMET_back_Final
                    del gmet_back_yr_grid
                    del nmme_final
                elif init_month > 4: 
                    gmet_yr_grid = gmet_front[np.isclose(gmet_front.Lat,lat[i]) & np.isclose(gmet_front.Lon,lon[j])]
                    gmet_back_yr_grid = gmet_back_yr[np.isclose(gmet_back_yr.Lat,lat[i]) & np.isclose(gmet_back_yr.Lon,lon[j])]
                    output_file = grid_filepath + '/'  + 'Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj]) + '_Merged_VIC_Forcings_init_month_' + str(init_month) + '_' + fullvariable_name[varnum] + '_' + str(y-3) + '_' + str(y+1) + '.nc' 
                    gmet_var = gridmet_header[varnum+4]
                    gridMET_final = gmet_yr_grid[gmet_var]
                    gridMET_back_Final = gmet_back_yr_grid[gmet_var]
                    if ((len(gridMET_final)> 0) & (len(gridMET_back_Final)> 0)):
                        all_data[:,i,j] = pd.concat([gridMET_final,nmme_final,gridMET_back_Final],axis=0)
                    else:
                        all_data[:,i,j] = np.nan
                        print("GridMET has NaN values",lat[i],lon[j],'\n')
                    
                    del gmet_yr_grid
                    del gridMET_final
                    del gridMET_back_Final
                    del gmet_back_yr_grid
                    del nmme_final


        toc = time.perf_counter()
        print(f"Saving GridWise 4-5 Years data took {toc - tic:0.4f} seconds")   

        #===============================================================================================
        #           Generating NetCDF files: Merged NMME and gridMET data
        #=============================================================================================== 
        try:
            ncfile = nc.Dataset(output_file,"w")
            ncfile.Conventions = "CF-1.6"
            ncfile.title = "VIC-CropSyst Input Forcing"
            ncfile.source = 'Downscaled NMME and GridMET merged together'
            ncfile.history = "Author Ashish Kondal and script is created on " + dt.date.today().isoformat()
            ncfile.date_created = str(dt.datetime.now())
            
            ncfile.start_date = inidate.isoformat()
            ncfile.end_date = enddate.isoformat()
            
            ncfile.createDimension("longitude", len(lon))
            ncfile.createDimension("latitude", len(lat))
            ncfile.createDimension("time", days)
            
         
            latvar = ncfile.createVariable("latitude", float, ("latitude",))
            latvar.long_name = "Latitude"
            latvar.units = "degrees_north"
            latvar[:] = lat

            lonvar = ncfile.createVariable("longitude", float, ("longitude",))
            lonvar.long_name = "Longitude"
            lonvar.units = "degrees_east"
            lonvar[:] = lon

            timevar = ncfile.createVariable("time", int, ("time",))
            timevar.long_name = "Time"
            timevar.units = "days since " + inidate.isoformat()
            timevar.calendar = 'gregorian'
            timevar[:] = range(0, days)

            data_var = ncfile.createVariable(fullvariable_name[varnum], float, ("time","latitude","longitude"))
            data_var.long_name = fullvariable_name[varnum]
            data_var.missing_value = -9999.0
            data_var[:,:,:] = all_data
            ncfile.close()
            del all_data
        except:
            print("Failed to Generate NetCDF")
        
#===============================================================================================
#           Reading the gridMET data of a particular gridblock
#===============================================================================================     

gridMET_filename = input_directory + 'GridMET/' + fullvariable_name[varnum] + '/GridMET_Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj])
print(gridMET_filename)
gmet_data = pd.read_csv(gridMET_filename, header=0, index_col=None, usecols = req_cols_gmet)
gmet_data = gmet_data.replace('None',np.nan)

#gmet_data.info(verbose=False,memory_usage="deep")  #check memory_usage
outpath = "/weka/data/lab/adam/ashish.kondal/NMME_ENSMEAN/Preprocessed_Data/VIC_CropSyst_Forcings/ASCII2NETCDF/" + basins[basin_id]
outpath1 = "/weka/data/lab/adam/ashish.kondal/NMME_ENSMEAN/Preprocessed_Data/VIC_CropSyst_Forcings/"
outpath2 = "/weka/data/lab/adam/ashish.kondal/NMME_ENSMEAN/Preprocessed_Data/VIC_CropSyst_Forcings/ASCII2NETCDF/"

# Create output directories
try:
    if os.path.exists(outpath1) == False:
        os.mkdir(outpath1)
    if os.path.exists(outpath2) == False:
        os.mkdir(outpath2)
    if os.path.exists(outpath) == False:
        os.mkdir(outpath)
except:
    pass

#===============================================================================================
#           Running the functions
#=============================================================================================== 

asc2nc(outpath,init_month,basin_id,gmet_data,varnum,ii,jj)

    
    

  