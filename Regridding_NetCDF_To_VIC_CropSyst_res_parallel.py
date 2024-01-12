#=================================================================
# Regriding NMME NetCDF Files to VIC-CropSyst's Resolution
# Author: Ashish Kondal
# Date: September 09, 2022
# Description: This script regrid Merged NMME + GridMET files to 1/16th resolution.
#=================================================================


#================================================================
#   IMPORTING LIBRARIES/PACKAGES
#================================================================

from __future__ import print_function
import sys
import os, string
import netCDF4 as nc
import scipy
import numpy as np
import pandas as pd
import datetime as dt
import itertools as itt

#================================================================
#   DEFINE VARIABLES AND INPUT FILES
#================================================================

init_month = sys.argv[1];                                       # Retrieving NMME's initialization month from sbatch script
init_month = int(init_month);

varnum = sys.argv[2]                                            # Retrieving "variables" number (goes from 1-8) from sbatch script
varnum = int(varnum)-1

basin_id = sys.argv[3]                                          # Retrieving "basins" number (goes from 1-4; see "basins") from batch script
basin_id = int(basin_id)-1

gridnum = sys.argv[4]                                          # Retrieving gridblock number (goes from 1-160) from batch script
gridnum = int(gridnum)-1
#f'basin = {basins[gridnum]}' # Goes from 0-159 in python counting

input_directory = 'Path_to_inputfiles/';
variables = ['_pr','_rhsmax','_rhsmin','_srad','_tasmax','_tasmin','_vs','_sph'];
fullvariable_name = ['precipitation_amount','max_air_temperature','min_air_temperature','wind_speed','specific_humidity','surface_downwelling_shortwave_flux_in_air','max_relative_humidity','min_relative_humidity'];
basins = ['okanogon','wallawalla','yakima','pnw']

f'Basin = {basins[basin_id]}'
f'Initialization = {init_month}'
f'Variable = {fullvariable_name[varnum]}'
f'Grid Number = {gridnum}'


# Latitudinal and Longitudinal Extent of Study Area
lat_blocks = np.arange(39.5,50.5,1)
lon_blocks = np.arange(-125,-108,1)


# Selecting coordinating based on current "gridnum"
cord_ind = pd.DataFrame(itt.product(range(len(lat_blocks)-1), range(len(lon_blocks)-1)))
ii = cord_ind.iloc[gridnum,0]
jj = cord_ind.iloc[gridnum,1]

outpath = input_directory + "VIC_CropSyst_Forcings/REGRIDDED/" 
outpath1 = input_directory + "VIC_CropSyst_Forcings/REGRIDDED/" + basins[basin_id] +  '/InitMonth_0'+ str(init_month) + '/' + fullvariable_name[varnum] +"/"
grid_filepath = outpath1 + 'Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj]) +'/'

# Check if output directory exists or not. If yes, delete previously generated files to avoid overwritting/appending. If not, then create output directory
try:
    if os.path.exists(outpath) == False:
        os.mkdir(outpath)
    if os.path.isdir(outpath + '/' + basins[basin_id])==False:
        os.mkdir(outpath + '/' + basins[basin_id])
    if os.path.isdir(outpath + '/' + basins[basin_id] + '/InitMonth_0' + str(init_month)) == False:
        os.mkdir(outpath + '/' + basins[basin_id] + '/InitMonth_0' + str(init_month))
    if os.path.exists(outpath1) == False:
        os.mkdir(outpath1)
    if os.path.exists(grid_filepath) == False:    
        os.mkdir(grid_filepath)    
except:
    try:
        for f in os.listdir(grid_filepath):
            os.remove(os.path.join(grid_filepath,f))
    except:
        pass
    
# Listing all files in input directory
path = input_directory + "VIC_CropSyst_Forcings/ASCII2NETCDF/" + basins[basin_id] + '/InitMonth_0'+ str(init_month) + '/' + fullvariable_name[varnum] +'/' + 'Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj]) +'/'
filelist = os.listdir(path)

# Retrieving Lat-Long from a input file to create mesh grid
test_filename = os.path.join(path,filelist[0])
test_data = nc.Dataset(test_filename)
lat = np.array(test_data['latitude'][:])
lon = np.array(test_data['longitude'][:])
time = np.array(test_data['time'][:])
test_data.close()
lat.sort()
lon.sort()
X, Y = np.meshgrid(lon,lat)

#========================================================================================================================
#                           READING SOIL PARAMETER FILE 
#   Reading Soil Parameter file (for VIC-Cropsyst runs) - coordinates in this file must match with our regridded coordinates. 
#   Thus, we extracted coordinates from this file and regridded our data for those coordinate to avoid any mismatch
#======================================================================================================================== 
soilfile_path = "soil_parameter_file.txt"
soil_data = pd.read_table(soilfile_path,header=None,index_col = None,sep=" ")
lat_vic1 = pd.unique(soil_data.iloc[:,2])               # Extracing lat's from soil file
lon_vic1 = pd.unique(soil_data.iloc[:,3])               # Extracing lon's from soil file
lat_vic1.sort()
lon_vic1.sort()

# Finding lat and lon in soil files which are closest or exact to the Lat and Lon in input files
lat_ind = np.zeros(shape=(len(lat),1))
for x in range(len(lat)):
   lat_ind[x] = (np.abs(lat_vic1-lat[x])).argmin()
   
lon_ind = np.zeros(shape=(len(lon),1))
for y in range(len(lon)):
   lon_ind[y] = (np.abs(lon_vic1-lon[y])).argmin()
   
lat_ind = lat_ind.tolist()
lon_ind = lon_ind.tolist()   
lat_ind = np.int64(lat_ind)
lon_ind = np.int64(lon_ind)

lat_vic = pd.unique(lat_vic1[lat_ind][:,0])
lon_vic = pd.unique(lon_vic1[lon_ind][:,0])

X1, Y1 = np.meshgrid(lon_vic,lat_vic)           # Creating meshgrid 

#==================================================================================================================
#                                   REGRIDDING INPUT FILES & SAVING IT INTO NETCDF FORMAT
# Coastal grid points are regridded with nearest neigbhour to avoid any nan's or zero in the final data.
# Inland grid points are regridded with linear interpolation 
#==================================================================================================================

for f in filelist:
    filename = os.path.join(path,f)
    data = nc.Dataset(filename)
    time = np.array(data['time'][:])
    nmme_var = fullvariable_name[varnum]
    all_data = np.zeros([len(time),len(lat_vic),len(lon_vic)], dtype=float)
    all_data[:,:,:] = -9999
    outfile = grid_filepath + 'Regridded_' + f
    for i in range(len(time)):
        vardata = np.array(data[nmme_var][i,:,:])
        if (sum(np.isnan(vardata).flatten()) > 0):
            all_data[i,:,:] = scipy.interpolate.griddata((X.flatten(),Y.flatten()),vardata.flatten(), (X1,Y1),method='nearest')
        else:
            all_data[i,:,:] = scipy.interpolate.griddata((X.flatten(),Y.flatten()),vardata.flatten(), (X1,Y1),method='linear')
        
        del vardata

    data.close()

    try:
        ncfile = nc.Dataset(outfile,"w")
        ncfile.Conventions = "CF-1.6"
        ncfile.title = "VIC-CropSyst Input Forcing - Regridded"
        ncfile.source = 'Downscaled NMME and GridMET merged together and regridded to 0.0625'
        ncfile.history = "Author Ashish Kondal and script is created on " + dt.date.today().isoformat()
        ncfile.date_created = str(dt.datetime.now())
        
        
        ncfile.createDimension("longitude", len(lon_vic))
        ncfile.createDimension("latitude", len(lat_vic))
        ncfile.createDimension("time", len(time))
        
     
        latvar = ncfile.createVariable("latitude", float, ("latitude",))
        latvar.long_name = "Latitude"
        latvar.units = "degrees_north"
        latvar[:] = lat_vic

        lonvar = ncfile.createVariable("longitude", float, ("longitude",))
        lonvar.long_name = "Longitude"
        lonvar.units = "degrees_east"
        lonvar[:] = lon_vic

        timevar = ncfile.createVariable("time", int, ("time",))
        timevar.long_name = "Time"
        timevar.units = "days since file title years"
        timevar.calendar = 'gregorian'
        timevar[:] = time

        data_var = ncfile.createVariable(fullvariable_name[varnum], float, ("time","latitude","longitude"))
        data_var.long_name = fullvariable_name[varnum]
        data_var.missing_value = -9999.0
        data_var[:,:,:] = all_data
        ncfile.close()
    except:
        print("Failed to Generate NetCDF")
        