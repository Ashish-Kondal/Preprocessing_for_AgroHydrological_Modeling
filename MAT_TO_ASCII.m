%===========================================================
% Title: Converting .mat format files into ASCII files.
% Author: Ashish Kondal (ashish.kondal@wsu.edu)
% Date Created: August 24, 2022
% Description: Create ASCII files for different initializations from mat file.
%===========================================================
function [] = MAT_TO_ASCII(basin_id)  

%===========================================================
%	DEFINE VARIABLES & INPUT DIRECTORY/SETTINGS
%===========================================================

	%disp("1 = Okanogon, 2 = wallawalla, 3 = Yakima, 4 = pnw") % Number denotes the basin_id
	basins = {'okanogon','wallawalla','yakima','pnw'};
	
	input_directory = 'inputpath/NMME_ENSMEAN/Data/';
	filename_prefix = 'nmme_bcsd_ENSMEAN_daily_';

	output_directory = 'outputpath/NMME_ENSMEAN/Preprocessed_Data/SubsetFiles/';

%===========================================================
%	READING VARIABLES STORED IN .MAT FILES 
%===========================================================

% 1. Reading Metadata File

	metadata_file = strjoin([input_directory,filename_prefix,basins(basin_id),'_metadata.mat'],"");
	metadata = load(metadata_file);
			
% 2. Convert datenum into actual dates
	if (basin_id == 1 || basin_id == 4)
		pseudo_metadata_file = strjoin([input_directory,filename_prefix,basins(3),'_metadata.mat'],"");
		pseudo_metadata = load(pseudo_metadata_file);
		[Forecast_Date] = datevec(pseudo_metadata.time_forecast + datenum(1900,1,1)); % Converting vector (datenum) to dates.
		[Initialization_Date] = datevec(pseudo_metadata.time_initialization + datenum(1900,1,1)); % Converting vector (datenum) to dates.
	else
		[Forecast_Date] = datevec(metadata.time_forecast + datenum(1900,1,1)); % Converting vector (datenum) to dates.
		[Initialization_Date] = datevec(metadata.time_initialization + datenum(1900,1,1)); % Converting vector (datenum) to dates.
	end	

% 3. Get Latitude and Longitude values
	lat = metadata.lat;
	lon = metadata.lon - 360;

% 3. Load Variable Files
	if (basin_id ==4)
		precip_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_pr_1982_2020.mat'],"");
		rhsmax_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_rhsmax_1982_2020.mat'],"");
		rhsmin_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_rhsmin_1982_2020.mat'],"");
		sph_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_sph_1982_2020.mat'],"");
		srad_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_srad_1982_2020.mat'],"");
		tasmax_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_tasmax_1982_2020.mat'],"");
		tasmin_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_tasmin_1982_2020.mat'],"");
		vs_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_vs_1982_2020.mat'],"");
	else 
		precip_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_pr.mat'],"");
		rhsmax_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_rhsmax.mat'],"");
		rhsmin_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_rhsmin.mat'],"");
		sph_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_sph.mat'],"");
		srad_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_srad.mat'],"");
		tasmax_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_tasmax.mat'],"");
		tasmin_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_tasmin.mat'],"");
		vs_file  = strjoin([input_directory,filename_prefix,basins(basin_id),'_vs.mat'],"");
	end
	precip = load(precip_file);
	rhsmax = load(rhsmax_file);
	rhsmin = load(rhsmin_file);
	sph = load(sph_file);
	srad = load(srad_file);
	tasmax = load(tasmax_file);
	tasmin = load(tasmin_file);
	vs = load(vs_file);
	varnames = {'Lat','Lon','Initialization_Year','Initialization_Month','Forecast_Year','Forecast_Month','Forecast_Day','Precipitation_mm', 'Max_Temperature_K', 'Min_Temperature_K', 'Wind_Speed_mps','Specific_Humidity_Kg_per_Kg','ShortWave_Radiation_W_per_m2','Max_Relative_Humidity','Min_Relative_Humidity'};

%===========================================================
%	WRITING/STORING VARIABLES IN CSV FILES 
%===========================================================
	
	for i = 1:length(lat)
		clat = lat(i);
		for j = 1:length(lon)
			data = NaN([length(Initialization_Date),15],'double');
			clon = lon(j);
			disp("Current Grid: %0.5f %0.5f \n", clat, clon)
			data(:,1:2) = repmat([clat clon],length(Forecast_Date),1);
			data(:,3:4) = Initialization_Date(:,1:2);
			data(:,5:7) = Forecast_Date(:,1:3);
			temp_data = precip.scale_factor * precip.daily_pr(i,j,:,:,:) + precip.add_offset;
			data(:,8) = temp_data(:);
			temp_data = double(NaN);
			temp_data = double(tasmax.scale_factor) * double(tasmax.daily_tasmax(i,j,:,:,:)) + tasmax.add_offset;
			data(:,9) = temp_data(:);
			temp_data = double(NaN);
			temp_data = double(tasmin.scale_factor) * double(tasmin.daily_tasmin(i,j,:,:,:)) + tasmin.add_offset;
			data(:,10) = temp_data(:);
			temp_data = NaN;
			temp_data = vs.scale_factor * vs.daily_vs(i,j,:,:,:) + vs.add_offset;
			data(:,11) = temp_data(:);
			temp_data = NaN;
			temp_data = sph.scale_factor * sph.daily_sph(i,j,:,:,:) + sph.add_offset;
			data(:,12) = temp_data(:);
			temp_data = NaN;
			temp_data = 0.1 * srad.daily_srad(i,j,:,:,:);
			data(:,13) = temp_data(:);
			temp_data = NaN;
			temp_data = rhsmax.scale_factor * rhsmax.daily_rhsmax(i,j,:,:,:) + rhsmax.add_offset;
			data(:,14) = temp_data(:);
			temp_data = NaN;
			temp_data = rhsmin.scale_factor * rhsmin.daily_rhsmin(i,j,:,:,:) + rhsmin.add_offset;
			data(:,15) = temp_data(:);
			temp_data = NaN;
			NaNrow = find(isnan(data(:,3))==false);
			data1 = data(NaNrow,:);
			clear data;
			for num_month = 1:1:12
				month_index = find(data1(:,4)== num_month);
				data2 = data1(month_index,:);
				output_filename = strjoin([output_directory,"BCSD_NMME_data_",basins(basin_id),"_ExtractedForInitializationMonth_",string(num_month),".txt"],"");
				clear month_index;
				if (i == 1 & j == 1)
					data2 = array2table(data2,'VariableNames',varnames);
					writetable(data2,output_filename,'WriteVariableNames',true);
					clear data2;
				else
					data2 = array2table(data2,'VariableNames',varnames);
					writetable(data2,output_filename,'WriteMode','append','WriteVariableNames',false);
					clear data2;
				end	
			end
			clear data1;
		end
	end
	clear all;	
end