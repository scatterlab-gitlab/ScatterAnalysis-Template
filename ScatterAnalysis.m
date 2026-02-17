%% README
%---------------------------------------------------------------------------------------------------
% Scatter Analysis script with code to run all types of experiments
%---------------------------------------------------------------------------------------------------

% NOTE: Conversion from Image to Real Life Measurement is 239.3 pix/mm for AAS

%% IMPORTANT PLOTTING NOTE!!!
%---------------------------------------------------------------------------------------------------
% Format all polots like this:

% f1 = figure('visible','off');
% ax1 = gca;
% 
% plot(x,y)
% 
% f1.Property = value
% ...
% 
% ax1.Property = value
% ...
%---------------------------------------------------------------------------------------------------
% Removing whitespace. This should be performed after all axis titles are added:

% % Get the OuterPosition and TightInset
% outerpos = ax1.OuterPosition;
% ti = ax1.TightInset;
% 
% % Calculate the new position to remove the whitespace
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% 
% % Set the new position of the axes
% ax1.Position = [left bottom ax_width ax_height];
%---------------------------------------------------------------------------------------------------

% Clear all variables, images, and command window to prepare for a new analysis 
clear; close all; clc;

%% FILE AND EXPERIMENT TYPE

% Mock data files for debugging
mock_AAS_data_path = '/Users/scatterlab/CSU Fullerton Dropbox/Scatter Lab/Shared/MiscFiles/Mock-AAS-Data/PL-13649-Mock-AAS-Data';
mock_ARS_data_path = '/Users/scatterlab/CSU Fullerton Dropbox/Scatter Lab/Shared/MiscFiles/Mock-ARS-Data/PL-13649-Mock-ARS-Data';

% String path to the folder will all the data
folder.data_path = '/Users/scatterlab/CSU Fullerton Dropbox/Scatter Lab/Shared/GWPAC_Lab_Data/ARS_TRS/2024_7_12_11_20_27_Crystaline_Silicon_SLED_ARS';

% AAS, ARS, CRYO, or TRS (Case Sensitive!!!)
experiment = 'ARS';

% Note to save as a text file about analysis specifics or changes
note_text = 'Trying to get good axies to the CCD images';

% Sample Name - NO UNDERSCORES!!! Use spaces
sample = 'Getting Axies';

%% DEBUGGING/PROCESSING VARIABLES
%---------------------------------------------------------------------------------------------------
% Logic variables to have the script only run certain parts to adjust parameters for the analysis
%---------------------------------------------------------------------------------------------------

% DEBUGGING
%---------------------------------------------------------------------------------------------------

% If DEBUG = 1, enters debug mode which means we are running the script to make adjustment
% If DEBUG = 0, we are running script to post a report on measurements
DEBUG = 1;

% If locateROI = 1, script will only display 'image_Selector' and let you adjust ROI
% If locateROI = 0, script runs normally
locate_ROI = 1;

% Number of images that will be analysed (DOES NOT effect ARS)
total_images = 3;

% Assign a value for n which will be the image number you want to evaluate
% before the entire script begins. This is good for checking the distance
% forward for the end images (and possibly claim values for farther images as
% well) Only works when locateROI == 1
image_selector = 3;

% PROCESSING
%---------------------------------------------------------------------------------------------------

% If printImages = 1, all images will be saved.
% If printImages = 0, no images will be saved; analysis graphs will be saved.
print_images = 1;

% If darkImageSave = 1 dark images will be processed and saved along with subtracted images
% If darkImageSave = 0 only subtracted images will be saved
save_dark_images = 0;

% If powermonitor = 1 the code will use Thorlabs power monitor data to calculate BRDF (Pickoff Power)
% If powermonitor = 0 it will use only incidentPower (initial measurement)
use_power_monitor = 1;

% If singledark = 1 only uses one dark image for all images (For TRS)
% If singledark = 0 uses all dark images to process subtracted image
single_dark = 0;

% If backgroundimage_texp = 1 it will use background based on exposure time
% If backgroundimage_texp = 0 it will chose background based on angle
use_background_exp = 0;

% If line removal = 1, it will replace CCD pixel line 1569 with 1570 for all images due to a hot line issue
% with the ARS camera 
% If lineremoval = 0, no additional post-processing will be done to ARS images
lineremoval = 1;

%% USER DEFINED IMAGE/ROI VARIABLES
%---------------------------------------------------------------------------------------------------
% Edit all these variables as necessary for the specific experiment being conducted
%---------------------------------------------------------------------------------------------------

% Number of ROIs wanted to be analized for background subtraction (typically 6)
numROIs = 6;

% Color Limit values min = 0, max = low enough to see everything
clim_Min = 0;
clim_Max = 1000;

% Color limit values for dark images only
dark_clim_Min = 0;
dark_clim_Max = 3000;

% Min x value of ROI
ROIxmin = 1602;
% Max x value of ROI
ROIxmax = 2721;

% Min y value of ROI
ROIymin = 1044;
% Max y value of ROI
ROIymax = 2050;

% How far forward is the observed surface of the optic from the center of the
% ARS optic holder on the table. [pixels]
df = 900;

% Dark region scaling factor mask parameters
% Find a spot on the optic near the ROI (but not inside of it) that will find black body radiation
Region_xcenter = 3453;
Region_ycenter = 431;
Region_xwidth = 200;
Region_ywidth = 200;

% OPTIONAL: Shift the region you found if its a good shape, but you don't want to find/retype all four values
% X: Shift left = NEGATIVE | Shift right = POSITIVE
% Y: Shift up = NEGATIVE | Shift down = POSITIVE
ROIxShift = 0;
ROIyShift = 0;

% Factor to scale up each outer ROI circle
switch experiment
  case 'AAS'
    ROIscaleFactor = 0.07;
  case 'TRS'
    ROIscaleFactor = 0.20;
  case 'CRYO'
    ROIscaleFactor = 0.4;
end

%% USER DEFINED EXPERIMENT VARIABLES
%---------------------------------------------------------------------------------------------------
% Edit all these variables as necessary for the specific experiment being conducted
%---------------------------------------------------------------------------------------------------

% Incident power [W] average before and after power measurement
incident_power = 0.0079;

% If using Neutral Density Filter (like for Spectralon) filter = 1/273, otherwise filter = 1)
filter = 1;

% Define legs for reflected light
switch experiment
  case 'ARS'
    % Enter the two legs of the right triangle made by the incident and reflected beams
    incident_leg = 222;
    separation_leg = 12;

    % Correction angle is calculated using arc tan() in degrees
    correction_angle = atand(separation_leg/incident_leg);
  case 'AAS'
    % Enter the two legs of the right triangle made by the incident and reflected beams
    incident_leg = 91.5;
    separation_leg = 10;

    % Correction angle is calculated using arc tan() in degrees
    correction_angle = atand(separation_leg/incident_leg);
  case 'CRYO'
    % Enter the two legs of the right triangle made by the incident and reflected beams
    incident_leg = 61;
    separation_leg = 8;

    % The correction angle for each experiment defined as half of the full angle between the incident laser light and the reflected laser light
    correction_angle = atand(separation_leg/incident_leg) / 2;
end

%% GRAPHICS DEFAULTS
%---------------------------------------------------------------------------------------------------
% Color Definitions for all grphs. This keeps formatting consistent throughout
%---------------------------------------------------------------------------------------------------

% Blue, Red, Orange, Green, Purple
colors = {[0 0 1] [0.81 0 0] [1 0.48 0] [0 0.75 0.2] [.87 .32 1]};

set(groot,'DefaultAxesFontName','Times New Roman');
set(groot,'DefaultAxesFontSize',13);
set(groot,'DefaultTextFontName', 'Times New Roman');
set(groot,'DefaultTextFontSize',13);

set(groot,'DefaultAxesColorOrder', ...
    [0.10 0.76 0.82;...
    0.94 0.66 0.07;...
    0.81 0 0;...
    0 0.75 0.16; ...
    .87 .32 1]);

%% EXPERIMENT SPECIFIC VARIABLES
%---------------------------------------------------------------------------------------------------
% The rest of these variables don't get changed that offten, but can be changed if necessary
%---------------------------------------------------------------------------------------------------
switch experiment
  case {'AAS','CRYO','TRS'}
    % Angular orrientation of optic
    theta_rot_stage = 0;
end

% Azimuthal angle (zero in the plane of the laser beam)
phi_s = 0;

% Calibration Factor found by performing a calibration based on the particular experiment
switch experiment
  case {'ARS','TRS'}
    muFcal = (4.86292e-13);
  case 'AAS'
    muFcal = 1.3867e-12;
  case 'CRYO'
    muFcal = (3.847e-13);
end

% Pickoff Power Correction Coefficients Ordered from highest order to lowest order coefficient
switch experiment
  case 'ARS'
    correction_coefficients = [4.275448819818776e+04,2.945571263780217e+04,-4.111329671854571e+03,...
                               3.765240633203828e+02,65.570193682492090,-0.003145133088967];
  case 'AAS'
    correction_coefficients = 9.2264;
  case 'CRYO'
    correction_coefficients = [19.5589,-0.0001];
  case 'TRS'
    % The calculated calibration constant for converting the TRS monitor power into
    % the incident power, used in line 437
    Pmon_const_trs = 101.5;
end

% Wavelength in [nm] - AAS(SLED-1041), ARS(SLED-1056, Mephisto-1064), CRYO(SLED-1900)
switch experiment
  case 'AAS'
    lambda = 1041e-9;
  case 'ARS'
    lambda = 1056e-9;
  case 'CRYO'
    lambda = 1999.9e-9;
end

 % Conversion factor used to convert pixels of CCD image to real life measurement in mm
 switch experiment
   case 'AAS'
     pix_per_mm = 239.3;
   case 'ARS'
     pix_per_mm = 239.3; % THIS IS VERY WRONG!!!! ONLY USING THIS FOR DEBUGGING PURPOSES!!!!!!!!!!
 end

% Uncertainty in power monitor measurement from Thorlab specs
uncert_pm = .05;

% Uncertainty in measuring the solid angle (see calibration script)
uncert_omega = .033;

% Uncertainty in laser power stability
uncert_laspower = .004;

% Calculate systematic uncertainty budget
sys_uncert = uncert_pm + uncert_omega + uncert_laspower;

%% FOLDER LOGIC CHECKS
%---------------------------------------------------------------------------------------------------
% Check if the slash is forwad for Mac or backward for Windows
%---------------------------------------------------------------------------------------------------

folder.slashCheck = contains(pwd, '/');

if folder.slashCheck == 1
  slash = '/';
elseif folder.slashCheck == 0
  slash = '\';
  % Check if the path has a C: or F: in front of it to determine Path Type
  folder.cCheck = contains(pwd,'C:');
  if folder.cCheck == 1
    folder.windows = 'C:';
  elseif folder.cCheck == 0
    folder.windows = 'F:';
  end
end

switch experiment
  case 'CRYO'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine what type of CRYO run this is (ColdRun or ColdTrap)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Look for 'Cold_Run'
    % If it returns 1, Cold Run
    % If it returns 0, Cold Trap
    folder.runLogic = contains(folder.data_path,'Cold_Run');
    
    switch folder.runLogic
      % Cold Run
      case 1
        folder.runType = 1;
        % Cold Trap
      case 0
        folder.runType = 0;
    end
end

% Add a slash at the end of the folder path for easy workability within code
folder.data_path = [folder.data_path,slash];

% When printing sample to a graph, underscores will subscript the next character, which is why we
% use spaces in its name. For folder paths, we replace spaces (' ') with underscores ('_')
folder.sample_name = strrep(sample,' ','_');

%% ROI VARIABLE ADJUSTMENTS
%---------------------------------------------------------------------------------------------------
% This section makes the ROI a perfect circle
%---------------------------------------------------------------------------------------------------

% Make ROI a circle
ROICircle = (max([ROIxmax-ROIxmin,ROIymax-ROIymin]) - min([ROIxmax-ROIxmin,ROIymax-ROIymin])) / 2;
if ROIxmax-ROIxmin > ROIymax-ROIymin
  ROIymax = ROIymax + ROICircle;
  ROIymin = ROIymin - ROICircle;
else
  ROIxmax = ROIxmax + ROICircle;
  ROIxmin = ROIxmin - ROICircle;
end


% ROI Shifting Calculations
ROIxmin = ROIxmin + ROIxShift;
ROIxmax = ROIxmax + ROIxShift;
ROIymin = ROIymin + ROIyShift;
ROIymax = ROIymax + ROIyShift;

% Width of the ROI
ROIwidth = ROIxmax - ROIxmin;
% Height of the ROI
ROIheight = ROIymax - ROIymin;

% The center x value for the ROI
ROIxcenter = (ROIwidth/2)+ ROIxmin;
% The center y value for the ROI
ROIycenter = (ROIheight/2)+ ROIymin;


%% PATHS AND FOLDERS
%---------------------------------------------------------------------------------------------------
% Create path variables to all data files from each experiment
%---------------------------------------------------------------------------------------------------

% Paths to Images/Dark Images
imagedir = [folder.data_path, 'images', slash];
darkimagedir = [folder.data_path, 'darkimages', slash];

% String that says when analysis was run
folder.timeAnalysis = ['Analyzed_',char(datetime('now','Format','yyyy_MM_dd_HH_mm'))];

% String that determines where analysis results will be stored
folder.analysis_folder = strrep(folder.data_path,'Data','Analysis');
folder.analysisPath = [folder.analysis_folder,folder.timeAnalysis,slash];

% Creates a separate folder if debug mode is enabled
if DEBUG == 1    
  folder.folder2replace = ['GWPAC_Lab_Analysis',slash];
  folder.debug_path = ['GWPAC_Lab_Analysis',slash,'Debug_'];
  folder.analysisPath = strrep(folder.analysisPath,folder.folder2replace,...
                                folder.debug_path);
  % Add the distance forward to the name for ARS analysis
  switch experiment
    case 'ARS'
      % If doing ARS debug, automatically adds df to end of analysis path.
      folder.distanceforward_path = [folder.timeAnalysis,'_df',sprintf('%.0f',df)];   
      folder.analysisPath = strrep(folder.analysisPath,folder.timeAnalysis,folder.distanceforward_path);
  end    
end

% creates structure that includes sample, time, and analysis path 
% directories (makes easier to use in functions)
directories = struct('spath',folder.sample_name,'tpath',folder.timeAnalysis,'apath',folder.analysisPath);

if locate_ROI == 1
  % Don't create any folders if locating ROI

elseif locate_ROI == 0
  %---------------------------------------------------------------------------------------------------
  % Creates all the analysis output directories.
  %---------------------------------------------------------------------------------------------------
  
  % Make a folder for the analysis results
  mkdir(folder.analysisPath);
  
  % Make a folder for the analysis PNG files
  mkdir([folder.analysisPath,'PNG']);
  
  % Make a folder for the analysis FIG files
  mkdir([folder.analysisPath,'FIG']);
  
  % Make a folder for the output plots
  mkdir([folder.analysisPath, folder.sample_name]);
  
  % Make a folder for the analysis video
  mkdir([folder.analysisPath,'Movie']);
  
  % Make a folder for camera image pngs for video later
  mkdir([folder.analysisPath,'CCD_PNG']);

  if save_dark_images == 1
    % Make a folder for dark camera image pngs
    mkdir([folder.analysisPath,'Dark_CCD_PNG']);
  end
end

%% EXPERIMENT DATA GATHERING
%---------------------------------------------------------------------------------------------------
% Import all of the measurement data
%---------------------------------------------------------------------------------------------------

% Read all data files from their paths
switch experiment
  case 'ARS'
    experiment_data = readtable([folder.data_path,'PowermeterData.txt'],'VariableNamesLine',1,'ReadVariableNames',true,'Delimiter','comma');
  case 'AAS'
    experiment_data = readtable([folder.data_path,'OvenData.txt'],'VariableNamesLine',1,'ReadVariableNames',true,'Delimiter','comma');
  case 'CRYO'
    experiment_data = readtable([folder.data_path,'CryoData.txt'],'VariableNamesLine',1,'ReadVariableNames',true,'Delimiter','comma');
    cam_data = readtable([folder.data_path, 'CamData.txt'],'VariableNamesLine',1,'ReadVariableNames',true,'Delimiter','comma');
    cam_setting = readtable([folder.data_path, 'CamInfo.txt'],'VariableNamesLine',1,'ReadVariableNames',true,'Delimiter','comma');
end

% Extract each variable from the data file
switch experiment
  case 'ARS'
    % DateTime, Theta_rs_Degrees, Pmon_W, Pic_ID, Temperature_C
    total_images = height(experiment_data);
    date_time = datetime(experiment_data.DateTime,'Format','yyyy_MM_dd_HH_mm_ss', 'TimeZone', 'America/Los_Angeles');
    theta_rot_stage = experiment_data.Theta_rs_Degrees;
    pickoff_power = experiment_data.Pmon_W;
    pic_ID = experiment_data.Pic_ID;
    % exposure_time

  case 'AAS'
    % DateTime,Oven_Temp_C,Oven_Power_percent,Working_Setpoint_C,Pmon1_W,Pmon2_W,Table_Temp_C,Door_Temp_C,Cam_Temp_C,Humidity
    date_time = datetime(experiment_data.DateTime,'Format','yyyy-MM-dd_HH_mm_ss', 'TimeZone', 'America/Los_Angeles');
    temperature = experiment_data.Oven_Temp_C;
    heater_power = experiment_data.Oven_Power_percent;
    temp_setpoint = experiment_data.Working_Setpoint_C;
    pickoff_power = experiment_data.Pmon1_W;
    transmitted_power = experiment_data.Pmon2_W;
    table_temp = experiment_data.Table_Temp_C;
    door_temp = experiment_data.Door_Temp_C;
    outside_cam_temp = experiment_data.Cam_Temp_C;
    humidity = experiment_data.Humidity;
    % image time
    % exposure_time

  case 'CRYO'
    % DateTime,ElapsedTime,LakeshoreSetpoint_K,ArraySetpoint_K,TemperatureA_K,TemperatureB_K,PressureA_mbar,PressureB_mbar,Heater_Power,Pickoff_W,Transmitted_W
    date_time = datetime(experiment_data.DateTime,'Format','yyyy-MM-dd_HH_mm_ss');
    elapsed_time = experiment_data.ElapsedTime;
    lakeshore_setpoint = experiment_data.LakeshoreSetpoint_K;
    array_setpoiont = experiment_data.ArraySetpoint_K;
    chamber_temp = experiment_data.TemperatureA_K;
    stinger_temp = experiment_data.TemperatureB_K;
    chamber_pressure = experiment_data.PressureA_mbar;
    stinger_pressure = experiment_data.PressureB_mbar;
    heater_percent = experiment_data.Heater_Power;
    pickoff_power = experiment_data.Pickoff_W;
    transmitted_power = experiment_data.Transmitted_W;
    image_time = datetime(cam_data.DateTime,'InputFormat','yyyy_MM_dd_HH_mm_ss');
    exposure_time = cam_setting.ExposureTime;
    image_duration_time = image_time - image_time(1);
end



%% (BETA) IMAGE VARIABLES
%---------------------------------------------------------------------------------------------------
% Interpolates the variabels that took place during image capture
%---------------------------------------------------------------------------------------------------

% image_temp(n) = interp1(date_time, chamber_temp, image_time(n),'nearest','extrap');
% image_power(n) = interp1(date_time, power_corrected, image_time(n), 'nearest','extrap');
% image_transmitted = interp1(date_time, transmitted_power, image_time(n),'nearest','extrap');


%% PREALLOCATIONS
%---------------------------------------------------------------------------------------------------
% Preallocations of variabels for better script performance
%---------------------------------------------------------------------------------------------------

% Determine amount of images being analyzed in the video
switch experiment
  case {'AAS','CRYO','TRS'}
      pic_ID = 1:1:total_images;
end

% Correct the pickoff power to the actual incident power
% Use polyval(p,x) to evaluate the polynomial p at each point in x.
% The argument p is a vector of length n+1 whose elements are the coefficients
% (in descending powers) of an nth-degree polynomial
if use_power_monitor == 1
  power_corrected = polyval(correction_coefficients,pickoff_power);
else
  power_corrected = incident_power .* ones(total_images,1);
end

switch experiment
  case {'CRYO'}
    image_temp = interp1(date_time, chamber_temp, image_time,'nearest','extrap');
    image_power = interp1(date_time, power_corrected, image_time, 'nearest','extrap');
    image_transmitted = interp1(date_time, transmitted_power, image_time,'nearest','extrap');
end

% We want dates to be made into a durations (hh:mm:ss) for easy reading of time of data
duration_time = date_time - date_time(1);

% theta_s, the scattering angle, is the polar angle from the normal of the 
% sample to the camera
theta_s = theta_rot_stage + correction_angle;

% Make incidentPower, filter, and phi_s vectors by multiplying by ones matrix
filter = filter .* ones(total_images,1);
phi_s = phi_s .* ones(total_images,1);

% States the steps of the theta used to calculate ellipse
ellipse_theta = 0:0.01:2*pi;

% n = numImages
% z = numROIs
% i = sllipse_theta
switch experiment
  case {'AAS','ARS','TRS'}
    exposure_time = zeros(1,total_images);
end

texpB = zeros(1,total_images);
switch experiment
  case {'AAS','ARS','TRS'}
    image_time = NaT(1,total_images,'TimeZone','America/Los_Angeles');
    image_temp = zeros(1,total_images);
    image_power = zeros(1,total_images);
    image_duration_time = duration(zeros(total_images,3));
end


scaledWidth = zeros(1,numROIs);
scaledHeight = zeros(1,numROIs);
widthROI = zeros(1,numROIs);
heightROI = zeros(1,numROIs);
ybottom = zeros(1,numROIs);
ytop = zeros(1,numROIs);
xleft = zeros(1,numROIs);
xright = zeros(1,numROIs);
xradius_ellipse = zeros(1,numROIs);
xcenter_ellipse = zeros(1,numROIs);
yradius_ellipse = zeros(1,numROIs);
ycenter_ellipse = zeros(1,numROIs);
x_ellipse = zeros(numROIs,length(ellipse_theta));
y_ellipse = zeros(numROIs,length(ellipse_theta));
area = zeros(total_images,numROIs);
mask = struct;
mask_crop = struct;
RoI = struct;
RoI_mask = struct;
Bccd = zeros(total_images,numROIs);
BRDFccd = zeros(total_images,numROIs);
ARBccd = zeros(numROIs,total_images);
time = zeros(1,total_images);

b = zeros(1,total_images);
SE = zeros(1,total_images);
ci = zeros(1,total_images);
tot_uncert = zeros(1,total_images);
uncert_BRDFyint = zeros(1,total_images);

BRDFyint = zeros(1,total_images);

% determine BRDF limit based on darkest small box in image
num_of_div=8;

% number of divisions per line
LIMwidth=100;
LIMheight=100;

LIMy = zeros(1, num_of_div);
LIMx = zeros(1, num_of_div);
LIMx1 = zeros(1, num_of_div);
LIMx2 = zeros(1, num_of_div);
LIMy1 = zeros(1, num_of_div);
LIMy2 = zeros(1, num_of_div);
LIMpic = struct;
LIM_ARB = zeros(1, num_of_div);
LIM_B = zeros(1, num_of_div);
LIM_BRDF = zeros(1, num_of_div);

PSD_fx_fy = zeros(total_images,1);
S_iso = zeros(total_images,1);
f_quad_sum = zeros(total_images,1);
minLIM = zeros(1,total_images);

Omega_int_ring = zeros(length(theta_s));
Omega_ext_ring = zeros(length(theta_s));
Omega_ring = zeros(length(theta_s));
meanBRDFyint = zeros(length(theta_s));
intscat_a = zeros(length(theta_s));
theta_diff = zeros(length(theta_s));


s = 0;

%% EXPERIMENT ANALYSIS
%-----------------------------------------------------------------------------------------------------------------------
% Each Ecperiment has their own case tructure for thier whole analysis
%-----------------------------------------------------------------------------------------------------------------------

switch experiment
  case 'AAS'
%% AAS ANALYSIS

    % Logic to determine if the loop should only go to selected image for
    % locateROI or go through the whole loop for the analysis
    if locate_ROI==1
      endloop = image_selector;
    else
      endloop = 1:total_images;
    end

    for n=endloop
      
      s = s + 1;
      clc;
      % Image status and timing
      disp([num2str(s),' out of ', num2str(total_images)])
      disp(toc)
      tic

%% ROI MASK
%---------------------------------------------------------------------------------------------------
% This sections grabs information from fit files and prepares ROIs for images
%---------------------------------------------------------------------------------------------------

      % Set file names in the following format -> Pic_ID.fit
      image_name = [num2str(pic_ID(n)),'.fit'];

      % Get exposure time from images
      exposure_time(n) = getexptime(imagedir,image_name);
      


      % determines if we are getting the background image using the exposure time of each image
      if use_background_exp==1
        % rewrites filename for the case where we are selecting dark image files based on their exposure time
        image_name = [num2str(exposure_time(n)*100),'.fit'];

        % gets the exposure time for the dark image
        [texpB(n)] = getexptime(darkimagedir, image_name);

      else
        if single_dark == 1
          % gets the exposure time for the one dark image 
          [texpB(n)] = getexptime(darkimagedir, image_name(1));
        else
          % gets the exposure time for the dark image
          [texpB(n)] = getexptime(darkimagedir, image_name);
        end
      end

      % Loops over set number of ROIs
      for z = 1:numROIs                                                      
        % The factor larger that the outer RoI radius will be, compared to the inner
        widthROI(1,z)=ROIwidth+ROIscaleFactor*ROIwidth*((z-1)/(numROIs/2));
        heightROI(1,z)=ROIheight+ROIscaleFactor*ROIheight*((z-1)/(numROIs/2));
        % determines the bottom used for each RoI
        ybottom(z) = round(ROIycenter - heightROI(1,z)/2);
        % determines the top used for each RoI
        ytop(z) = round(ROIycenter + heightROI(1,z)/2);

        % determines the left value used for each RoI and each image
        xleft(z) = round(ROIxcenter - (widthROI(1,z)/2));
        % determines the right value used for each RoI and each image
        xright(z) = round(ROIxcenter + (widthROI(1,z)/2));
        xradius_ellipse(z) = (xright(z) - xleft(z))/2;
        % determines the xcenter used for each ellipse RoI and each image
        xcenter_ellipse(z) = (xright(z) + xleft(z))/2;
        % determines the yradius used for each ellipse RoI and each image
        yradius_ellipse(z) = (ytop(z) - ybottom(z))/2;
        % determines the ycenter used for each ellipse RoI and each image
        ycenter_ellipse(z) = (ybottom(z) + ytop(z))/2;
        

        % states the steps of the theta used to calculate ellipse
        ellipse_theta = 0:0.01:2*pi;

        for i = 1 : length(ellipse_theta)
          % determines the x value of the ellipse used for the RoI
          x_ellipse(z,i) = xradius_ellipse(z) * cos(ellipse_theta(i)) + xcenter_ellipse(z);
          % determines the x value of the ellipse used for the RoI
          y_ellipse(z,i) = yradius_ellipse(z) * ...
              sin(ellipse_theta(i)) + ycenter_ellipse(z);
        end
      end
   

  
      % Function that gets the date and time the image was taken from the image header
      image_time(n) = getimagetime(imagedir,image_name);
      image_duration_time(n) = image_time(n) - image_time(1);
      image_temp(n) = interp1(duration_time, temperature, image_duration_time(n),'nearest','extrap');

      if use_power_monitor==1
        image_power(n) = interp1(duration_time, power_corrected, image_duration_time(n), 'nearest','extrap') ;
        % takes the Pmon values from the thorlabs power monitor and chooses the 
        % ones measured closest to when the image was taken using interp1
        % power_corrected(n) = polyval(correction_coefficients,power_time(n));
      end
      
      % Read scattering image pixel values as doubles
      image.fit = double(fitsread([imagedir,char(image_name)]));           
      if lineremoval == 1
          image.fit(:,1570) = image.fit(:,1569);
      end
      if use_background_exp==1
        % Read dark image pixel values as doubles
        image_dark.fit = double(fitsread([darkimagedir,num2str(exposure_time(n)*100),char('.fit')]));

      else
        % Read dark image pixel values as doubles
        image_dark.fit = double(fitsread([darkimagedir,char(image_name)]));

        % The AAS experiment type will take a region of the image and
        % calculate a scaling factor for the dark images 
        % The line below is the sum of the region for the dark and
        % bright images which will be used in the subtracted image
        % calculations below

        Region_x1 = Region_xcenter - Region_xwidth;
        Region_x2 = Region_xcenter + Region_xwidth;
        Region_y1 = Region_ycenter - Region_ywidth;
        Region_y2 = Region_ycenter + Region_ywidth;

        % Get Region counts
        % Crop background image size
        dark_Region.fit = image_dark.fit(Region_y1:Region_y2,Region_x1:Region_x2,:);
        bright_Region.fit = image.fit(Region_y1:Region_y2,Region_x1:Region_x2,:);

        % bright_mubRegion = mean(mean(bright_Region.fit));
        bright_Region_Counts = sum(sum(bright_Region.fit));

        % dark_mubRegion = mean(mean(dark_Region.fit)); 
        dark_Region_Counts = sum(sum(dark_Region.fit));

        Factor = bright_Region_Counts/dark_Region_Counts;
      end

      img.fit = image.fit - (image_dark.fit*Factor);

      for z = 1:numROIs
        % creates a mask that will be used to overlay with the image for analysis
        mask(z).fit = double(poly2mask(squeeze(x_ellipse(z,:)), squeeze(y_ellipse(z,:)), 4096, 4096));
        % Create a cropped image of mask
        mask_crop(z).fit = mask(z).fit(ybottom(z):ytop(z),xleft(z):xright(z),:);
        % Create a cropped image of only the Region Of Interest
        RoI(z).fit = img.fit(ybottom(z):ytop(z),xleft(z):xright(z),:);
        % By multiplying the RoI and mask_crop it makes anything outside of the ellipse 0
        RoI_mask(z).fit = RoI(z).fit .* mask_crop(z).fit;
        % Finds the number of non-zero (nnz) pixels or cells in the RoI_mask which is used as the area
        area(n,z) = nnz(RoI_mask(z).fit);
      end

%% BRDF CALCULATIONS
%---------------------------------------------------------------------------------------------------
% This section does the calculations for BRDF
%---------------------------------------------------------------------------------------------------

      % This statement decides how background will be subtracted. From
      % background images or from defining a background ROI.

      for z = 1:numROIs
        % Calculate ARBccd summing over ROI pixel values and divide by exposure time
        ARBccd(z,n) = sum(sum(RoI_mask(z).fit))./(exposure_time(n));
      end
        
      % Calculate BRDF,PSD,and Sigma for each image BRDF calculation is based on Pmon_cal. Pmon_cal is the calibrated
      % monitor power calculated from power monitor measurements and adjusted to an incident power using ratios.
      % For TRS if powermonitor==0, ie there is no power monitor data, Pmon_cal is just incidentPower, the initial
      % measured incident power
      for z = 1:numROIs
        % B [1/str] Transpose of ARBccd was taken so that it could be divided by incidentPower. transmission of room light filter is accounted for here
        Bccd(n,z) = ARBccd(z,n)'.*muFcal./power_corrected(n)./filter(n);
        
       % BRDF [1/str]
       BRDFccd(n,z) = Bccd(n,z)./cosd(theta_s);
      end

      uncert_BRDFccd = BRDFccd*sys_uncert;
      
      % Fit the counts enclosed by each RoI vs the Area of each RoI with a
      % line y=mx+b where y is the counts enclosed, and x is the area
      % enclosed
      % Enclosed counts (actually BRDF associated with) that RoI
      y = BRDFccd(n,:);
      % Area enclosed by RoI
      x = area(n,:);
      
      % Uses MATLAB curve-fitting toolbox and creates a least-squares fit model to fit data values
      mdl = fitlm(x',y');
      % Estimated y-intercept from fit
      b(n) = mdl.Coefficients{1,{'Estimate'}};
      % Standard error for y-intercept from fit
      SE(n) = mdl.Coefficients{1,{'SE'}};
      % Determines 95% confidence bounds for the y-intercept value (1.96 is the width of the std dev curve) [ref script: linear_fit_errors.m]
      ci(n) = SE(n).*1.96;
      
      % Set the BRDF value to the y-intercept of the linear fit
      if numROIs==1
          BRDFyint(n) = y;
      else
          BRDFyint(n) = b(n);
      end
      
      % Calculates total uncertainty in measurement by taking quadrature sum of systematic and statistical errors (assumes errors are uncorrelated)
      tot_uncert(n) = sqrt(sys_uncert^2 + SE(n).^2);
      
      % Calculates the absolute uncertainty in BRDF measurement
      uncert_BRDFyint(n) = b(n)*tot_uncert(n);
      
      switch experiment
        case 'ARS'
          % Calculate PSD_fx_fy for each point, brdf_to_psd expects wavelength in nm
          [PSD_fx_fy(n),S_iso(n),f_quad_sum(n)] = brdf_to_psd(theta_s(n),phi_s(n),correction_angle,lambda*1e9,BRDFyint(n)');
          
          % calculate sigma (can only do so if we have second point already)
          if n>1
            [sigma_rel, sigma_squared, sigma_cumulative] = psd_to_sigma(f_quad_sum, S_iso);
          end
      end
      
      % Determine BRDF limit based on darkest small box in image number of divisions per line
      %---------------------------------------------------------------------------------------------

      % Finds the x and y coordinates for all test box centers
      for N = 1:num_of_div
        LIMy(N)=N*(length(img.fit(:,1))-100)/num_of_div;
        LIMx(N)=N*(length(img.fit(1,:))-100)/num_of_div;
      end
      
      % Creates a counter to store BRDF information for each box
      LIMcount = 1;
      
      
      for Y = 1:num_of_div
        % make a box
        for X = 1:num_of_div
          LIMx1(LIMcount) = round(LIMx(X)-LIMwidth/2);
          LIMx2(LIMcount) = round(LIMx(X)+LIMwidth/2);
          LIMy1(LIMcount) = round(LIMy(Y)-LIMheight/2);
          LIMy2(LIMcount) = round(LIMy(Y)+LIMheight/2);

          % Crop beam image size
          LIMpic(LIMcount).fit = img.fit(LIMy1(LIMcount):LIMy2(LIMcount),LIMx1(LIMcount):LIMx2(LIMcount),:);

          % sum over ROI and correct
          LIM_ARB(LIMcount) = sum(sum(LIMpic(LIMcount).fit))./(exposure_time(n));

          % B [1/sr] Transpose of ARBccd was taken so that it could be divided
          % by incidentPower. transmission of room light filter is accounted for here
          LIM_B(LIMcount) = LIM_ARB(LIMcount)'.*muFcal./power_corrected(n)./filter(n);

          % BRDF [1/sr]
          LIM_BRDF(LIMcount) = LIM_B(LIMcount) ./ cosd(theta_s) * ((xright(1)-xleft(1)) * (ytop(1)-ybottom(1)) / (LIMwidth*LIMheight));
          
          LIMcount = LIMcount+1;
        end
      end
      
      %Find smallest LIM_BRDF and plot its box in image.
      [minLIM(n),minLIMpos] = min(LIM_BRDF);
      LIMITx1 = LIMx1(minLIMpos);
      LIMITx2 = LIMx2(minLIMpos);
      LIMITy1 = LIMy1(minLIMpos);
      LIMITy2 = LIMy2(minLIMpos);
      %END of BRDF limitations

%% (DEBUG) LOCATE ROI IMAGE PLOTTING
%---------------------------------------------------------------------------------------------------
% This section plots the images that are uses for locating the ROI before the analysis starts
%---------------------------------------------------------------------------------------------------

      if locate_ROI == 1
        % Display the corrected image to determine wehre the ROI min's and
        % max's are, and input them accordingly

        % Get the width and height of the CCD image (pixel integer)
        img_x = 1:1:width(img.fit);
        img_y = 1:1:height(img.fit);

        % Pixel integer shifted so zero in center (integer)
        img_x_centered = img_x-length(img_x)/2; 
        img_y_centered = img_y-length(img_y)/2;

        % Convert integer scales into mm scales using calibrated value
        img_x_mm = img_x_centered / pix_per_mm;
        img_y_mm = img_y_centered / pix_per_mm;
    
        % Plot the image and RoIs
        figure;
        ax1 = gca;

        % Plot the CCD image with the converted x and y values
        imagesc(img_x_mm,img_y_mm,img.fit);
        
        hold on;

        boxcolors = {[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1]};
        
        % sized by the number of ROIs currently 6   
        linewidth = {.5, .1, .1, .1, .1, .1};
        linestyle = {'-',':',':',':',':',':'};
          

        for z = 1:numROIs
          plot(squeeze(x_ellipse(z,:)), squeeze(y_ellipse(z,:)),'LineStyle',linestyle{z},'color',boxcolors{z}, ...
               'LineWidth',linewidth{z});

          plot([Region_x1 Region_x1 Region_x2 Region_x2 Region_x1], [Region_y1 Region_y2 Region_y2 Region_y1 Region_y1],'color','m','linewidth',1);
        end

        text(200, 3886, [['Image Number: ',num2str(image_selector)],' \newline',['Exposure Time = ', ...
             num2str(exposure_time(n)),' s'], ' \newline', ['BRDF = ',num2str(sprintf('%.2s',BRDFyint(n))),' 1/str'],...
             ' \newline', ['Background Factor = ', num2str(sprintf('%.2f',Factor)),' (Bright Count / Dark Count)']],...
             'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10, 'Color', [1-eps 1 1]);

        % CCD image properties
        colorbar;
        colormap('gray');
        ax1.CLim = [clim_Min, clim_Max];
        % ax1.XTick = -8:1:8;
        % ax1.YTick = -8:1:8;
        ax1.DataAspectRatio = [1 1 1];

        return
      end

%% FIGURE 1: IMAGE, HISTOGRAM, & BRDF SUBPLOTS

%---------------------------------------------------------------------------------------------------
% This section plots the image, histogram, and BRDF subplots
%---------------------------------------------------------------------------------------------------

      % Figure with image, BRDF plot, Histogram, PSD, and info.
      %---------------------------------------------------------------------------------------------------
    
      f1 = figure('Position',[1 1000 1200 1000],'Visible','off');
      set(0,'defaulttextfontname','Times New Roman')
      set(0,'defaultaxesfontname','Times New Roman')
        
      
      % Plot the CCD image with its ROIs
      %---------------------------------------------------------------------------------------------------
      ax1 = nexttile([2,2]);

      % Get the width and height of the CCD image (pixel integer)
      img_x = 1:1:width(img.fit);
      img_y = 1:1:height(img.fit);

      % Pixel integer shifted so zero in center (integer)
      img_x_centered = img_x-length(img_x)/2;
      img_y_centered = img_y-length(img_y)/2;

      % Convert integer scales into mm scales using calibrated value
      img_x_mm = img_x_centered / pix_per_mm;
      img_y_mm = img_y_centered / pix_per_mm;
      
      imagesc(img_x_mm,img_y_mm,img.fit);
      
      hold on;
      ax1.CLim = [clim_Min, clim_Max];
      colormap('gray');
      boxcolors = {[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1]};
      
      % sized by the number of ROIs currently 6   
      linewidth = {.5, .1, .1, .1, .1, .1};
      linestyle = {'-',':',':',':',':',':'};
      

      for z = 1:numROIs
         plot(((squeeze(x_ellipse(z,:))) ./ pix_per_mm) - pixel_x_shift,...
               ((squeeze(y_ellipse(z,:))) ./ pix_per_mm) - pixel_y_shift,...
               'LineStyle',linestyle{z},'color',boxcolors{z},'LineWidth',linewidth{z});
      end

      title([sample,', image ',char(image_name)], 'FontSize',20,'FontName','Times New Roman','Interpreter', 'none');
      axis square;
      axis off;

      % Plot histogram
      ax2 = nexttile([2,2]);

      % take matrix RoI(n).fit and turn into a giant column vector
      imagevector = reshape(RoI_mask(1).fit,numel(RoI_mask(1).fit),1);
         
      % Make the bin center values where 65600 is calculated how????????????
      bincentervalue = 0:100:65600;
      
      % Make the histogram
      h = histogram(imagevector,0:100:65700);
      
      % Extract the counts from the histogram
      numinbin = h.Values;
      % make a bar plot of log10 of the number for each bin and gets handle
      bh = bar(bincentervalue,numinbin);
      
      hold on;
      % some made up power law trend to be understood later
      histtheory = 3e8*bincentervalue.^(-2);
      
      plot(bincentervalue,histtheory,'k--');
      % x-axis limit
      ax2.XLim = [0 65600];
      % y-axis limit
      ax2.YLim = [min(bincentervalue) max(bincentervalue)];
      ax2.FontSize = 18;
      ax2.XScale = 'log';
      ax2.YScale = 'log';
      bh.FaceColor = [1 0 0];
      bh.EdgeColor = [1 0 0];
      ylabel(ax2,'# of pixels','FontName','Times New Roman');
      xlabel(ax2,'Pixel value [counts]','FontName','Times New Roman');
      title('Pixel values in region of interest','FontSize',20,'FontName','Times New Roman')
      
      
      % Create text information to add to plot
      ELAPSE_TIME = ['Elapsed Time = ',char(image_duration_time(n))];
      POWERINC = ['Incident power = ',num2str(sprintf('%.2f',power_corrected(n)*1000)), ' mW'];
      POWERSCAT = ['Scattered power = ',num2str(sprintf('%.2s',ARBccd(1).*muFcal.*1000)),' {\mu}W'];
      BRDFVALUE = ['BRDF = ',num2str(sprintf('%.2s',BRDFyint(n))),' 1/str'];
      TEMP = ['Temperature = ', num2str(round(image_temp(n),1)), '{\circ}C' ];
      
      % Add text information to plot
      text(65600, max(bincentervalue), [' \newline',ELAPSE_TIME,...
          '\newline',POWERINC,'\newline',POWERSCAT,'\newline',BRDFVALUE],...
          'VerticalAlignment','top','HorizontalAlignment','right','FontSize',12)

      % Plot BRDF vs Time and Temperature vs Time on 3-pane subplot
      ax3 = nexttile([2,4]);
      yyaxis left

      plot(image_duration_time(1:n),BRDFyint(1:n),'LineStyle','-');
      set(ax3, 'YScale', 'log');
      hold on;
      
      ylabel('BRDF [1/str]','FontSize',18,'FontName','Times New Roman');
      xlabel('Time Elapsed','FontSize',18,'FontName','Times New Roman');
      
      yyaxis right
      plot(image_duration_time(1:n),image_temp(1:n));

      ylabel('Temperature [{\circ}C]','FontSize',18);
      grid on
      grid minor
      set(ax3, 'FontSize',16);            
      title(['Normalized Scatter and Temperature Profile ', strrep(folder.sample_name,'_','\_')] , 'FontSize',24)
      
    
      % Save 3-pane subplot to fig only if last image
      if n == total_images
          saveas(f1,[folder.analysisPath,'FIG',slash,num2str(pic_ID(n)),'.fig']);
      end

      figure_size = get(gcf, 'position');
      set(f1,'PaperPosition',figure_size/100, 'visible', 'off');

      if print_images == 1
        % Save 3-pane subplot frame as png for movie making
        print(f1,[folder.analysisPath,'PNG',slash,num2str(pic_ID(n)),'.png'], '-dpng','-r300');
      end
      delete(f1);
                
      

%% FIGURE 2: CCD IMAGES
%---------------------------------------------------------------------------------------------------
% Create PNGs of just the CCD images
%---------------------------------------------------------------------------------------------------
      if print_images == 1

        f2 = figure('Visible','off');
        ax5= gca;

        imagesc(img_x_mm,img_y_mm,img.fit);
        
        hold on;
        % Adjusts Colormap limits [cMin,cMax] on output images, change Brightness
        ax5.CLim = [clim_Min, clim_Max];
        colormap('gray');
        boxcolors = {[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1]};
        linewidth = {.5, .1, .1, .1, .1, .1};
        linestyle = {'-',':',':',':',':',':'};
        
        for z = 1:numROIs
           plot(((squeeze(x_ellipse(z,:))) ./ pix_per_mm) - pixel_x_shift,...
               ((squeeze(y_ellipse(z,:))) ./ pix_per_mm) - pixel_y_shift,...
               'LineStyle',linestyle{z},'color',boxcolors{z},'LineWidth',linewidth{z});
        end

        TITLE = [sample, ' Image ', char(image_name)];

        text((2850 / pix_per_mm) - pixel_x_shift, (3886 / pix_per_mm) - pixel_y_shift, [' \newline',TITLE,' \newline',ELAPSE_TIME,...
            '\newline',POWERINC,'\newline',POWERSCAT,'\newline',BRDFVALUE, '\newline', TEMP],...
            'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',11, 'Color', [1-eps 1 1]);

        title([folder.sample_name,', image ',char(image_name)], 'FontSize',20,'FontName','Times New Roman','Interpreter', 'none');

        set(f2,'units','points','position',[0 0 480 480]); % make it a square
        
        set(ax5,'position',[0 0 1 1]) % make square axes fill the figure
        axis square;
        axis on;

  
        print(f2,[folder.analysisPath,'CCD_PNG',slash,num2str(pic_ID(n)),'.png'], '-dpng','-r300');

        delete(f2);
      
  
        for z = 1:numROIs
          % Sets the RoI value in the structure to null, replacing the 4096x4096
          RoI(z).fit = 0;
          
          % Sets the mask value in the structure to null, replacing the 4096x4096
          mask(z).fit = 0;
          
          % Overwrites the values of each element to null
          mask_crop(z).fit = 0;
          
          % Overwrites the values of each element to null
          RoI_mask(z).fit = 0;
        end
        
        % Overwrites the values of each element to null
        image.fit = 0;
        
        % Overwrites the values of each element to null
        image_dark.fit = 0;
        
        % Overwrites the values of each element to null
        img.fit = 0;
      end
    end
    
%% FIGURE 3: BRDF VS TTIME
%---------------------------------------------------------------------------------------------------
% Plot and save BRDF vs theta_s and BRDFlimit
%---------------------------------------------------------------------------------------------------
    
    if print_images == 1
      f3 = figure('Visible','off');
      ax3 = gca;
     
      semilogy(image_duration_time,minLIM','v-','DisplayName','BRDF limit');
      hold on;

      semilogy(image_duration_time,BRDFyint,'s','MarkerFaceColor','g','LineStyle','none','Color','g','DisplayName',sample);
      xlabel('Elapsed Time (hh:mm:ss)');
      
      ylabel('BRDF [1/str]');
      BRDF_title = ['BRDF ',sample];
      title(BRDF_title);
  
  %---------------------------------------------------------------------------------------------------
  % Figure 3 Properties
  %---------------------------------------------------------------------------------------------------
      ax3.FontName = 'Times New Roman';
      ax3.FontSize = 30;
      f3.Position = [0,0,1618,1000];
  
      saveas(f3,[folder.analysisPath,folder.sample_name,slash,'BRDF_',folder.timeAnalysis,'_', folder.sample_name,'.fig']);
      saveas(f3,[folder.analysisPath,folder.sample_name,slash,'BRDF_',folder.timeAnalysis,'_', folder.sample_name,'.png']);
      delete(f3);
    end

%% FIGURE 4: INCIDENT POWER

    if print_images == 1
      % Plots of incident power over scattering angle
      f4 = figure('Position',[0 0 1000 1000],'Visible','off');
      hold on;
      plot(pic_ID, power_corrected,'s-','MarkerFaceColor','b');
      hold on;
      plot(pic_ID, power_corrected,'s-','MarkerFaceColor','r');
      xlabel('Time', 'FontSize', 20);
      title('Incident power over time','FontSize',30,'FontName','Times New Roman','Interpreter','none');
      ylabel('Power [W]', 'FontSize', 20);
      legend('Incident Power (measured twice, before and after)', 'Calibrated Monitor Power');
      grid('on');
      box('on');
      orient landscape;
      
      saveas(f4, [folder.analysisPath,folder.sample_name,slash,'Calibrated_monitor_power_',folder.timeAnalysis,'_', folder.sample_name,'.fig']);
      saveas(f4, [folder.analysisPath,folder.sample_name,slash,'Calibrated_monitor_power_',folder.timeAnalysis,'_', folder.sample_name,'.png']);
      delete(f4);
    end

%% FIGURE 5: HEATER POWER
    
    % Create Temperature/Power Graph from excel file

    f5_a = figure('Visible','off');
    yyaxis left
    plot(duration_time,temperature,'b-')
    hold on
    % ylim([0 max(temperature)+50]) ;  %edit limit depending on max temp, 0 to 500 for 450C exp
    ylabel('Temperature [{\circ}C]','FontSize',14)
    set(gca,'FontSize',14)%20
    
    yyaxis right
    plot(duration_time,heater_power,'r-')
    ylabel('Output Power [%]','FontSize',14)%18
    % ylim([0 max(power)+10]) ;
    set(gca,'FontSize',14)%20
    
    xlabel('Time (hr)','FontSize',18)
    % legend('Temperature [{\circ}C]','Output Power [%]','FontSize',14)
    title(['Temperature-Power Profile ', strrep(folder.sample_name,'_','\_')],'FontSize',18)%22
    grid on
    grid minor
    saveas(f5_a, [folder.analysisPath,folder.sample_name,slash, 'Temp_Power_', folder.sample_name,'_', folder.timeAnalysis, '.fig']);
    saveas(f5_a, [folder.analysisPath,folder.sample_name,slash, 'Temp_Power_', folder.sample_name,'_', folder.timeAnalysis, '.png']);
    delete(f5_a);
    
    % create a monitored power vs time plot that shows the entire power
    % monitor measurements, not just those from the interpolation
    
    if use_power_monitor==1
      f5_b = figure('Visible','off');
      plot(date_time, power_corrected, 's-', 'MarkerFaceColor', 'r')
      xlabel('Elapsed Time', 'FontSize', 20);
      ylabel('Power [W]', 'FontSize', 20);
      title('Monitored Power Over Time','FontSize',30,'FontName','Times New Roman','Interpreter','none');
      saveas(f5_b, [folder.analysisPath,folder.sample_name,slash,'Monitor_power_',folder.timeAnalysis,'_', folder.sample_name,'.fig']);
      saveas(f5_b, [folder.analysisPath,folder.sample_name,slash,'Monitor_power_',folder.timeAnalysis,'_', folder.sample_name,'.png']);
      delete(f5_b);
           
      % Make a graph of power monitor measurements alongside heater power percentage to check correlation
      f5_c = figure('Visible','off');
      yyaxis right
      plot(image_duration_time, power_corrected, 's-','MarkerFaceColor','r')
      ylabel('Power [W]', 'FontSize', 20)
      yyaxis left
      plot(duration_time, heater_power, 'b-')
      ylabel('Output Power [%]', 'FontSize', 20)
      legend('Calibrated Monitor Power', 'Heater Power')
      xlabel('Elapsed Time',  'FontSize', 20)
      title('Laser and Heater Power During Annealing Run', 'FontSize',20)
      saveas(f5_c, [folder.analysisPath,folder.sample_name,slash,'heater_laser_power_',folder.timeAnalysis,'_', folder.sample_name,'.fig']);
      saveas(f5_c, [folder.analysisPath,folder.sample_name,slash,'heater_laser_power_',folder.timeAnalysis,'_', folder.sample_name,'.png']);
      delete(f5_c);
    end

%% FIGURE 6: LASER POWER PROFILE

    % Laser Power Profile plots
    f6 = figure('Visible','off');
    yyaxis left
    plot(image_duration_time,BRDFyint,'b-')
    hold on
    ylabel('BRDF [1/str]','FontSize',18)
    set(gca,'FontSize',12)
    yyaxis right
    plot(image_duration_time,power_corrected,'r-')
    ylabel('Laser Power [W]','FontSize',18)
    set(gca,'FontSize',12)
    xlabel('Time (hh:mm:ss)','FontSize',14)

    legend('BRDF [1/str]','Laser Power [W]','FontSize',14)
    title(['BRDF-LaserPower Profile ', strrep(folder.sample_name,'_','\_')],'FontSize',18)
    grid on
    grid minor
    saveas(f6, [folder.analysisPath,folder.sample_name,slash,'BRDF_LaserPower_', folder.sample_name,'_', folder.timeAnalysis, '.fig']);
    saveas(f6, [folder.analysisPath,folder.sample_name,slash,'BRDF_LaserPower_', folder.sample_name,'_', folder.timeAnalysis, '.png']);
    delete(f6);


  case 'ARS'
%% ARS ANALYSIS

    % Logic to determine if the loop should only go to selected image for
    % locateROI or go through the whole loop for the analysis
    if locate_ROI==1
      endloop = image_selector;
    else
      endloop = 1:total_images;
    end

    for n=endloop
      
      s = s + 1;
      clc;
      % Image status and timing
      disp([num2str(s),' out of ', num2str(total_images)])
      disp(toc)
      tic

%% ROI MASK
%---------------------------------------------------------------------------------------------------
% This sections grabs information from fit files and prepares ROIs for images
%---------------------------------------------------------------------------------------------------

      % Set file names in the following format -> Pic_ID.fit
      image_name = [num2str(pic_ID(n)),'.fit'];

      % Get exposure time from images
      exposure_time(n) = getexptime(imagedir,image_name);
                                                            
      % Read exposure times for the dark images

      % determines if we are getting the background image using the exposure time of each image
      if use_background_exp==1
        % rewrites filename for the case where we are selecting dark image files based on their exposure time
        image_name = [num2str(exposure_time(n)*100),'.fit'];

        % gets the exposure time for the dark image
        [texpB(n)] = getexptime(darkimagedir, image_name);
      else

      % gets the exposure time for the dark image
      [texpB(n)] = getexptime(darkimagedir, image_name);
      end

      % Loops over set number of ROIs
      for z = 1:numROIs                                                      
        if z==1
          % determines the width for each RoI
          widthROI(1,z) = round(2*ROIwidth*(z)/(z + 1));
          % determines the height used for each RoI
          heightROI(1,z) = round(2*ROIheight*(z)/(z + 1));
          % determines the bottom used for each RoI
          ybottom(z) = round(ROIycenter - heightROI(1,z)/2);
          % determines the top used for each RoI
          ytop(z) = round(ROIycenter + heightROI(1,z)/2);
        else
          % determines the width for each RoI
          widthROI(1,z) = round(2*ROIwidth*(z)/(z + 1.5));
          % determines the height used for each RoI
          heightROI(1,z) = round(2*ROIheight*(z)/(z + 1.5));
          % determines the bottom used for each RoI
          ybottom(z) = round(ROIycenter - heightROI(1,z)/2);
          % determines the top used for each RoI
          ytop(z) = round(ROIycenter + heightROI(1,z)/2);
        end

        % Define the RoI that contains the beam (changing RoI)
        % determines the left value used for each RoI and each image
        xleft(z) = round((ROIxcenter + df*sind(theta_rot_stage(n))) - ((widthROI(1,z)*cosd(theta_s(n)))/2));

        % determines the right value used for each RoI and each image
        xright(z) = round((ROIxcenter + df*sind(theta_rot_stage(n)))+((widthROI(1,z)*cosd(theta_s(n)))/2));

        % Determines the xradius, xcenter, yradius, ycenter used for
        % each ellipse RoI and each image
        xradius_ellipse(z) = (xright(z) - xleft(z))/2;
        xcenter_ellipse(z) = (xright(z) + xleft(z))/2;
        yradius_ellipse(z) = (ytop(z) - ybottom(z))/2;
        ycenter_ellipse(z) = (ybottom(z) + ytop(z))/2;


        % states the steps of the theta used to calculate ellipse
        ellipse_theta = 0:0.01:2*pi;

        for i = 1 : length(ellipse_theta)
          % determines the x value of the ellipse used for the RoI
          x_ellipse(z,i) = xradius_ellipse(z) * cos(ellipse_theta(i)) + xcenter_ellipse(z);
          % determines the x value of the ellipse used for the RoI
          y_ellipse(z,i) = yradius_ellipse(z) * ...
              sin(ellipse_theta(i)) + ycenter_ellipse(z);
        end
      end
   
      % Read scattering image pixel values as doubles
      image.fit = double(fitsread([imagedir,char(image_name)]));           
      if lineremoval == 1
          image.fit(:,1570) = image.fit(:,1569);
      end
      if use_background_exp==1
        % Read dark image pixel values as doubles
        image_dark.fit = double(fitsread([darkimagedir,num2str(exposure_time(n)*100),char('.fit')]));
      else
  

      % Read dark image pixel values as doubles
      image_dark.fit = double(fitsread([darkimagedir,char(image_name)]));

      % The AAS experiment type will take a region of the image and
      % calculate a scaling factor for the dark images 
      % The line below is the sum of the region for the dark and
      % bright images which will be used in the subtracted image
      % calculations below

      Region_x1 = Region_xcenter - Region_xwidth;
      Region_x2 = Region_xcenter + Region_xwidth;
      Region_y1 = Region_ycenter - Region_ywidth;
      Region_y2 = Region_ycenter + Region_ywidth;

      % Get Region counts
      % Crop background image size
      dark_Region.fit = image_dark.fit(Region_y1:Region_y2,Region_x1:Region_x2,:);
      bright_Region.fit = image.fit(Region_y1:Region_y2,Region_x1:Region_x2,:);

      % bright_mubRegion = mean(mean(bright_Region.fit));
      bright_Region_Counts = sum(sum(bright_Region.fit));

      % dark_mubRegion = mean(mean(dark_Region.fit)); 
      dark_Region_Counts = sum(sum(dark_Region.fit));

      Factor = bright_Region_Counts/dark_Region_Counts;
          
      end

      img.fit = image.fit - (image_dark.fit*Factor);

      for z = 1:numROIs
        % creates a mask that will be used to overlay with the image for analysis
        mask(z).fit = double(poly2mask(squeeze(x_ellipse(z,:)), squeeze(y_ellipse(z,:)), 4096, 4096));

        % Create a cropped image of mask
        mask_crop(z).fit = mask(z).fit(ybottom(z):ytop(z),xleft(z):xright(z),:);

        % Create a cropped image of only the Region Of Interest
        RoI(z).fit = img.fit(ybottom(z):ytop(z),xleft(z):xright(z),:);

        % By multiplying the RoI and mask_crop it makes anything outside of the ellipse 0
        RoI_mask(z).fit = RoI(z).fit .* mask_crop(z).fit;

        % Finds the number of non-zero (nnz) pixels or cells in the RoI_mask which is used as the area
        area(n,z) = nnz(RoI_mask(z).fit);
      end

%% BRDF CALCULATIONS
%---------------------------------------------------------------------------------------------------
% This section does the calculations for BRDF
%---------------------------------------------------------------------------------------------------

      % This statement decides how background will be subtracted. From
      % background images or from defining a background ROI.
        
      for z = 1:numROIs
        % Calculate ARBccd summing over ROI pixel values and divide by exposure time
        ARBccd(z,n) = sum(sum(RoI_mask(z).fit))./(exposure_time(n));
      end
        
      % Calculate BRDF,PSD,and Sigma for each image
      % BRDF calculation is based on Pmon_cal. Pmon_cal is the calibrated
      % monitor power calculated from power monitor measurements and adjusted
      % to an incident power using ratios. For TRS if powermonitor==0, ie
      % there is no power monitor data, Pmon_cal is just incidentPower, the initial
      % measured incident power
      for z = 1:numROIs
        % B [1/str] Transpose of ARBccd was taken so that it could be divided by incidentPower. transmission of room light filter is accounted for here
        Bccd(n,z) = ARBccd(z,n)'.*muFcal./power_corrected(n)./filter(n);

        % BRDF [1/str]
        BRDFccd(n,z) = Bccd(n,z)./cosd(theta_s(n));

      end
      uncert_BRDFccd = BRDFccd*sys_uncert;
      
      % Fit the counts enclosed by each RoI vs the Area of each RoI with a
      % line y=mx+b where y is the counts enclosed, and x is the area
      % enclosed
      % Enclosed counts (actually BRDF associated with) that RoI
      y = BRDFccd(n,:);
      % Area enclosed by RoI
      x = area(n,:);
      
      % Uses MATLAB curve-fitting toolbox and creates a least-squares fit model to fit data values
      mdl = fitlm(x',y');
      % Estimated y-intercept from fit
      b(n) = mdl.Coefficients{1,{'Estimate'}};
      % Standard error for y-intercept from fit
      SE(n) = mdl.Coefficients{1,{'SE'}};
      % Determines 95% confidence bounds for the y-intercept value (1.96 is the width of the std dev curve) [ref script: linear_fit_errors.m]
      ci(n) = SE(n).*1.96;
      
      % Set the BRDF value to the y-intercept of the linear fit
      if numROIs==1
          BRDFyint(n) = y;
      else
          BRDFyint(n) = b(n);
      end
      
      % Calculates total uncertainty in measurement by taking quadrature sum of systematic and statistical errors (assumes errors are uncorrelated)
      tot_uncert(n) = sqrt(sys_uncert^2 + SE(n).^2);
      
      % Calculates the absolute uncertainty in BRDF measurement
      uncert_BRDFyint(n) = b(n)*tot_uncert(n);
      
      % Calculate PSD_fx_fy for each point, brdf_to_psd expects wavelength in nm
      [PSD_fx_fy(n),S_iso(n),f_quad_sum(n)] = brdf_to_psd(theta_s(n),phi_s(n),correction_angle,lambda*1e9,BRDFyint(n)');
      
      % calculate sigma (can only do so if we have second point already)
      if n>1
        [sigma_rel, sigma_squared, sigma_cumulative] = psd_to_sigma(f_quad_sum, S_iso);
      end
      
      % Determine BRDF limit based on darkest small box in image number of divisions per line
      %---------------------------------------------------------------------------------------------

      % Finds the x and y coordinates for all test box centers
      for N = 1:num_of_div
        LIMy(N)=N*(length(img.fit(:,1))-100)/num_of_div;
        LIMx(N)=N*(length(img.fit(1,:))-100)/num_of_div;
      end
      
      % Creates a counter to store BRDF information for each box
      LIMcount = 1;
      
      
      for Y = 1:num_of_div
        % make a box
        for X = 1:num_of_div
          LIMx1(LIMcount) = round(LIMx(X)-LIMwidth/2);
          LIMx2(LIMcount) = round(LIMx(X)+LIMwidth/2);
          LIMy1(LIMcount) = round(LIMy(Y)-LIMheight/2);
          LIMy2(LIMcount) = round(LIMy(Y)+LIMheight/2);
          % Crop beam image size
          LIMpic(LIMcount).fit = img.fit(LIMy1(LIMcount):...
              LIMy2(LIMcount),LIMx1(LIMcount):LIMx2(LIMcount),:);
          % sum over ROI and correct
          LIM_ARB(LIMcount) = sum(sum(LIMpic(LIMcount).fit))./(exposure_time(n));
          % B [1/sr] Transpose of ARBccd was taken so that it could be divided
          % by incidentPower. transmission of room light filter is accounted for here
          LIM_B(LIMcount) = LIM_ARB(LIMcount)'.*muFcal...
              ./power_corrected(n)./filter(n);

          % BRDF [1/sr]
          LIM_BRDF(LIMcount) = LIM_B(LIMcount)./cosd(theta_s(n)) * ((xright(1)-xleft(1)) * (ytop(1)-ybottom(1))/(LIMwidth*LIMheight));
          
          LIMcount = LIMcount+1;
        end
      end
      
      %Find smallest LIM_BRDF and plot its box in image.
      [minLIM(n),minLIMpos] = min(LIM_BRDF);
      LIMITx1 = LIMx1(minLIMpos);
      LIMITx2 = LIMx2(minLIMpos);
      LIMITy1 = LIMy1(minLIMpos);
      LIMITy2 = LIMy2(minLIMpos);
      %END of BRDF limitations

%% (DEBUG) LOCATE ROI IMAGE PLOTTING
%---------------------------------------------------------------------------------------------------
% This section plots the images that are uses for locating the ROI before the analysis starts
%---------------------------------------------------------------------------------------------------

      if locate_ROI == 1
        % Display the corrected image to determine wehre the ROI min's and
        % max's are, and input them accordingly

        % Get the width and height of the CCD image (pixel integer)
        img_x = 1:1:width(img.fit);
        img_y = 1:1:height(img.fit);

        % Pixel integer shifted so zero in center (integer)
        img_x_centered = img_x-length(img_x)/2; 
        img_y_centered = img_y-length(img_y)/2;

        % Convert integer scales into mm scales using calibrated value
        img_x_mm = img_x_centered / pix_per_mm;
        img_y_mm = img_y_centered / pix_per_mm;
    
        % Plot the image and RoIs
        figure;
        ax1 = gca;

        % Plot the CCD image with the converted x and y values
        imagesc(img_x_mm,img_y_mm,img.fit);
        
        hold on;

        boxcolors = {[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1]};
        
        % sized by the number of ROIs currently 6   
        linewidth = {.5, .1, .1, .1, .1, .1};
        linestyle = {'-',':',':',':',':',':'};
          

        for z = 1:numROIs
          plot(squeeze(x_ellipse(z,:)), squeeze(y_ellipse(z,:)),'LineStyle',linestyle{z},'color',boxcolors{z},'LineWidth',linewidth{z})
          plot([Region_x1 Region_x1 Region_x2 Region_x2 Region_x1], [Region_y1 Region_y2 Region_y2 Region_y1 Region_y1],'color','m','linewidth',1)
        end

        text(200, 3886, [['Image Number: ',num2str(image_selector)],...
                         ' \newline',['Exposure Time = ',num2str(exposure_time(n)),' s'], ...
                         ' \newline', ['BRDF = ',num2str(sprintf('%.2s',BRDFyint(n))),' 1/str'], ...
                         ' \newline', ['Background Factor = ', num2str(sprintf('%.2f',Factor)),' (Bright Count / Dark Count)']],...
                         'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10, 'Color', [1-eps 1 1]);

        % CCD image properties
        colorbar;
        colormap('gray');
        ax1.CLim = [clim_Min, clim_Max];
        % ax1.XTick = -8:1:8;
        % ax1.YTick = -8:1:8;
        ax1.DataAspectRatio = [1 1 1];

        return
      end

%% FIGURE 1: IMAGE, HISTOGRAM, & BRDF SUBPLOTS

%---------------------------------------------------------------------------------------------------
% This section plots the image, histogram, and BRDF subplots
%---------------------------------------------------------------------------------------------------

      % Figure with image, BRDF plot, Histogram, PSD, and info.
      %---------------------------------------------------------------------------------------------------
    
      f1 = figure('Position',[1 1000 1200 1000],'Visible','off');
      set(0,'defaulttextfontname','Times New Roman')
      set(0,'defaultaxesfontname','Times New Roman')
        
      
      % Plot the CCD image with its ROIs
      %---------------------------------------------------------------------------------------------------
      ax1 = nexttile([2,2]);

      % Get the width and height of the CCD image (pixel integer)
      img_x = 1:1:width(img.fit);
      img_y = 1:1:height(img.fit);

      % Pixel integer shifted so zero in center (integer)
      img_x_centered = img_x-length(img_x)/2; 
      img_y_centered = img_y-length(img_y)/2;

      % Convert integer scales into mm scales using calibrated value
      img_x_mm = img_x_centered / pix_per_mm;
      img_y_mm = img_y_centered / pix_per_mm;

      % Pixel shift value
      pixel_x_shift = (length(img_x)/2) / pix_per_mm;
      pixel_y_shift = (length(img_y)/2) / pix_per_mm;

      % Plot image with proper scale
      imagesc(img_x_mm, img_y_mm, img.fit);

      hold on;
      ax1.CLim = [clim_Min, clim_Max];
      colormap('gray');
      boxcolors = {[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1]};
      
      % sized by the number of ROIs currently 6   
      linewidth = {.5, .1, .1, .1, .1, .1};
      linestyle = {'-',':',':',':',':',':'};
      
        for z = 1:numROIs
          plot(((squeeze(x_ellipse(z,:))) ./ pix_per_mm) - pixel_x_shift,...
                 ((squeeze(y_ellipse(z,:))) ./ pix_per_mm) - pixel_y_shift,...
                 'LineStyle',linestyle{z},'color',boxcolors{z},'LineWidth',linewidth{z});
        end

      title([sample,', image ',char(image_name)], 'FontSize',20,'FontName','Times New Roman','Interpreter', 'none');
      axis square;
      axis off;

      %---------------------------------------------------------------------------------------------------
      % Plot the BRDF vs theta_s and BRDF limit vs theta_s
      %---------------------------------------------------------------------------------------------------
      ax2 = nexttile([2,2]);

      % Plot errorbars on y-intercept BRDF value
      errorbar(theta_s(1:n),BRDFyint(1:n),uncert_BRDFyint(1:n),'s','Color',[0 0.5 0]);

      % Chooses the x-axis display limit
      ax2.XLim = [(min(theta_s(1:n))-5) (max(theta_s(1:n))+5)];
      ax2.FontSize = 18;
      ax2.YScale = 'log';
      grid on; hold on;

      % Displays the grid on the BRDF PLOT
      ylabel('BRDF [1/str]','FontSize',18,'FontName','Times New Roman');
      xlabel('\theta_s { [degs.]}','FontSize',18,'FontName','Times New Roman');
      hold on;
      semilogy(theta_s(1:n),minLIM(1:n),'v-');
      legend({'BRDF','BRDF Limit'}, 'FontSize',18,'FontName','Times New Roman');
      title('Normalized scatter','FontSize',20,'FontName','Times New Roman')
      
      %---------------------------------------------------------------------------------------------------
      % Plot histogram
      %---------------------------------------------------------------------------------------------------
      ax3 = nexttile([2,2]);

      % take matrix RoI(n).fit and turn into a giant column vector
      imagevector = reshape(RoI_mask(1).fit,numel(RoI_mask(1).fit),1);

      % Make the bin center values where 65600 is calculated how????????????
      bincentervalue = 0:100:65600;
      
      % Make the histogram
      h = histogram(imagevector,0:100:65700);
      
      % Extract the counts from the histogram
      numinbin = h.Values;

      % make a bar plot of log10 of the number for each bin and gets handle
      bh = bar(bincentervalue,numinbin);
      
      hold on;
      % some made up power law trend to be understood later
      histtheory = 3e8*bincentervalue.^(-2);
      
      plot(bincentervalue,histtheory,'k--');
      % gets handles of axes
      
      % x-axis limit
      ax3.XLim = [0 65600];
      % y-axis limit
      ax3.YLim = [min(bincentervalue) max(bincentervalue)];
      ax3.FontSize = 18;
      ax3.XScale = 'log';
      ax3.YScale = 'log';
      bh.FaceColor = [1 0 0];
      bh.EdgeColor = [1 0 0];
      ylabel(ax3,'# of pixels','FontName','Times New Roman');
      xlabel(ax3,'Pixel value [counts]','FontName','Times New Roman');
      title('Pixel values in region of interest','FontSize',20,'FontName','Times New Roman')
      
      % Create text information to add to plot
      ANGLEtheta_s=['Angle = ',num2str(theta_s(n)), '\circ'];
      POWERINC=['Incident power = ',num2str(sprintf('%.2f',power_corrected(n).*1000)), ' mW'];
      % Text information on ARBccd of smallest box
      POWERSCAT=['Scattered power = ',num2str(sprintf('%.2s',ARBccd(1).*muFcal.*1000)),' {\mu}W'];
      BRDFVALUE=['BRDF = ',num2str(sprintf('%.2s',BRDFyint(n))),' 1/str'];
      SPATIALFREQ=['Spatial Frequency =  ',num2str(sprintf('%.2s',f_quad_sum(n))),' 1/nm '];
      PSD=['PSD =  ',num2str(PSD_fx_fy(n),'%10.2e'),' nm^3'];
      EXPOSURETIME = ['Exposure Time = ',num2str(exposure_time(n)),' s'];
      
      % Can only report sigma if we have second point already)
      if n>1
          SIGMA=['\sigma = ',num2str(sigma_cumulative(n-1),2), ' nm'];
      else
          SIGMA='\sigma = undefined';
      end
      
      % Add text information to plot
      text(65600,max(bincentervalue),[' \newline',ANGLEtheta_s,...
          '\newline',POWERINC,'\newline',POWERSCAT,'\newline',BRDFVALUE,...
          '\newline',SPATIALFREQ,'\newline',PSD,'\newline',SIGMA],...
          'VerticalAlignment','top','HorizontalAlignment','right','FontSize',12)
      
      %---------------------------------------------------------------------------------------------------
      % Plot the PSD and cumulative sigma versus spatial frequency
      %---------------------------------------------------------------------------------------------------
      ax4 = nexttile([2,2]);
      
      if n>1

        yyaxis left
        
        loglog(f_quad_sum(1:n),PSD_fx_fy(1:n),'LineStyle','none',...
               'Marker','square','MarkerSize',5,'MarkerEdgeColor',[0 0 0.8]);
        ylabel('Power spectral density [nm^4]');
        xlabel('Spatial frequency [1/nm]','FontSize',20,'FontName','Times New Roman');
        ax4.XGrid = 'off';
        ax4.XMinorGrid = 'off';
        ax4.XMinorTick = 'off';
        ax4.Box = 'on';
        ax4.YGrid = 'on';
        ax4.YMinorGrid = 'on';
        ax4.YMinorTick = 'off';
        ax4.YTickMode = 'auto';
        ax4.FontName = 'Times New Roman';
        ax4.FontSize = 18;
        ax4.YColor = [0 0 0.8];
        
        yyaxis right
        semilogx(f_quad_sum(1:n-1),sigma_cumulative(1:n-1),'LineStyle',...
                 'none','Marker','square','MarkerSize',5,...
                 'MarkerEdgeColor',[1 0.7 0]);
        ylabel('\sigma [nm]');
        ax4.Box = 'off';
        ax4.YGrid = 'on';
        ax4.YMinorGrid = 'off';
        ax4.YMinorTick = 'off';
        ax4.YTickMode = 'auto';
        ax4.FontName = 'Times New Roman';
        ax4.FontSize = 18;
        ax4.YColor = [1 0.7 0];

        legend({'PSD','\sigma'},'FontSize',18,'FontName','Times New Roman');
        title('Power spectrum and Surface (rms) roughness','FontSize',20,...
              'FontName','Times New Roman','Interpreter','none');
      else
          loglog(f_quad_sum,PSD_fx_fy,'s','DisplayName',folder.sample_name,...
              'Color',[0.5 0 0])
      end
      
      % Save 3-pane subplot to fig and png files
      if n==total_images
          saveas(gcf,[folder.analysisPath,'FIG',slash,num2str(theta_rot_stage(n)),'.fig']);
      end
      figure_size = get(gcf, 'position');
      set(gcf,'PaperPosition',figure_size/100);

      if print_images == 1
        saveas(gcf,[folder.analysisPath,'PNG',slash,num2str(pic_ID(n)),'.png']);
      elseif print_images == 0
        if n == total_images
          saveas(gcf,[folder.analysisPath,'PNG',slash,num2str(pic_ID(n)),'.png']);
        end
      end
      close(gcf);

%% FIGURE 2: CCD IMAGES
%---------------------------------------------------------------------------------------------------
% Create PNGs of just the CCD images
%---------------------------------------------------------------------------------------------------
      if print_images == 1

        f2 = figure('Visible','off');
        ax5= gca;

        % Plot image with proper scale
        imagesc(img_x_mm, img_y_mm, img.fit);

        hold on;
        % Adjusts Colormap limits [cMin,cMax] on output images, change Brightness
        ax5.CLim = [clim_Min, clim_Max];
        colormap('gray');
        boxcolors = {[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1]};
        linewidth = {.5, .1, .1, .1, .1, .1};
        linestyle = {'-',':',':',':',':',':'};
        
        for z = 1:numROIs
          plot(squeeze(x_ellipse(z,:)), squeeze(y_ellipse(z,:)),'LineStyle',linestyle{z},'color',boxcolors{z}...
               ,'LineWidth',linewidth{z});
        end

        text((200/pix_per_mm) - pixel_x_shift, (3886/pix_per_mm) - pixel_y_shift, [' \newline',EXPOSURETIME,'\newline',ANGLEtheta_s,...
            '\newline',POWERINC,'\newline',POWERSCAT,'\newline',BRDFVALUE,...
            '\newline',SPATIALFREQ,'\newline',PSD,'\newline',SIGMA],...
            'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10, 'Color', [1-eps 1 1])
        
        title([folder.sample_name,', image ',char(image_name)], 'FontSize',20,'FontName','Times New Roman','Interpreter', 'none');
        set(gcf,'units','points','position',[0 0 480 480]) % make it a square
        set(ax5,'position',[0 0 1 1]) % make square axes fill the figure
        axis square;
        axis off;

        print(gcf,[folder.analysisPath,'CCD_PNG',slash,num2str(n),'.png'], '-dpng','-r300');
        
        delete(f2);

        for z = 1:numROIs
          % Sets the RoI value in the structure to null, replacing the 4096x4096  
          RoI(z).fit = 0;
            
          % Sets the mask value in the structure to null, replacing the 4096x4096
          mask(z).fit = 0;
          
          % Overwrites the values of each element to null
          mask_crop(z).fit = 0;
          
          % Overwrites the values of each element to null
          RoI_mask(z).fit = 0;
        end

        % Overwrites the values of each element to null
        image.fit = 0;
        
        % Overwrites the values of each element to null
        image_dark.fit = 0;
      end
    end

%% Calculate TIS from BRDF data

    % Rayleigh-Rice TIS Approx for small scatter angles in (nm)
    [TIS_Rayleigh] = (((4 * pi * sigma_rel * cosd(correction_angle)) / lambda))^2;
    
    % Calculate TIS from integrating BRDF over Solid Angle
    % using equations 6,7 in Magana-Sandoval et al 2012

    for g = 1:length(theta_s)-1
      % Set start and stop polar angles
      theta_1 = theta_s(g);
      theta_2 = theta_s(g+1);

      % Calcute solid angle of the ring
      Omega_int_ring(g) = 2*pi*(1-cosd(theta_1));
      Omega_ext_ring(g) = 2*pi*(1-cosd(theta_2));
      Omega_ring(g) = Omega_ext_ring(g) - Omega_int_ring(g);
      
      % Get mean of mean of y-intercept BRDF [1/str]
      % for that angle range
      meanBRDFyint(g) = mean([BRDFyint(n), BRDFyint(g+1)]);
      
      % Get integral for this angle range
      % To get integrated scattering across hemisphere we multiply by the solid angle and use
      % cosine-corrected BRDF (multiply by cosine)
      intscat_a(g) = Omega_ring(g) * meanBRDFyint(g) * cosd(theta_s(g));
      
      theta_diff(g) = theta_2 - theta_1;
    end   

    Omega_tot = sum(Omega_ring);
    RTIS = sum(intscat_a);
    RTIS_cumulative = cumsum(intscat_a);
    
    TIS_diff = zeros(length(theta_diff)-1);
    TIS_deriv = length(theta_diff)-1;

    for k = 1:length(theta_diff)-1
      
      TIS_diff(k) = RTIS_cumulative(k+1) - RTIS_cumulative(k);
      TIS_deriv(k) = TIS_diff(k)./theta_diff(k);
    end

%% FIGURE 3: BRDF VS TTIME
%---------------------------------------------------------------------------------------------------
% Plot and save BRDF vs theta_s and BRDFlimit
%---------------------------------------------------------------------------------------------------
    
    if print_images == 1
      f3 = figure('Visible','off');
      ax3 = gca;
      
      semilogy(theta_s,minLIM,'v-','DisplayName','BRDF limit');
      hold on;
      grid on;
      
      errorbar(theta_s,BRDFyint,uncert_BRDFyint,'s','MarkerFaceColor','g','LineStyle','none','Color', 'g', 'DisplayName', sample);
      xlabel('\theta { [degs.]}');
  
      ylabel('BRDF [1/str]');
      BRDF_title = ['BRDF ',sample];
      title(BRDF_title);
  
  %---------------------------------------------------------------------------------------------------
  % Figure 3 Properties
  %---------------------------------------------------------------------------------------------------
      ax3.FontName = 'Times New Roman';
      ax3.FontSize = 30;
      f3.Position = [0,0,1618,1000];
  
      saveas(f3,[folder.analysisPath,folder.sample_name,slash,'BRDF_',folder.timeAnalysis,'_', folder.sample_name,'.fig']);
      saveas(f3,[folder.analysisPath,folder.sample_name,slash,'BRDF_',folder.timeAnalysis,'_', folder.sample_name,'.png']);
      delete(f3);
    end
    
%% FIGURE 4: POWER SPECTRAL DENSITY PLOT

% Plot and save PSD and Sigma Plots as .fig, .png
%---------------------------------------------------------------------------------------------------
    if print_images == 1
      f4 = figure('Visible','off');
      ax_PSD = gca;
  
      yyaxis left
      plot_PSD = loglog(f_quad_sum, PSD_fx_fy);
  
      ax_PSD.YLabel.String = 'Power Spectral Density [nm^4]';
      ax_PSD.XLabel.String = 'spatial frequency [1/nm]';
  
      ax_PSD.Box = 'on';
      ax_PSD.YGrid = 'on';
      ax_PSD.YMinorGrid = 'on';
      ax_PSD.YMinorTick = 'off';
      ax_PSD.XGrid = 'off';
      ax_PSD.XMinorGrid = 'on';
      ax_PSD.YTickMode = 'auto';
      ax_PSD.FontName = 'Times New Roman';
      ax_PSD.FontSize = 20;
      ax_PSD.YLim = [min(PSD_fx_fy) max(PSD_fx_fy)];
      ax_PSD.YColor = [0 0 0.8];    
  
      yyaxis right
      plot_sigma = loglog(f_quad_sum(1:end-1), sigma_cumulative);
  
      ax_PSD.YLabel.String = '\sigma [nm]';
  
      ax_PSD.YGrid = 'on';
      ax_PSD.YMinorGrid = 'off';
      ax_PSD.YMinorTick = 'off';
      ax_PSD.YTickMode = 'auto';
      ax_PSD.FontName = 'Times New Roman';
      ax_PSD.FontSize = 20;
      ax_PSD.XLim = [min(f_quad_sum) max(f_quad_sum)];
      ax_PSD.YColor = [1 0.7 0];
  
      plot_PSD.LineStyle = 'none';
      plot_PSD.Marker = 'square';
      plot_PSD.MarkerSize = 7;
      plot_PSD.MarkerFaceColor = [0 0 0.8];
  
      plot_sigma.LineStyle = 'none';
      plot_sigma.Marker = 'square';
      plot_sigma.MarkerSize = 7;
      plot_sigma.MarkerFaceColor = [1 0.7 0];
  
      legend([plot_PSD,plot_sigma],'PSD','\sigma');
      title('Power spectrum and surface (rms) roughness' ,'FontSize',30,'FontName','Times New Roman','Interpreter','none');
  
      saveas(f4, [directories.apath,directories.spath,'/','PSD_and_Sigma',directories.tpath,'.fig']);
      saveas(f4, [directories.apath,directories.spath,'/','PSD_and_Sigma',directories.tpath,'.png']);
      delete(f4);
    end
    
%% FIGURE 5: TIS PLOT
%---------------------------------------------------------------------------------------------------
    if print_images == 1
      TIS_title = strcat('Cumulative TIS ', sample);
      
      f5 = figure('Visible','off','Position',[1 1000 1000 1000]);
      ARS_ax = gca;
  
      set(ARS_ax,'YMinorTick','on','YMinorGrid','on',...
          'Position',[0.13 0.148923319327731 0.775 0.729067111294278],...
          'FontSize',30,...
          'FontName','Times New Roman');
      box(ARS_ax,'on');
      grid(ARS_ax,'on');
      hold(ARS_ax,'on');
      hold on;
      
      h = plot(theta_s(1:end-1), RTIS_cumulative);
      
      set(h,'LineStyle','none','Marker','square',...
          'MarkerSize',7,'MarkerFaceColor',[51 153 255]./255);
      
      ylabel('TIS','FontSize',30,'FontName','Times New Roman');
      xlabel('\theta_s {[deg.]}','FontSize',30,'FontName','Times New Roman');
      legend(TIS_title);
      title(TIS_title,'FontSize',36,'FontName','Times New Roman','Interpreter','none');
      grid('on');
      box('on');
      
      text(50,max(RTIS_cumulative),num2str(RTIS,3),'FontSize',20,...
          'Color',[51 153 255]./255);
      orient landscape;
      
      saveas(f5, [folder.analysisPath,folder.sample_name,slash,'TIS_',folder.timeAnalysis,'_', folder.sample_name,'.fig']);
      saveas(f5, [folder.analysisPath,folder.sample_name,slash,'TIS_',folder.timeAnalysis,'_', folder.sample_name,'.png']);
      delete(f5);
    end
      
%% FIGURE 6: TIS DERIVATIVE VS THETA PLOT

      % Plots value of TIS_derivative as a function of theta_s
    if print_images == 1
      f6 = figure('Position', [1 1000 1000 1000],'Visible','off');
      hold on;
      
      h = plot(theta_s(1:end-2), TIS_deriv,'Visible','off');
      
      set(h,'LineStyle','none','Marker','square',...
          'MarkerSize',7,'MarkerFaceColor',[51 255 249]./255);
      
      ylabel('Derivative dTIS/d\theta_s [ppm/deg.]','FontSize',30,'FontName','Times New Roman');
      xlabel('\theta_s {[deg.]}','FontSize',30,'FontName','Times New Roman');
      legend(h,'dTIS/d\theta_s');
      title('Derivative of TIS','FontSize',36,'FontName','Times New Roman','Interpreter','none');
      grid('on');
      box('on');
      orient landscape;
      
      saveas(f6, [folder.analysisPath,folder.sample_name,slash,'Derivative_TIS_',folder.timeAnalysis,'_', folder.sample_name,'.fig']);
      saveas(f6, [folder.analysisPath,folder.sample_name,slash,'Derivative_TIS_',folder.timeAnalysis,'_', folder.sample_name,'.png']);
      delete(f6);
    end

%% FIGURE 7: INCIDENT POWER

    if print_images == 1

      % Plots of incident power over scattering angle
      f7 = figure('Position',[0 0 1000 1000],'Visible','off');
      hold on;
      
      plot(theta_s, power_corrected,'s-','MarkerFaceColor','b');
      hold on;
      plot(theta_s, power_corrected,'s-','MarkerFaceColor','r');
      xlabel('\theta_s [deg]', 'FontSize', 20);
      title('Incident power over scattering angle','FontSize',30,'FontName','Times New Roman','Interpreter','none');
      
      ylabel('Power [W]', 'FontSize', 20);
      legend('Incident Power (measured twice, before and after)', 'Calibrated Monitor Power');
      grid('on');
      box('on');
      orient landscape;
      
      saveas(f7, [folder.analysisPath,folder.sample_name,slash,'Calibrated_monitor_power_',folder.timeAnalysis,'_', folder.sample_name,'.fig']);
      saveas(f7, [folder.analysisPath,folder.sample_name,slash,'Calibrated_monitor_power_',folder.timeAnalysis,'_', folder.sample_name,'.png']);
      delete(f7);
    end

%% FIGURE 8: LASER POWER PROFILE

    if print_images == 1
        % Laser Power Profile plots
        f8 = figure('Visible','off');
        yyaxis left
        plot(theta_s,BRDFyint,'b-')
        hold on
        ylabel('BRDF [1/str]','FontSize',18)
        set(gca,'FontSize',20)
        yyaxis right
        plot(theta_s,power_corrected,'r-')
        ylabel('Laser Power [W]','FontSize',18)
        set(gca,'FontSize',20)
        xlabel('Theta (deg)','FontSize',18)
    
        legend('BRDF [1/str]','Laser Power [W]','FontSize',14)
        title(['BRDF-LaserPower Profile ', strrep(folder.sample_name,'_','\_')],'FontSize',18)
        grid on
        grid minor
        saveas(f8, [folder.analysisPath,folder.sample_name,slash,'BRDF_LaserPower_', folder.sample_name,'_', folder.timeAnalysis, '.fig']);
        saveas(f8, [folder.analysisPath,folder.sample_name,slash,'BRDF_LaserPower_', folder.sample_name,'_', folder.timeAnalysis, '.png']);
        delete(f8);
    end
    
  case 'CRYO'
%% CRYO ANALYSIS
%-----------------------------------------------------------------------------------------------------------------------
% Prepare to open CCD images (get filenames, exposures, and RoIs)
%-----------------------------------------------------------------------------------------------------------------------

    % Logic to determine if the loop should only go to selected image for
    % locateROI or go through the whole loop for the analysis
    if locate_ROI==1
      endloop = image_selector;
    else
      endloop = 1:total_images;
    end

    for n=endloop

      s = s + 1;
      clc;
      % Image status and timing
      disp([num2str(s),' out of ', num2str(total_images)])
      disp(toc)
      tic

      % Appends our image data (pixel values) to the cell array
      image_name = readtable([imagedir,num2str(pic_ID(n)),'.csv']);
      
      % Grabs the exposure time saved in CamSetting.txt for the corresponding image
      
    
      % Read exposure times for the dark images. Determines if we are getting the
      % background image using the exposure time of each image
      if use_background_exp==1                                         
        % Rewrites filename for the case where we are selecting dark image files
        % based on their exposure time
        image_name = [num2str(exposure_time*100),'.csv'];
        darkfnm = readtable([darkimagedir,num2str(pic_ID(n)),'.csv']);
      else
        darkfnm = readtable([darkimagedir,num2str(pic_ID(n)),'.csv']);
      end
    
    
      % Creating ROIs
      for z = 1:numROIs
        
        scaledWidth(1,z)=ROIwidth+ROIscaleFactor*ROIwidth*((z-1)/(numROIs/2));
        scaledHeight(1,z)=ROIheight+ROIscaleFactor*ROIheight*((z-1)/(numROIs/2));
    
        % Determines the bottom used for each RoI
        ybottom(z) = round(ROIycenter - scaledHeight(1,z)/2);
    
        % Determines the top used for each RoI
        ytop(z) = round(ROIycenter + scaledHeight(1,z)/2);

        % determines the left value used for each RoI and each image
        xleft(z) = round(ROIxcenter - (scaledWidth(1,z)/2));

        % determines the right value used for each RoI and each image
        xright(z) = round(ROIxcenter + (scaledWidth(1,z)/2));
        xradius_ellipse(z) = (xright(z) - xleft(z))/2;

        % determines the xcenter used for each ellipse RoI and each image
        xcenter_ellipse(z) = (xright(z) + xleft(z))/2;

        % determines the yradius used for each ellipse RoI and each image
        yradius_ellipse(z) = (ytop(z) - ybottom(z))/2;

        % determines the ycenter used for each ellipse RoI and each image
        ycenter_ellipse(z) = (ybottom(z) + ytop(z))/2;   

        for i = 1:length(ellipse_theta)
          % determines the x value of the ellipse used for the RoI
          x_ellipse(z,i) = xradius_ellipse(z) * ...         
             cos(ellipse_theta(i)) + xcenter_ellipse(z);

          % determines the x value of the ellipse used for the RoI
          y_ellipse(z,i) = yradius_ellipse(z) * ...          
             sin(ellipse_theta(i)) + ycenter_ellipse(z);
        end
    
      end
   
    
    %% CRYO ROI MASK
    %---------------------------------------------------------------------------------------------------
    % This section makes the ROI mask for each image
    %---------------------------------------------------------------------------------------------------

      % Subtracted Image
      image = table2array(image_name - darkfnm);
      
      for z = 1:numROIs
        % This section of code creates an overlay of the RoI on top
        % of the images using a cropped "mask". This mask is made
        % with the poly2mask function. The function works as the
        % following: poly2mask(xi,yi,m,n) xi,yi are the verticies of
        % the closed region of the mask while m,n are the size of
        % the mask.

        % Create a mask to overlay with the image for analysis
        mask.fit = double(poly2mask(squeeze(x_ellipse(z,:)),...
                                    squeeze(y_ellipse(z,:)), 840, 840));

        % Create a cropped image of mask
        mask_crop.fit = mask.fit(ybottom(z):ytop(z),xleft(z):xright(z),:);

        % Create a cropped image of only the Region Of Interest
        RoI.fit = image(ybottom(z):ytop(z),xleft(z):xright(z));

        % Multiply RoI and mask_crop to make anything outside the ellipse 0
        RoI_mask.fit = RoI.fit .* mask_crop.fit;
        area(n,z) = nnz(RoI_mask.fit);

        % Calculate ARBccd summing over ROI pixel values and divide by exposure time
        ARBccd(z,n) = sum(sum(RoI_mask.fit))./(exposure_time);

        % Take matrix RoI(n).fit and turn into a giant column vector
        imagevector = reshape(RoI_mask.fit,numel(RoI_mask.fit),1);
      end
    
    %% CRYO BRDF CALCULATIONS
    %---------------------------------------------------------------------------------------------------
    % This section does the calculations of the data
    %---------------------------------------------------------------------------------------------------
      
      % This statement decides how background will be subtracted. From
      % background images or from defining a background ROI.
      
      for z = 1:numROIs
         % moved above for memory saving   
      end
      
      % Calculate BRDF,PSD,and Sigma for each image
      % BRDF calculation is based on Pmon_cal. Pmon_cal is the calibrated
      % monitor power calculated from power monitor measurements and adjusted
      % to an incident power using ratios. For TRS if powermonitor==0, ie
      % there is no power monitor data, Pmon_cal is just incidentPower, the initial
      % measured incident power
      for z = 1:numROIs
        % B [1/str] Transpose of ARBccd was taken so that it could be divided by 
        % incidentPower. Transmission of room light filter is accounted for here
        Bccd(z) = ARBccd(z,n)' .*muFcal./power_corrected(n)./filter(n);
        
        % BRDF [1/str]
        BRDFccd(n,z) = Bccd(z)./cosd(theta_s);
        
      end
      uncert_BRDFccd = BRDFccd*sys_uncert;
      
      % Fit the counts enclosed by each RoI vs the Area of each RoI with a
      % line y=mx+b where y is the counts enclosed, and x is the area
      % enclosed
    
      % Enclosed counts (actually BRDF associated with) that RoI
      y = BRDFccd(n,:);
    
      % Area enclosed by RoI
      x = area(n,:);
    
      % Uses MATLAB curve-fitting toolbox and creates a least-squares fit model to fit data values
      mdl = fitlm(area(n,:)',y');
    
      % Estimated y-intercept from fit
      b(n) = mdl.Coefficients{1,{'Estimate'}};
    
      % Standard error for y-intercept from fit
      SE(n) = mdl.Coefficients{1,{'SE'}};
      
      % Determines 95% confidence bounds for the y-intercept value (1.96 is the width of the std dev curve) [ref script: linear_fit_errors.m]
      ci(n) = SE(n).*1.96;
      
      % Set the BRDF value to the y-intercept of the linear fit
      if numROIs==1
          BRDFyint(n) = y;
      else
          BRDFyint(n) = b(n);
      end
      
      % Calculates total uncertainty in measurement by taking quadrature sum of
      % systematic and statistical errors (assumes errors are uncorrelated)
      tot_uncert(n) = sqrt(sys_uncert^2 + SE(n).^2);
      
      % Calculates the absolute uncertainty in BRDF measurement
      uncert_BRDFyint(n) = b(n)*tot_uncert(n);

      % Finds the x and y coordinates for all test box centers
      for N = 1:num_of_div
        LIMy(N)=N*(length(image(:,1))-100)/num_of_div;
        LIMx(N)=N*(length(image(1,:))-100)/num_of_div;
      end
      
      % Creates a counter to store BRDF information for each box
      LIMcount = 1;
    
      for Y = 1:num_of_div
        for X = 1:num_of_div % make a box
          LIMx1(LIMcount) = round(LIMx(X)-LIMwidth/2);
          LIMx2(LIMcount) = round(LIMx(X)+LIMwidth/2);
          LIMy1(LIMcount) = round(LIMy(Y)-LIMheight/2);
          LIMy2(LIMcount) = round(LIMy(Y)+LIMheight/2);
          % Crop beam image size
          LIMpic(LIMcount).fit = image(LIMy1(LIMcount):...
              LIMy2(LIMcount),LIMx1(LIMcount):LIMx2(LIMcount),:);
          % Sum over ROI and correct
          LIM_ARB(LIMcount) = sum(sum(LIMpic(LIMcount).fit))./(exposure_time);
          % B [1/sr] Transpose of ARBccd was taken so that it could be divided by
          % incidentPower. Transmission of room light filter is accounted for here
          LIM_B(LIMcount) = LIM_ARB(LIMcount)'.*muFcal./power_corrected(n)./filter(n);
          
          LIM_BRDF(LIMcount) = LIM_B(LIMcount)./cosd(theta_s) * n*((xright(1)-xleft(1))*(ytop(1)-ybottom(1))/(LIMwidth*LIMheight));
          LIMcount = LIMcount+1;
        end
      end
        
    
      %Find smallest LIM_BRDF and plot its box in image.
      [minLIM(n),minLIMpos] = min(LIM_BRDF);
      LIMITx1 = LIMx1(minLIMpos);
      LIMITx2 = LIMx2(minLIMpos);
      LIMITy1 = LIMy1(minLIMpos);
      LIMITy2 = LIMy2(minLIMpos);
      %END of BRDF limitations

%% CRYO (DEBUG) LOCATE ROI IMAGE PLOTTING

      if locate_ROI == 1
        % Display the corrected image to determine wehre the ROI min's and
        % max's are, and input them accordingly
    
        % Plot the image and RoIs
        ax1 = gca;
        imagesc(image);
        colorbar;
        % Keep the aspect ration [1 1 1] to avoid warping the image
        daspect([1 1 1]);
        
        hold on;
        ax1.CLim = [clim_Min, clim_Max];
        colormap('gray');
        boxcolors = {[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1]};
        
        % sized by the number of ROIs currently 6   
        linewidth = {.5, .1, .1, .1, .1, .1};
        linestyle = {'-',':',':',':',':',':'};
        
        for z = 1:numROIs
              plot(squeeze(x_ellipse(z,:)), squeeze(y_ellipse(z,:)),'LineStyle',linestyle{z},'color',boxcolors{z} ...
                   ,'LineWidth',linewidth{z});
        end

        text(200, 700, [' \newline',['Exposure Time = ',num2str(exposure_time),' s']...
                        ,' \newline',['BRDF = ',num2str(sprintf('%.2s',BRDFyint(n))),' 1/str']] ...
                        ,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10, 'Color', [1-eps 1 1]);
        return
      end

%% CRYO FIGURE 1: IMAGE, HISTOGRAM, BRDF, PSD
      if print_images == 1

        % Title of FIGURE 1
        cSiTitle = [sample,': Image ',num2str(n)];
      
        f1 = figure('Visible','off');
  
        % Plot the CCD image with its ROIs
        %---------------------------------------------------------------------------------------------------
        ax1a = nexttile([2,2]);
        imagesc(image);
        
        hold on;
        
        % Colors of each ROI ring
        boxcolors = {[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1]};
        
        % sized by the number of ROIs currently 6   
        linewidth = {.5, .1, .1, .1, .1, .1};
        linestyle = {'-',':',':',':',':',':'};
          

        for z = 1:numROIs
          plot(squeeze(x_ellipse(z,:)), squeeze(y_ellipse(z,:)),...
               'LineStyle',linestyle{z},'color',boxcolors{z},...
               'LineWidth',linewidth{z});
        end
  
        colormap('gray');
  
        ax1a.FontSize = 18;
        ax1a.FontName = 'Times New Roman';
        ax1a.CLim = [clim_Min, clim_Max];
        ax1a.DataAspectRatio = [1 1 1];
  
        axis off;
      
        % Plot histogram
        %---------------------------------------------------------------------------------------------------
        ax1b = nexttile([2,2]);
      
        % Make the bin center cvalues where 65600 is ???????????????????????????
        bincentervalue = 0:100:65600;
      
        % Make the histogram
        h = histogram(imagevector,0:100:65700);
      
        % Extract the counts from the histogram
        numinbin = h.Values;
      
        % make a bar plot of log10 of the number for each bin and gets handle
        bh = bar(bincentervalue,numinbin);
        
        hold on;
  
        % Some made up power law trend to be understood later
        histtheory = 3e8*bincentervalue.^(-2);
        plot(bincentervalue,histtheory,'k--')
        
        ax1b.XLim = [0 65600];
        ax1b.YLim = [min(bincentervalue) max(bincentervalue)];
        ax1b.FontSize = 18;
        ax1b.XScale = 'log';
        ax1b.YScale = 'log';
        ax1b.XLabel.String = 'Pixel value [counts]';
        ax1b.YLabel.String = '# of pixels';
        ax1b.Title.String = 'Pixel values in region of interest';
  
        bh.FaceColor = [1 0 0];
        bh.EdgeColor = [1 0 0];
        
        % Create text information to add to plot
        ELAPSE_TIME = ['Elapsed Time = ', char(image_duration_time(n))];
        % ELAPSE_TIME = ['Elapsed Time = ', num2str(durationElapsedate_time(n)/86400),' Days'];
        POWERINC = ['Incident power = ',num2str(sprintf('%.2f',power_corrected(n)*1000)), ' mW'];
        POWERSCAT = ['Scattered power = ',num2str(sprintf('%.2s',ARBccd(1).*muFcal.*1000)),' {\mu}W'];
        BRDFVALUE = ['BRDF = ',num2str(sprintf('%.2s',BRDFyint(n))),' 1/str'];
        TEMP = ['Temperature = ', num2str(round(image_temp(n),1)), '{\circ}K' ];
        
        % Add text information to plot
        text(65600, max(bincentervalue), [' \newline',ELAPSE_TIME,...
            '\newline',POWERINC,'\newline',POWERSCAT,'\newline',BRDFVALUE],...
            'VerticalAlignment','top','HorizontalAlignment','right','FontSize',12)
        
        % Plot BRDF vs Time and Temperature vs Time on 3-pane subplot
        %---------------------------------------------------------------------------------------------------
      
        colororder([0 0.75 0.75;1 0.5 0]);
      
        ax1c = nexttile([2,4]);
        yyaxis left
        
        semilogy(image_duration_time(1:n),BRDFyint,'LineStyle','-');
      
        hold on;
  
        ax1c.XLabel.String = 'Time Elapsed [Hours:Minutes:Seconds]';
        ax1c.YLabel.String = 'BRDF [1/str]';
  
        yyaxis right
  
        semilogy(image_duration_time,image_temp(1:n));
  
        ax1c.YLabel.String = 'Temperature [K]';
      
        grid on
        grid minor
            
        ax1c.Title.String = ['Normalized Scatter and Temperature Profile ', strrep(folder.sample_name,'_','\_')];
        
        f1.Position = [0 0 120 100];
  
        % Only save the four panel figure file if it's the last loop itteration
        if n == total_images
          saveas(f1,[folder.analysisPath,'FIG',slash,num2str(pic_ID(n)),'.fig']);
        end

        print(f1,[folder.analysisPath,'PNG',slash,num2str(pic_ID(n)),'.png'], '-dpng','-r300');

        delete(f1);
      end
%% CRYO FIGURE 2: CCD IMAGES

    % Create separate set of pngs with just the CCD image
    %---------------------------------------------------------------------------------------------------
      if print_images == 1
        % These PNGs will be used to make the video separately with FFMPEG
        f2 = figure('Visible','off');
        ax2= gca;
  
        imagesc(image);
  
        [imageHeight,imageWidth] = size(image);
      
        hold on;
        

        boxcolors = {[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1],[0 1 1]};
        linewidth = {.5, .1, .1, .1, .1, .1};
        linestyle = {'-',':',':',':',':',':'};
        
        for z = 1:numROIs
              plot(squeeze(x_ellipse(z,:)), squeeze(y_ellipse(z,:)),...
                   'LineStyle',linestyle{z},'color',boxcolors{z},...
                   'LineWidth',linewidth{z});
        end

        text((imageWidth-(0.05*imageWidth)), (0.05*imageHeight), ['\newline',ELAPSE_TIME,...
            '\newline',POWERINC,'\newline',POWERSCAT,'\newline',BRDFVALUE, '\newline', TEMP],...
            'HorizontalAlignment','right','VerticalAlignment','top','FontSize',11, 'Color', [1-eps 1 1])
  
        f2.Units = 'points';
        f2.Position = [0 0 imageWidth imageHeight];
  
        colormap('gray');

        ax2.DataAspectRatio = [1 1 1];
        ax2.CLim = [clim_Min, clim_Max];
        ax2.Position = [0 0 1 1];
  
        axis off;

        print(f2,[folder.analysisPath,'CCD_PNG',slash,num2str(pic_ID(n)),'.png'], '-dpng','-r300');

        delete(f2);
      end
      
      % Overwrites the values of each element to null
      img = zeros;
      image = zeros;
    end
    
%% CRYO FIGURE 3: 
    
    % BRDF vs theta_s and BRDFlimit
    %---------------------------------------------------------------------------------------------------

    BRDF_title = ['BRDF ',folder.sample_name];

    f3 = figure('Visible','off');

    ax3 = gca;

    semilogy(image_duration_time(1:total_images),minLIM,'v-','DisplayName','BRDF limit');

    box on
    grid on
    hold on;

    semilogy(image_duration_time(1:total_images),BRDFyint,'s','MarkerFaceColor','g','Marker','square', ...
             'LineStyle','none','Color', 'g', 'DisplayName', sample);

    legend('minLim', 'BRDFyint')

    f3.Position = [1 1000 1000 1000];

    ax3.Title.String = BRDF_title;
    ax3.XLabel.String = 'Elapsed Time (hh:mm:ss)';
    ax3.YLabel.String = 'BRDF [1/str]';
    ax3.YMinorTick = 'on';
    ax3.YMinorGrid = 'on';
    % ax3.Position = [0.13 0.148923319327731 0.775 0.729067111294278];

    saveas(gcf,[folder.analysisPath,folder.sample_name,slash,'BRDF_',...
           folder.timeAnalysis,'_', folder.sample_name,'.fig']);
    saveas(gcf,[folder.analysisPath,folder.sample_name,slash,'BRDF_',...
           folder.timeAnalysis,'_', folder.sample_name,'.png']);

    delete(f3);
    
%% CRYO FIGURE 4: GET RID OF?

    % Incident power over scattering angle
    %---------------------------------------------------------------------------------------------------
    
    f4 = figure('Visible','off');
    hold on;
    plot(image_duration_time, image_power,'s-','MarkerFaceColor','b');
    hold on;
    plot(image_duration_time, image_power,'s-','MarkerFaceColor','r');

    legend('Incident Power (measured twice, before and after)', 'Calibrated Monitor Power');
    
    grid on
    box on

    f4.Position = [1 1000 1000 1000];

    ax4.Title.String = 'Incident power over time';
    ax4.XLabel.String = 'Time (Days)';
    ax4.YLabel.String = 'Power [W]';
    
    saveas(f4, [folder.analysisPath,folder.sample_name,slash,'Calibrated_monitor_power_',folder.timeAnalysis,'_', folder.sample_name,'.fig']);
    saveas(f4, [folder.analysisPath,folder.sample_name,slash,'Calibrated_monitor_power_',folder.timeAnalysis,'_', folder.sample_name,'.png']);
    delete(f4);
    
%% CRYO FIGURE 5: TEMP POWER GRAPH
    % Temperature/Power Graph
    %---------------------------------------------------------------------------------------------------%%
    
    f5 = figure('Visible','off');
    ax5 = gca;

    yyaxis left

    plot(image_duration_time,image_temp,'b-')

    hold on
    
    ax5.YLabel.String = 'Temperature [{\circ}K]';
    
    yyaxis right

    plot(image_duration_time,image_transmitted,'r-')

    legend('Temperature', 'Transmitted Power')

    ax5.Title.String = ['Temperature-Power Profile ', strrep(folder.sample_name,'_','\_')];
    ax5.XLabel.String = 'Time (days)';
    ax5.YLabel.String = 'Transmitted Power';

    grid on
    grid minor

    saveas(f5, [folder.analysisPath,folder.sample_name,slash, 'Temp_Power_', folder.sample_name,'_', folder.timeAnalysis, '.fig']);
    saveas(f5, [folder.analysisPath,folder.sample_name,slash, 'Temp_Power_', folder.sample_name,'_', folder.timeAnalysis, '.png']);
    close(f5);
    
%% CRYO FIGURE 6: MONITORED POWER PLOT
    % Monitored power vs time plot that shows the entire power
    %-------------------------------------------------------------------------------------------------------------------
    
    if use_power_monitor==1
      figure('Visible','off');
      plot(image_duration_time, image_transmitted, 's-', 'MarkerFaceColor', 'r')
      xlabel('Elapsed Time (days)', 'FontSize', 20);
      ylabel('Power [W]', 'FontSize', 20);
      title('Monitored Power Over Time','FontSize',30,'FontName','Times New Roman','Interpreter','none');
      saveas(gcf, [folder.analysisPath,folder.sample_name,slash,'Monitor_power_',folder.timeAnalysis,'_', folder.sample_name,'.fig']);
      saveas(gcf, [folder.analysisPath,folder.sample_name,slash,'Monitor_power_',folder.timeAnalysis,'_', folder.sample_name,'.png']);
      delete(gcf);
    
    end
    
    figure('Visible','off');
    yyaxis left
    plot(image_duration_time(1:total_images),BRDFyint,'b-')
    hold on
    ylabel('BRDF [1/str]','FontSize',18)
    set(gca,'FontSize',12)
    
    yyaxis right
    plot(image_duration_time,image_power,'r-')
    ylabel('Laser Power [W]','FontSize',18)
    set(gca,'FontSize',12)
    hold;
    
    plot(image_duration_time,image_transmitted,'g-')
    set(gca,'FontSize',14)
    
    %legend('BRDF [1/str]','Laser Power [W]','FontSize',14)
    title(['BRDF-LaserPower Profile ', strrep(folder.sample_name,'_','\_')],'FontSize',18)
    grid on
    grid minor
    saveas(gcf, [folder.analysisPath,folder.sample_name,slash,'BRDF_LaserPower_', folder.sample_name,'_', folder.timeAnalysis, '.fig']);
    saveas(gcf, [folder.analysisPath,folder.sample_name,slash,'BRDF_LaserPower_', folder.sample_name,'_', folder.timeAnalysis, '.png']);
    delete(gcf);
end

%% DATA TABLE SAVING
%---------------------------------------------------------------------------------------------------
% Making Data table for analyzed BRDF Data vs Temperature
%---------------------------------------------------------------------------------------------------

switch experiment
  case 'ARS'
    % Custom MakeTable() function to make and save table. Enter cell of variable_names, cell of variables, and string of save path
    MakeTable({'Angle','BRDF','Temperature','TIS'},{theta_s,BRDFyint,image_temp,RTIS_cumulative},[folder.analysisPath,slash,'BRDF_Data.txt']);
    
  case {'CRYO','AAS'}
    % Custom MakeTable() function to make and save table. Enter cell of variable_names, cell of variables, and string of save path
    MakeTable({'Image_Time','Image_BRDF','Image_Temperature','Full_Time','Full_Chamber_Temp','Full_Stinger_Temp'},...
              {image_duration_time,BRDFyint,image_temp,duration_time,chamber_temp,stinger_temp},[folder.analysisPath,slash,'BRDF_Data.txt']);
end

%---------------------------------------------------------------------------------------------------
% Note of the type of Analysis
%---------------------------------------------------------------------------------------------------

file_name = 'Note.txt';
file_ID = fopen([folder.analysisPath,slash,file_name], 'w');
fprintf(file_ID, note_text);
fclose(file_ID);

%% Working on some things here. DO NOT DELETE! Let me COOK
%---------------------------------------------------------------------------------------------------

% figure; title('Relative Humidity Plot');
% hold on;

% p1 = plot(propper_time,humid.Sensor_1,'Marker','v');
% p2 = plot(propper_time,humid.Sensor_2,'Marker','v');

% xlabel('Time (hh:mm:ss)');
% ylabel('Relative Humidity (%)');
% legend('Door RH%','Camera RH%');
% p1.MarkerSize = 9; p2.MarkerSize = 9;
% p1.LineWidth = 4; p2.LineWidth = 4;
% fontsize(20,'points');
% pbaspect([1,1,1]);
% set(gcf,'Position',[0,0,1000,1000]); % Set position and size of plot that will be saved
% 
% saveas(gcf,'/Users/scatterlab/Desktop/PL_13649_AAS_Relative_Humidity.png');