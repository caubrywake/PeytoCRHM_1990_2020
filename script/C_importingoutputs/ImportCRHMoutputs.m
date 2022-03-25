%% Import simulations results 

%% parameters
hruelev = [2949 2802 2701 2654 2560 2449 2252 2211 2176 2141 2956 2799 2660 2552 2460 2405 ...
2251 2705 2545 2445 2200 2741 2501 2554 2273 2215 2771 2502 2285 2702 2533 2375 ...
2393 2413 2482 2727 2847 ];
hruarea = [0.3231 1.786 0.2819 0.8719 1.459 0.8744 0.3325 0.2331 0.2971 0.0863 0.2594 1.104 0.9737 0.3394 0.2731 0.2081 ...
0.0719 0.5981 0.3425 0.2275 0.1306 0.965 0.7519 0.6163 0.4619 0.26 1.395 0.5512 0.4463 0.4806 0.4644 0.9269 ...
0.1144 0.1912 0.09313 0.35 0.2625 ];
save ('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\HRUelev_area.mat', 'hruelev', 'hruarea')

%% Setup
clear all
close all

addpath(genpath('D:\PeytoCRHM_1990_2020\chrm'))
addpath('D:\PeytoCRHM_1990_2020\function')
savedir = 'D:\PeytoCRHM_1990_2020\data_process\chrmoutput\'
figdir = 'D:\PeytoCRHM_1990_2020\fig\chrm\'

%% Load model result
files = dir('D:\PeytoCRHM_1990_2020\chrm\output\*.txt');   % load all the images in that directory
nfiles = length(files); 

for i = 1:nfiles
    filenames{:, i} = files(i,:).name
end 

filenames= filenames';

%% Load specific variable

[SWE] = ImportOutput(filenames{3}, 'SWE');
[timeCRHM] = ImportOutputTime(filenames{3});
figure; plot(timeCRHM, SWE)
timeCRHM = dateshift(timeCRHM,'start', 'hour', 'nearest');
save (strcat(savedir, 'SWE.mat'), 'SWE', 'timeCRHM')

[basinflow, basingw] = ImportOutput(filenames{4}, 'basinflow', 'basingw');
[timeCRHM] = ImportOutputTime(filenames{4});
timeCRHM = dateshift(timeCRHM,'start', 'hour', 'nearest');
save (strcat(savedir, 'basinflow.mat'), 'basinflow','basingw', 'timeCRHM')

[ice, firn] = ImportOutput(filenames{5}, 'ice', 'firn');
[timeCRHM] = ImportOutputTime(filenames{5});
timeCRHM = dateshift(timeCRHM,'start', 'hour', 'nearest');
save (strcat(savedir, 'firn_ice.mat'), 'firn', 'ice','timeCRHM')
figure; plot(timeCRHM, ice(:, [8,9,16,20]));
legend ('8','9','16','20')


[glacierh2o] = ImportOutput(filenames{3}, 'glacier_h2o');
[timeCRHM] = ImportOutputTime(filenames{3});
timeCRHM = dateshift(timeCRHM,'start', 'hour', 'nearest');
figure; plot(timeCRHM, glacierh2o(:, [8,9,16,20]));
legend ('8','9','16','20')
save (strcat(savedir, 'glacierh2o.mat'), 'glacierh2o', 'timeCRHM')

%% Import all files
close all
clear all
files = dir('D:\PeytoCRHM_1990_2020\chrm\output\*.txt');   % load all the images in that directory
savedir = 'D:\PeytoCRHM_1990_2020\data_process\chrmoutput\'

fileList = files; % get all the CRHM outputs starting with PeytoCUR_1OBS
 
for i = 1:numel(fileList)
fn = fileList(i).name; % file to import
H = importdata(fn,' ',2); %  % import headers 
D = importdata(fn) ; % import data
headers = regexp(H(1, 1), '\s+', 'split'); % split headers
headers = string(vertcat(headers{:})); % split headers
idxvar = [1, find(contains(headers,'(1)'))]; % select the number of variable in the file
% time is always the first one, followed by 2 or 3 variables
numvar = numel(idxvar); % number of variables

for ii= 1:numvar % for each variable, select all the column with data corresponding to that name
    varname =char(headers(idxvar(ii)));
    varname = strcat(varname(1:end-3)); % remove hru number from name
    if varname == 't' % excpetion is time - it does not have a number
        varname = 'time';
    end 
Index = strfind(headers, varname);
Index = find(not(cellfun('isempty', Index)));
varname = strcat(varname);
assignin('base',varname,D.data(:, Index));
end 

end 
%convert the time to a datetime format
time= datetime(datevec(time+ 693960));
time = dateshift(time,'start','hour', 'nearest');

save (strcat(savedir, 'output_20210323.mat'),'time','basinflow','basingw', ...
    'cumSWE_in','cumSWE_out', 'firn','firnmelt','hru_actet','hru_drift','hru_rain','hru_snow',...
    'hru_subl','hru_t','ice','icemelt','infil','meltrunoff','runoff','snowinfil','SWE','SWEmelt');