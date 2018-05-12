clc        % clear command window
clear      % clear all variables
%close all  % close all figures

%%
% The data file (AirQualityUCI.csv) is published 
%by UCI Repository and can be downloaded from: 
% https://archive.ics.uci.edu/ml/datasets/Air+quality

%%
% 
file_name = 'AirQualityUCI.mat';
load(file_name);
data = AirQualityUCI;

% get only benzen concentration
x= table2array(data(:,6)); % benzen concentration


% remove NaN values from the column
indx = ~isnan(x); % get indices of NaN values
x = x(indx);      % get actual values of ind

% Missing values are tagged with -200 value, so -200 is removed
x = x(x ~= -200);

% column 6 (benzene concentration is in microg/m^3)
% so to convert it to ppm (Part-per-million)
% 1 microg/m^3 = 0.001 ppm, to do so
x= x/1000;
