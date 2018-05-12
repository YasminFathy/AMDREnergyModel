clc        % clear command window
clear      % clear all variables
%close all  % close all figures

%%
% The data file (energydata_complete.csv) is published 
%by UCI Repository and can be downloaded from: 
% https://archive.ics.uci.edu/ml/datasets/Appliances+energy+prediction

%%
% 
file_name = 'energydata_complete.mat';
load(file_name);
data = energydatacomplete;
% get only T1 temperature in kitchen
x= table2array(data(:,4));