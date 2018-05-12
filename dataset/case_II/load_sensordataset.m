%% 
clc, clear, close all


file_name = 'data.mat';
load(file_name);
% remove any rows with missing values in any column
data=data(~any(ismissing(data),2),:); 

sensors_dataset=table2array(data(:,4:end)); % columns in resultant matrix are
%moteid, temperature, humidity, light, voltage

% it seems the first row has invalid data, so just drop it
sensors_dataset = sensors_dataset(2:end,:);
% get unique mote ids
mote_ids = unique(sensors_dataset(:,1));

no_of_sensors = length(mote_ids); % number of sensors

% % create a matix for each mote
% for id=1:dataset_size
%    % get rows' indices of specific mote id
%    indx = find(sensors_dataset(:,1) == mote_ids(id));
%    % save all asociated info. for a given mote id into a matrix
%    sensor{id} = sensors_dataset(indx,:);    
% end
%no_of_sensors = length(sensor{1});% no. of sensor readings =43047
 
sensors = struct([]);% create a struct of sensors
i =1; % initialise an incremental variable

for id = 1:no_of_sensors % loop on the entire data set
   % get rows' indices of specific mote id
   indx = find(sensors_dataset(:,1) == mote_ids(id));
   if length(indx) >= 4000 
       
       %if all(sensors_dataset(indx,2) > 0) && all(sensors_dataset(indx,2) < 30)           
           indx = indx(1:4000);
           % save all asociated info. for a given mote id into the struct
           sensors(i).moteid = sensors_dataset(indx,1);
           sensors(i).temperature = sensors_dataset(indx,2);
           sensors(i).humidity = sensors_dataset(indx,3);
           sensors(i).light = sensors_dataset(indx,4);
           sensors(i).voltage = sensors_dataset(indx,5);
           no_of_sensors = i; % to make sure we have only the actual number of 
                              % selected sensors
           mote_ids_lst(i) = mote_ids(id);
       %end
       i = i + 1;
   end
end

% no of observations per sensor, all sensors have the same no of
% observations
no_of_observations = length(sensors(1).moteid);