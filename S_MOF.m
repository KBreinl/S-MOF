
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in and prepare all data
disp('Read in and prepare all data...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in daily data
[rain_date_D, rain_D, rain_NOS_D] = read_file('Daily_rain.csv');
[temp_date_D, temp_D, temp_NOS_D] = read_file('Daily_temp.csv');

% Read in hourly data
[rain_date_H, rain_H, rain_NOS_H] = read_file('Hourly_rain.csv');
[temp_date_H, temp_H, temp_NOS_H] = read_file('Hourly_temp.csv');

% Date and matrices for the simulation (according to observations)
rain_date_D=datenum(rain_date_D);

% Reads in assignments between daily and hourly sites
links_rain=csvread('links_rain.csv');
links_temp=csvread('links_temp.csv');

% Extend the number of hourly sites to the number of daily sites
rain_H=rain_H(:,links_rain(:,2));
temp_H=temp_H(:,links_temp(:,2));

% Read in measuring interval for hourly rainfall and temperature
wind = read_param('wind');
wind = wind+1;
neighbor = read_param('nn');
historical=read_param('historical');

% Standardize daily data for distance calculation
rain_D_agg=rain_D;
for i=1:rain_NOS_D;
    k=find(rain_D_agg(:,i)>0);
    rain_D_agg(k,i)=zscore(rain_D_agg(k,i));
end

temp_D_agg=zscore(temp_D);

% Preallocate matrix for disaggregated simulations (w/o first and last 24h)
rain_H_sim=zeros(24,rain_NOS_D,length(rain_D)-2);
temp_H_sim=zeros(24,temp_NOS_D,length(temp_D)-2);

% Convert the hours to daily IDs
rain_date_H=datenum(year(rain_date_H),month(rain_date_H),day(rain_date_H));
temp_date_H=datenum(year(temp_date_H),month(temp_date_H),day(temp_date_H));

% Shift by one hour
rain_date_H=circshift(rain_date_H,1);
temp_date_H=circshift(temp_date_H,1);

rain_date_H(1)=[];
temp_date_H(1)=[];
rain_H(1,:)=[];
temp_H(1,:)=[];

% Preallocate matrix for nearest neighbor sampling
distance_matrix_NN=zeros(neighbor,3);

% Identify incomplete days for the rain with less than 24 hours
ident_r(1:length(rain_date_H),1)=9999;
for i=min(rain_date_H):max(rain_date_H);
    k=find(rain_date_H==i);
    if length(k)==24;
        ident_r(k,1)=1;
    else
        ident_r(k,1)=0;
    end
end

% Identify incomplete days for the temp with less than 24 hours
ident_t(1:length(temp_date_H),1)=9999;
for i=min(temp_date_H):max(temp_date_H);
    k=find(temp_date_H==i);
    if length(k)==24;
        ident_t(k,1)=1;
    else
        ident_t(k,1)=0;
    end
end

% Only keep days with 24 observation values per day
rain_date_H=rain_date_H(ident_r==1);
temp_date_H=temp_date_H(ident_t==1);
rain_H=rain_H(ident_r==1,:);
temp_H=temp_H(ident_t==1,:);

% Identify days with simultaneous rainfall and temp observations
k=intersect(rain_date_H,temp_date_H);
ident_r=ismember(rain_date_H,k);
ident_t=ismember(temp_date_H,k);
rain_date_H=rain_date_H(ident_r);
temp_date_H=temp_date_H(ident_t);
rain_H=rain_H(ident_r,:);
temp_H=temp_H(ident_t,:);

rain_H_agg(1:length(rain_H)/24,1:rain_NOS_D)=9999;
temp_H_agg(1:length(temp_H)/24,1:temp_NOS_D)=9999;

% Aggregate hourly values into daily for comparison to daily observations
[rain_date_D_agg,~,c] = unique(rain_date_H);
for i=1:rain_NOS_D
    agg_out = [rain_date_D_agg, accumarray(c,rain_H(:,i))];
    rain_H_agg(:,i)=agg_out(:,2);
end;

[rain_date_D_agg,~,c] = unique(temp_date_H);
for i=1:temp_NOS_D
    agg_out = [rain_date_D_agg, accumarray(c,temp_H(:,i))];
    temp_H_agg(:,i)=agg_out(:,2)/24;
end;

% Standardize data for distance calculation
for i=1:rain_NOS_D;
    k=find(rain_H_agg(:,i)>0);
    rain_H_agg(k,i)=zscore(rain_H_agg(k,i));
end

temp_H_agg=zscore(temp_H_agg);

% Convert day in day of the year (doy) for disaggregation
v0 = datevec(rain_date_D_agg);
v0(:,2:3) = 1;
rain_date_D_doy=datenum(rain_date_D_agg) - datenum(v0) + 1;

% Convert string date of daily simulations into numeric value
rain_date_D_agg=datenum(rain_date_D_agg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disaggregate the daily time series
disp('Disaggregate the daily time series...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through all simulated days and ignore first and last (no day
% before/after available)

h = waitbar(0,'Disaggregating daily data...');
steps=length(rain_date_D_agg)-1;

for i=2:length(rain_date_D_agg)-1;
    
    waitbar(i/steps)
    
    % Define days for disaggregation within a window
    wind_da=rain_date_D_agg(i)-(wind+1):rain_date_D_agg(i)+(wind+1);
    v1=datevec(wind_da);
    v1(:,2:3) = 1;
    v2=datenum(wind_da') - datenum(v1) + 1;
    
    % Define day to disaggregate (with previous and following day and only
    % sites with rainfall (for comparison based on distances)
    rain_D_comp=rain_D_agg(i-1:i+1,:);
    temp_D_comp=temp_D_agg(i,:);
    
    % Keep day to disaggregate non-standardized (for actual disaggregation)
    rain_D_disagg=rain_D(i-1:i+1,:);
    temp_D_disagg=temp_D(i,:);
    
    % Identify days of actual day to be disaggregated with rainfall
    id_wet=find(rain_D_comp(2,:)~=0);
    rain_D_comp=rain_D_comp(:,id_wet);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %Keep day to be disaggregated with original/non-standardized values
    rain_D_disagg=rain_D_disagg(:,id_wet);
    
    % Identify days for disaggregation in observed hourly (aggregated daily)
    % values (except first and last day)
    k=ismember(rain_date_D_doy,v2);
    rain_date_D_window=rain_date_D_agg(k);
    rain_date_D_window=rain_date_D_window(2:end-1);
    
    % rainfall matrix of neighbors to compare with
    z=rain_H_agg(k,:);
    
    % Preallocate
    rain_H_agg_comp=zeros(length(z)-2,length(id_wet));
    
    for ii=2:length(z)-1;
        rain_H_agg_comp(ii-1,:,1)=z(ii-1,id_wet);
        rain_H_agg_comp(ii-1,:,2)=z(ii,id_wet);
        rain_H_agg_comp(ii-1,:,3)=z(ii+1,id_wet);
    end
    
    % temperature matrix of neighbors to compare with
    z=temp_H_agg(k,:);
    
    % Preallocate temperature matrix and fill it
    temp_H_agg_comp=zeros(length(z)-2,temp_NOS_D);
    
    for ii=2:length(z)-1;
        temp_H_agg_comp(ii-1,:)=z(ii,:);
    end
    
    % Check if more than one station has rainfall or not
    if size(rain_H_agg_comp,2)>1;
        [a,b]=find(sum((rain_H_agg_comp(:,:,2)~=0)')==length(id_wet));
    elseif size(rain_H_agg_comp,2)==1;
        [b,a]=find(rain_H_agg_comp(:,:,2)~=0);
    end
    
    % Preallocate first distance matrix
    distance_matrix_1=zeros(length(z)-2,2)+9999;
    
    % Add rainfall distances and type of distance measure
    dist='cityblock';
    
    if isempty(id_wet)==0
        % Compare distances in regard to actual rain and temp day plus
        % precedent and consecutive rain day
        distance_matrix_1(b,1)=pdist2([rain_H_agg_comp(b,:,1),rain_H_agg_comp(b,:,2),rain_H_agg_comp(b,:,3),temp_H_agg_comp(b,:)],[rain_D_comp(1,:),rain_D_comp(2,:),rain_D_comp(3,:),temp_D_comp],dist);
    else
        % Add temperature distances if day is completely dry
        distance_matrix_1(:,1) = pdist2(temp_H_agg_comp,temp_D_comp,dist);
    end
    
    % Add corresponding day to distance matrix
    distance_matrix_1(:,2)=rain_date_D_window;
    
    if isempty(id_wet)==0
        % Filter out day with rainfall
        distance_matrix_2=distance_matrix_1(distance_matrix_1(:,1)~=9999,:);
        % Sort rows in distances according to minimum distance
        distance_matrix_2=sortrows(distance_matrix_2,1);
    elseif isempty(id_wet)==1
        distance_matrix_2=distance_matrix_1;
        distance_matrix_2=sortrows(distance_matrix_2,1);
    end
    
    % Take out same year if resampling of historical values to avoid
    % recreation of observations
    if isequal(historical, 'on')
        distance_matrix_2=distance_matrix_2(year(distance_matrix_2(:,2))~=year(rain_date_D_agg(i)),:);
    end;
    
    % Nearest neighbor (NN) algorithm sampling
    distance_matrix_NN(:,1:2)=distance_matrix_2(1:neighbor,:);
    distance_matrix_NN(distance_matrix_NN(:,1)==0,1)=0.000001;
    distance_matrix_NN(:,3)=cumsum((1:length(distance_matrix_NN(:,1)))'./distance_matrix_NN(:,1)/sum((1:neighbor)'./distance_matrix_NN(:,1)));
    ran=rand;
    [~,I] = min(abs(distance_matrix_NN(:,3)-ran));
    k=find(rain_date_H==distance_matrix_NN(I,2));
    
    % Disaggregate the rain
    rain_disagg=bsxfun(@times,rain_D_disagg(2,:),bsxfun(@rdivide,rain_H(k,id_wet),sum(rain_H(k,id_wet))));
    % Disaggregate the temperature
    temp_disagg=bsxfun(@plus,temp_D_disagg,bsxfun(@minus,temp_H(k,:),mean(temp_H(k,:))));
    
    if isempty(id_wet)==0
        rain_H_sim(:,id_wet,i-1)=rain_disagg;
    end
    
    temp_H_sim(:,:,i-1)=temp_disagg;
    
end

% Convert array of simulation into matrix
rain_H_sim=reshape(permute(rain_H_sim,[1,3,2]),[],rain_NOS_D);
temp_H_sim=reshape(permute(temp_H_sim,[1,3,2]),[],temp_NOS_D);

fclose('all');

% Remove first and last day as those are not disaggregated (to have the
% daily equivalent to the hourly data)
rain_D(1,:)=[];
rain_D=rain_D(1:end-1,:);
temp_D(1,:)=[];
temp_D=temp_D(1:end-1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disaggregate the daily time series
disp('Write out disaggregated data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% t = (datetime(year(rain_date_D(2)),month(rain_date_D(2)),day(rain_date_D(2)),0,0,0):...
%     hours:datetime(year(rain_date_D(end)),month(rain_date_D(end)),day(rain_date_D(end)),0,0,0))';
% 
% writetable(table(t(1:end-1),rain_H_sim), 'Hourly_rain_sim.csv', 'WriteVariableNames', false);
% writetable(table(t(1:end-1),sim_temp_h), 'Hourly_temp_sim.csv', 'WriteVariableNames', false);

clearvars -except rain_D temp_D rain_H_sim temp_H_sim 

