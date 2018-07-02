
clear
close all

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

% Read coordinates of sites
rain_D_LL=csvread('Daily_rain_LL.csv');
temp_D_LL=csvread('Daily_temp_LL.csv');
rain_H_LL=csvread('Hourly_rain_LL.csv');
temp_H_LL=csvread('Hourly_temp_LL.csv');

% Read in measuring interval for hourly rainfall and temperature
wind = read_param('wind');
wind = wind+1;
neighbor = read_param('nn');
historical=read_param('historical');
int_rain=read_param('int_rain');
int_temp=read_param('int_temp');

% Check and conduct interpolation (rainfall)
if size(rain_D_LL,1)~=size(rain_H_LL,1)
disp('Interpolation of rainfall...')
    dist_rain=zeros(rain_NOS_D,rain_NOS_H);
    dist_calc=zeros(rain_NOS_D,rain_NOS_H);
    for i=1:length(rain_D_LL)
        for ii=1:length(rain_H_LL)
            dist_calc(i,ii)=haversine(rain_D_LL(i,:),rain_H_LL(ii,:));
        end
        [sortedX, sortedInds] = sort(dist_calc(i,:),'ascend');
        dist_rain(i,:)=sortedInds;
        dist_calc(i,:)=sortedX;
    end
    
    % Chose closest sites (simple) or apply advanced interpolation
    if isequal(int_rain, 'simple')
        rain_H=rain_H(:,dist_rain(:,1));
        rain_NOS_H=rain_NOS_D;
    else
        rain_H=rain_H(:,dist_rain(:,1));
        rain_NOS_H=rain_NOS_D;
        id_rec=1:rain_NOS_D;
        id_rec=id_rec(dist_calc(:,1)~=0);
        rain_H(:,id_rec)=0;
        
        id_empty=1:rain_NOS_D;
        id_empty=id_empty(dist_calc(:,1)==0);
        for i=1:rain_NOS_D
            dist_rain(i,:)=id_empty(dist_rain(i,:));
        end
        
        rain_H=interpolation_rain(rain_H,rain_D,rain_date_D,rain_date_H,id_rec,dist_rain,dist_calc);
    end
end

% Check and conduct interpolation (temperature)
if size(temp_D_LL,1)~=size(temp_H_LL,1)
disp('Interpolation of temperature...')
    dist_temp=zeros(temp_NOS_D,temp_NOS_H);
    dist_calc=zeros(temp_NOS_D,temp_NOS_H);
    for i=1:length(temp_D_LL)
        for ii=1:length(temp_H_LL)
            dist_calc(i,ii)=haversine(temp_D_LL(i,:),temp_H_LL(ii,:));
        end
        [sortedX, sortedInds] = sort(dist_calc(i,:),'ascend');
        dist_temp(i,:)=sortedInds;
        dist_calc(i,:)=sortedX;
    end
    
    % Chose closest sites (simple) or apply advanced interpolation
    if isequal(int_temp, 'simple')
        temp_H=temp_H(:,dist_temp(:,1));
        temp_NOS_H=temp_NOS_D;
    else
        temp_H=temp_H(:,dist_temp(:,1));
        temp_NOS_H=temp_NOS_D;
        id_rec=1:temp_NOS_D;
        id_rec=id_rec(dist_calc(:,1)~=0);
        temp_H(:,id_rec)=0;
        
        id_empty=1:temp_NOS_D;
        id_empty=id_empty(dist_calc(:,1)==0);
        for i=1:temp_NOS_D
            dist_temp(i,:)=id_empty(dist_temp(i,:));
        end
        
        temp_H=interpolation_temp(temp_H,temp_D,temp_date_D,temp_date_H,id_rec,dist_temp,dist_calc);
    end
end

% Standardize daily data for distance calculation
rain_D_agg=rain_D;
for i=1:rain_NOS_D;
    k=find(rain_D_agg(:,i)>0);
    rain_D_agg(k,i)=sqrt(rain_D_agg(k,i));
end

temp_D_agg=(temp_D);

% Preallocate matrix for disaggregated simulations (w/o first and last 24h)
rain_H_sim=zeros(24,rain_NOS_D,length(rain_D)-2);
temp_H_sim=zeros(24,temp_NOS_D,length(temp_D)-2);

% Convert the hours to daily IDs
rain_date_H=datenum(year(rain_date_H),month(rain_date_H),day(rain_date_H));
temp_date_H=datenum(year(temp_date_H),month(temp_date_H),day(temp_date_H));

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
    rain_H_agg(k,i)=sqrt(rain_H_agg(k,i));
end

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
steps=length(rain_date_D_agg)-1;

for i=2:length(rain_date_D_agg)-1;
    % Define days for disaggregation within a window
    wind_da=rain_date_D_agg(i)-(wind+1):rain_date_D_agg(i)+(wind+1);
    v1=datevec(wind_da);
    v1(:,2:3) = 1;
    v2=datenum(wind_da') - datenum(v1) + 1;
    
    % Define day to disaggregate (with previous and following day and only
    % sites with rainfall (for comparison based on distances)
    rain_D_comp=rain_D_agg(i-1:i+1,:);
    temp_D_comp=temp_D_agg(i-1:i+1,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    % Reduce number of neighbors if the default cannot be fullfilled
    if isequal(historical, 'off')
        compl_test=(z(2:end-1,id_wet)>0)+0;
        if length(id_wet)>1
            compl_test=(sum((compl_test~=0)'))';
        end
        compl_test(:,2)=rain_date_D_window;
    elseif isequal(historical, 'on')
        compl_test=(z(2:end-1,id_wet)>0)+0;
        if length(id_wet)>1
            compl_test=(sum((compl_test~=0)'))';
        end
        compl_test(:,2)=year(rain_date_D_window);
        compl_test(compl_test(:,2)==year(rain_date_D(i)),1)=0;
        compl_test(:,2)=rain_date_D_window;
    end
    
    if sum(compl_test(:,1)==size(id_wet,2))<neighbor
        % Define correction of number if nearest neighbord (nn_c)
        nn_c=neighbor-sum(compl_test(:,1)==size(id_wet,2));
        % disp(nn_c)
    else
        nn_c=0;
    end
    
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
        temp_H_agg_comp(ii-1,:,1)=z(ii-1,:);
        temp_H_agg_comp(ii-1,:,2)=z(ii,:);
        temp_H_agg_comp(ii-1,:,3)=z(ii+1,:);
    end
    
    % Check if more than one station has rainfall or not
    if size(rain_H_agg_comp,2)>1;
        [a,b]=find(sum((rain_H_agg_comp(:,:,2)~=0)')==length(id_wet));
    elseif size(rain_H_agg_comp,2)==1;
        [b,a]=find(rain_H_agg_comp(:,:,2)~=0);
    end
    
    % Preallocate first distance matrix
    distance_matrix_1=zeros(length(z)-2,2)+9999;
    distance_matrix_1t=zeros(length(z)-2,2)+9999;
    
    % Add rainfall distances and type of distance measure
    dist='cityblock';
    
    % Distances for rainfall
    if isempty(id_wet)==0
        % Compare distances in regard to actual rain and temp day plus
        % precedent and consecutive rain day
        distance_matrix_1(b,1)=pdist2([rain_H_agg_comp(b,:,1),rain_H_agg_comp(b,:,2),rain_H_agg_comp(b,:,3)],[rain_D_comp(1,:),rain_D_comp(2,:),rain_D_comp(3,:)],dist);
    end
    
    % Distances for temperature
    distance_matrix_1t(:,1) = pdist2([temp_H_agg_comp(:,:,1),temp_H_agg_comp(:,:,2),temp_H_agg_comp(:,:,3)],[temp_D_comp(1,:),temp_D_comp(2,:),temp_D_comp(3,:)],dist);
    
    % Add corresponding day to distance matrix
    distance_matrix_1(:,2)=rain_date_D_window;
    distance_matrix_1t(:,2)=rain_date_D_window;
    
    distance_matrix_2t=distance_matrix_1t;
    distance_matrix_2t=sortrows(distance_matrix_2t,1);
    
    % Take out same year if resampling of historical values to avoid
    % recreation of observations
    if isempty(id_wet)==0
        distance_matrix_2=distance_matrix_1(distance_matrix_1(:,1)~=9999,:);
        distance_matrix_2=sortrows(distance_matrix_2,1);
        
        if isequal(historical, 'on')
            distance_matrix_2=distance_matrix_2(year(distance_matrix_2(:,2))~=year(rain_date_D_agg(i)),:);
        end
        
        % Nearest neighbor (NN) algorithm sampling
        distance_matrix_NN=zeros(neighbor-nn_c,2);
        distance_matrix_NN(:,1:2)=distance_matrix_2(1:neighbor-nn_c,:);
        distance_matrix_NN(distance_matrix_NN(:,1)==0,1)=0.000001;
        distance_matrix_NN(:,3)=cumsum((1:length(distance_matrix_NN(:,1)))'./distance_matrix_NN(:,1)/sum((1:neighbor-nn_c)'./distance_matrix_NN(:,1)));
        ran=rand;
        [~,I] = min(abs(distance_matrix_NN(:,3)-ran));
        k=find(rain_date_H==distance_matrix_NN(I,2));
        rain_H_disagg=rain_H(k,:);
        
        % Disaggregate the rain
        rain_disagg=bsxfun(@times,rain_D_disagg(2,:),bsxfun(@rdivide,rain_H_disagg(:,id_wet),sum(rain_H_disagg(:,id_wet))));
        rain_H_sim(:,id_wet,i-1)=rain_disagg;
    end
    
    if isequal(historical, 'on')
        distance_matrix_2t=distance_matrix_2t(year(distance_matrix_2t(:,2))~=year(rain_date_D_agg(i)),:);
    end
    
    distance_matrix_NNt(:,1:2)=distance_matrix_2t(1:neighbor,:);
    distance_matrix_NNt(distance_matrix_NNt(:,1)==0,1)=0.000001;
    distance_matrix_NNt(:,3)=cumsum((1:length(distance_matrix_NNt(:,1)))'./distance_matrix_NNt(:,1)/sum((1:neighbor)'./distance_matrix_NNt(:,1)));
    ran=rand;
    [~,I] = min(abs(distance_matrix_NNt(:,3)-ran));
    k=find(temp_date_H==distance_matrix_NNt(I,2));
    
    % Disaggregate the temperature
    temp_disagg=bsxfun(@plus,temp_D_disagg,bsxfun(@minus,temp_H(k,:),mean(temp_H(k,:))));
    temp_H_sim(:,:,i-1)=temp_disagg;
end;

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

% Smooth the temperature between the 24 hour blocks
for i=1:size(temp_H_sim,2);
    t=24;
    for ii=1:length(temp_H_sim)/24-2
        temp_H_sim(t-2:t+2,i)=smooth(temp_H_sim(t-2:t+2,i));
        t=t+24;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct time series 
disp('Reducing additional wet hours (correction routine)...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cut off too many rainfall values below length of observations and
% redistribute cut off rain first over values below threshold and the over
% all
for i=1:size(rain_H_sim,2)
    
    k_obs=find(rain_H(:,i)>0);
    obs=rain_H(rain_H(:,i)>0,i);
    k_sim=find(rain_H_sim(:,i)>0);
    sim=rain_H_sim(rain_H_sim(:,i)>0,i);
    
    dif=length(sim)-length(obs);
    
    if dif>0;
        k_sim(:,2)=rain_H_sim(rain_H_sim(:,i)>0,i);
        k_sim=sortrows(k_sim,2);
        
        data=(1:length(rain_H_sim))';
        data(ismember(data(:,1),k_sim(1:dif,1)),2)=data(ismember(data(:,1),k_sim(1:dif,1)),1);
        
        idx=ismember(rain_date_H,rain_date_H(k_sim(1:dif,1)));
        data2=rain_date_H(ismember(rain_date_H,unique(rain_date_H(idx))));
        data2(:,2)=data(idx,1);
        data2(:,3)=data(idx,2);
        data2(:,4)=rain_H_sim(idx,i);
        
        % Redistribute cut off rainfall to other amounts (more emphasis on
        % small values)
        rain_cut=sum(data2(data2(:,3)~=0,4));
        k=find(data2(:,3)==0 & data2(:,4)>0);
        data3=data2(k,4);
        
        % Put more emphasis on the small values
        data3=abs(data3-max(data3));
        
        %Redistribute
        data2(k,4)=data2(k,4)+data3/sum(data3)*rain_cut;
        rain_H_sim(k_sim(1:dif,1),i)=0;
        rain_H_sim(data2(k,2),i)=data2(k,4);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disaggregate the daily time series
disp('Write out disaggregated data...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rain_date_H, rain_H, rain_NOS_H] = read_file('Hourly_rain.csv');
[temp_date_H, temp_H, temp_NOS_H] = read_file('Hourly_temp.csv');

rain_date_H(1:24,:)=[];
rain_date_H((length(rain_date_H)-23):end,:)=[];

temp_date_H(1:24,:)=[];
temp_date_H((length(temp_date_H)-23):end,:)=[];

rain_H(1:24,:)=[];
rain_H((length(rain_H)-23):end,:)=[];

temp_H(1:24,:)=[];
temp_H((length(temp_H)-23):end,:)=[];

fclose('all');

disp('Simulation finished...')

clearvars -except rain_H_sim temp_H_sim rain_H temp_H 
