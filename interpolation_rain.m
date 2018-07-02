function [interp_rain] = interpolation_rain(interp_rain,rain_D,rain_date_D,rain_date_H,id_rec,dist_rain,dist_calc)

% Loops through sites without hourly recording
for i=1:length(id_rec)
    for ii=1:length(rain_date_D);
        if rain_D(ii,id_rec(i))>0
            k=find(rain_date_H==rain_date_D(ii));
            % Build matrix with closest sites, rainfall amount and distance
            id=sum(interp_rain(k,dist_rain(id_rec(i),:)))';
            id(:,2)=dist_rain(id_rec(i),:);
            id(:,3)=dist_calc(id_rec(i),:);
            
            % Remove sites without rainfall
            id(id(:,1)==0,:)=[];
            % If not enough neighbor sites duplicate the same site
            if size(id,1)<3
                id=repmat(id,3,1);
            end
            
            if isempty(id)==0
                % Limit to three neighboring sites
                id=id(1:3,:);
                % Inverse distance weighing
                id(:,3)=1./id(:,3);
                id(:,3)=id(:,3)/sum(id(:,3));
                
                % Build raw hyetograph from three sites (weighted by distance)
                hyeto=interp_rain(k,id(1,2))*id(1,3)+interp_rain(k,id(2,2))*id(2,3)+interp_rain(k,id(3,2))*id(3,3);
                % Calculate average duration (again weighted by distance)
                idx=round((sum(interp_rain(k,id(1,2))>0)*id(1,3)+sum(interp_rain(k,id(2,2))>0)*id(2,3)+sum(interp_rain(k,id(3,2))>0)*id(3,3)));
                % Cut out fraction (in sequence of overlapping hyetographs)
                if idx>1;
                    [a,~]=find(hyeto>0);
                    idx_range=randi([1 length(a)-idx+1],1);
                    idx_range=a(idx_range:idx_range+idx-1);
                    hyeto(ismember(1:24,idx_range)==0)=0;
                end
                % build fragments
                hyeto=hyeto./sum(hyeto);
                
                interp_rain(k,id_rec(i))=hyeto*rain_D(ii,id_rec(i));
                
            else
                % If there are not three stations with rainfall,take any
                % suitable day from the nearest site
                id2=rain_D(:,dist_rain(id_rec(i),1));
                id2(:,2)=datenum(rain_date_D);
                id2=sortrows(id2,1);
                id2(id2(:,1)==0,:)=[];
                id2(:,3)=cumsum((1:length(id2(:,1)))'./id2(:,1)/sum((1:length(id2))'./id2(:,1)));
                ran=rand;
                [~,idx] = min(abs(id2(:,3)-ran));
                idx=id2(idx,2);
                
                hyeto=interp_rain(datenum(rain_date_H)==idx,dist_rain(id_rec(i),1));
                hyeto=hyeto./sum(hyeto);

                interp_rain(k,id_rec(i))=hyeto*rain_D(ii,id_rec(i));
            end
        end
    end
end

