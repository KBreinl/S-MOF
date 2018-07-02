function [interp_temp] = interpolation_temp(interp_temp,temp_D,temp_date_D,temp_date_H,id_rec,dist_temp,dist_calc)

% Loops through sites without hourly recording
for i=1:length(id_rec)
    for ii=1:length(temp_date_D);
          k=find(temp_date_H==temp_date_D(ii));
            % Build matrix with closest sites, rainfall amount and distance
            id=sum(interp_temp(k,dist_temp(id_rec(i),:)))';
            id(:,2)=dist_temp(id_rec(i),:);
            id(:,3)=dist_calc(id_rec(i),:);
            
            % If not enough neighbor sites duplicate the same site
            if size(id,1)<3
                id=repmat(id,3,1);
            end
            % Limit to three neighboring sites
            id=id(1:3,:);
            % Inverse distance weighing
            id(:,3)=1./id(:,3);
            id(:,3)=id(:,3)/sum(id(:,3));
            
            % Build raw circle from three sites (weighted by distance)
            circle_temp=interp_temp(k,id(1,2))*id(1,3)+interp_temp(k,id(2,2))*id(2,3)+interp_temp(k,id(3,2))*id(3,3);
            % build fragments
            circle_temp=circle_temp-mean(circle_temp);
            circle_temp=circle_temp+temp_D(ii,id_rec(i));
            interp_temp(k,id_rec(i))=circle_temp;
    end
end

