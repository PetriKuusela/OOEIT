function [sigma] = GenerateEllipseInGrid(gF,HF,min_sigma,max_sigma,width_x,width_y,...
cp_x,cp_y,band,CONDUCTIVE_BLOB)
%This function is used to generate targets for simulations

 gtmp = gF;
 x1a = linspace(cp_x-width_x/2,cp_x+width_x/2,500);
 x1b = x1a(2:end-1);
 y1a = cp_y-(width_y/2)*sqrt(1-(1-1e-5)*(x1a-cp_x).^2/(width_x/2)^2);
 y1b = cp_y+(width_y/2)*sqrt(1-(1-1e-5)*(x1b-cp_x).^2/(width_x/2)^2);
 x1 = [x1a(:);x1b(:)];
 y1 = [y1a(:);y1b(:)];
 tmp_pts1 = [x1,y1];
 tmp_pts1 = [tmp_pts1];

 x2a = linspace(cp_x-(width_x/2+band),cp_x+(width_x/2+band),500);
 x2b = x2a(2:end-1);
 y2a = cp_y-(width_y/2+band)*sqrt(1-(x2a-cp_x).^2/(width_x/2+band)^2+1e-12);
 y2b = cp_y+(width_y/2+band)*sqrt(1-(x2b-cp_x).^2/(width_x/2+band)^2+1e-12);
 x2 = [x2a(:);x2b(:)];
 y2 = [y2a(:);y2b(:)];
 tmp_pts2 = [x2,y2];

 rmind = find((1/(width_x/2+band)^2)*(gtmp(:,1)-cp_x).^2 ...
           + (1/(width_y/2+band)^2)*(gtmp(:,2)-cp_y).^2 < 1 & ...
           (1/(width_x/2)^2)*(gtmp(:,1)-cp_x).^2 ...
	      + (1/(width_y/2)^2)*(gtmp(:,2)-cp_y).^2 > 1);
 gtmp(rmind,:) = [];
 %lastmaxind = length(gtmp);
 gtmp = [gtmp;tmp_pts1;tmp_pts2];
 %minind = [lastmaxind+(1:length(tmp_pts1))];
 minind = find((1/(width_x/2+1e-2)^2)*(gtmp(:,1)-cp_x).^2 ...
	       + (1/(width_y/2+1e-2)^2)*(gtmp(:,2)-cp_y).^2 < 1);
 sigma_tmp = max_sigma*ones(size(gtmp,1),1);
 sigma_tmp(minind) = min_sigma;
 sigma = griddata(gtmp(:,1),gtmp(:,2),sigma_tmp,gF(:,1),gF(:,2));

 if CONDUCTIVE_BLOB
   sigma = min(min(sigma))+max(max(sigma)) - sigma;
 end

