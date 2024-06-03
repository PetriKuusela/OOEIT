function [sigma] = GenerateEllipse(g,sigma_out,sigma_in,width_x,width_y,...
cp_x,cp_y,band)
%This function is used to generate targets for simulations

if width_x > width_y
    maxdist = width_x;
    fpdist = sqrt(width_x^2-width_y^2);
    fp1 = [cp_x-fpdist cp_y];
    fp2 = [cp_x+fpdist cp_y];
else
    maxdist = width_y;
    fpdist = sqrt(width_y^2-width_x^2);
    fp1 = [cp_x cp_y-fpdist];
    fp2 = [cp_x cp_y+fpdist];
end

dist1 = sqrt(sum((g - fp1).^2,2));
dist2 = sqrt(sum((g - fp2).^2,2));
dist = dist1 + dist2;

inell = dist < 2*maxdist;
outell = dist > 2*maxdist + 2*band;
inband = ~inell & ~outell;

sigma = zeros(size(g,1),1);
sigma(inell) = sigma_in;
sigma(outell) = sigma_out;
t = (dist(inband)-2*maxdist)/(2*band);
sigma(inband) = sigma_out*t + sigma_in*(1-t);

