function nodeVals = GenerateEllipse(g,backgroundVal,inclusionVal,width_x,width_y,...
cp_x,cp_y,band)
%This function is used to generate targets for simulations. It creates a
%distribution with homogeneous background with value backgroundVal, and a
%single elliptic inclusion with value inclusionVal. The ellipse has width
%and height defined by width_x and width_y, respectively, and centerpoint
%defined by cp_x and cp_y. band defines how wide the transition zone from
%inclusion to the background is.

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

nodeVals = zeros(size(g,1),1);
nodeVals(inell) = inclusionVal;
nodeVals(outell) = backgroundVal;
t = (dist(inband)-2*maxdist)/(2*band);
nodeVals(inband) = backgroundVal*t + inclusionVal*(1-t);

