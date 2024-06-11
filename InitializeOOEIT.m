function InitializeOOEIT()

%Find the path of the current folder, which is assumed to be the home
%folder for OOEIT
p = mfilename('fullpath');
%Remove the filename from the end of the path, so we get the folder path
p = p(1:end-length(mfilename));

%Add the relevant subfolders containing all the necessary OOEIT files
addpath([p 'ForwardProblemSolvers\']);
addpath([p 'InverseSolvers\']);
addpath([p 'MiscClasses\']);
addpath([p 'MiscFunctions\']);
addpath([p 'Priors\']);

end