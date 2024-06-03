function initialize_OOEIT()

p = mfilename('fullpath');
p = p(1:end-length(mfilename));

addpath([p 'ForwardProblemSolvers\']);
addpath([p 'InverseSolvers\']);
addpath([p 'MiscClasses\']);
addpath([p 'MiscFunctions\']);
addpath([p 'Priors\']);

end