function [] = set_paths()
%--------------------------------------------------------------------------
%
% Copyright (c) 2014 Jeffrey Byrne <jebyrne@cis.upenn.edu>
%
%--------------------------------------------------------------------------

%% Disable name conflicts during
warning('off','MATLAB:dispatcher:nameConflict');  


%% Unpack Dependencies
deps = {'mdaisy-v1.0'};
for k=1:length(deps)
  % Unpack
  depdir = fullfile(pwd,'deps',deps{k});
  if ~exist(depdir,'dir')
    fprintf('[%s]: Unpacking %s\n', mfilename, depdir);
    zipfile = strcat(depdir,'.zip');
    unzip(zipfile, depdir);
  end

  % Paths
  fprintf('[%s]: Adding path %s\n', mfilename, depdir');
  addpath(genpath(depdir),'-begin');
end

%% Compile MEX

%% Restore
warning('on','MATLAB:dispatcher:nameConflict');  

