function [] = demo_pvc(indir, outdir)
%------------------------------------------------------------------------
%
% Photorealistic virtual world (PVW) correspondences
%
% http://people.csail.mit.edu/biliana/projects/iccv2011/
%
% Jeffrey Byrne <jebyrne@cis.upenn.edu>
% 
%------------------------------------------------------------------------

camera = 1:4;  % camera index
location = {1:60, 1:61, 1:38, 1:41};  % location given camear
orientation = 1:3;  % camera orientation (yaw)
timeofday = 9:2:17;  % 9am - 5pm

% Read imagery from original PVC dataset 
im_obs = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', camera(1), location{1}(1), orientation(2), timeofday(1)))));
im_ref = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', camera(1), location{1}(2), orientation(2), timeofday(1)))));

% Load derived correspondences generated from pvc_correspondence
mat = load(fullfile(outdir, sprintf('asgn_%d_%dto%d_%d.mat', camera(1), location{1}(2), location{1}(1), orientation(2))));

% FILENAME FORMAT
% asgn_{CameraIndex}_{Location2Location}_{orientation}.mat
% asgn_{CameraIndex}_{LocationIndex}_{orientation2orientation}.mat
% asgn_{CameraIndex}_{Location2Location}_{orientation2orientation}.mat

% FRAME FORMAT:
% Pixel ij stored in mat.fr_obs(1:2,k) in image im_obs corresponds to pixel
% ij stored in mat.fr_obs2ref(1:2,k) in image im_ref.  

% Display correspondence scatterplot using the nsd toolbox
nsd.show.matching(im_ref, im_obs, mat.fr_obs(:,1:100:end)', mat.fr_obs2ref(:,1:100:end)', mat.fr_ref(:,1:100:end)');            

