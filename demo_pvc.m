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
do_plot = false;


%% Translation (camera 1)
% Read imagery from original PVC dataset
im_obs = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', camera(1), location{1}(1), orientation(2), timeofday(1)))));
im_ref = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', camera(1), location{1}(2), orientation(2), timeofday(1)))));

% Load derived correspondences generated from pvc_correspondence
mat = load(fullfile(outdir, sprintf('asgn_%d_%dto%d_%d.mat', camera(1), location{1}(2), location{1}(1), orientation(2))));

% Display correspondence scatterplot using the nsd toolbox
k = randperm(size(mat.fr_obs,2), 500);
nsd.show.matching(im_ref, im_obs, mat.fr_obs(:, k)', mat.fr_obs2ref(:, k)', []);
figure(1); subplot(1,2,1); title('Right'); subplot(1,2,2); title('Left');


%% Orientation (camera 2)
im_obs = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', camera(2), location{2}(1), orientation(3), timeofday(1)))));
im_ref = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', camera(2), location{2}(1), orientation(2), timeofday(1)))));

% Load derived correspondences generated from pvc_correspondence
mat = load(fullfile(outdir, sprintf('asgn_%d_%d_%dto%d.mat', camera(2), location{2}(1), orientation(2), orientation(3))));

k = randperm(size(mat.fr_obs,2), 500);
nsd.show.matching(im_ref, im_obs, mat.fr_obs(:, k)', mat.fr_obs2ref(:, k)', []);
figure(2); subplot(1,2,1); title('Right'); subplot(1,2,2); title('Left');


%% Translation + Orientation (camera 4)
im_obs = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', camera(4), location{4}(1), orientation(3), timeofday(1)))));
im_ref = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', camera(4), location{4}(2), orientation(1), timeofday(1)))));

% Load derived correspondences generated from pvc_correspondence
mat = load(fullfile(outdir, sprintf('asgn_%d_%dto%d_%dto%d.mat', camera(4), 2, 1, 1, 3)));

k = randperm(size(mat.fr_obs,2), 500);
nsd.show.matching(im_ref, im_obs, mat.fr_obs(:, k)', mat.fr_obs2ref(:, k)', []);
figure(3); subplot(1,2,1); title('Right'); subplot(1,2,2); title('Left');


%% Translation + Orientation (camera three)
im_obs = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', camera(3), location{4}(1), orientation(3), timeofday(1)))));
im_ref = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', camera(3), location{4}(2), orientation(1), timeofday(1)))));

% Load derived correspondences generated from pvc_correspondence
mat = load(fullfile(outdir, sprintf('asgn_%d_%dto%d_%dto%d.mat', camera(3), 2, 1, 1, 3)));

k = randperm(size(mat.fr_obs,2), 500);
nsd.show.matching(im_ref, im_obs, mat.fr_obs(:, k)', mat.fr_obs2ref(:, k)', []);
figure(4); subplot(1,2,1); title('Right'); subplot(1,2,2); title('Left');


%% Dense stride

% Read imagery from original PVC dataset 
im_obs = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', camera(1), location{1}(4), orientation(2), timeofday(1)))));
im_ref = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', camera(1), location{1}(5), orientation(2), timeofday(1)))));

% Load derived correspondences generated from pvc_correspondence
mat = load(fullfile(outdir, sprintf('asgn_%d_%dto%d_%d.mat', camera(1), location{1}(5), location{1}(4), orientation(2))));

% Strided correspondence
stride = [32 64];
for s=1:length(stride)
  fr_obs = mat.fr_obs(:,1:100:end);
  fr_obs2ref = mat.fr_obs2ref(:,1:100:end);
  fr_ref = mat.fr_ref(:,1:100:end);

  fr_obs = stride(s)*round(fr_obs ./ stride(s));  fr_obs(3,:)=1; fr_obs(4,:)=0;
  fr_obs2ref = stride(s)*round(fr_obs2ref ./ stride(s));  fr_obs2ref(3,:)=1; fr_obs2ref(4,:)=0;
  fr_ref = stride(s)*round(fr_ref ./ stride(s));  fr_ref(3,:)=1; fr_ref(4,:)=0;
  k_validobs = nsd.util.inmat(size(im_obs), fr_obs(1,:), fr_obs(2,:));
  k_validobs2ref = nsd.util.inmat(size(im_obs), fr_obs2ref(1,:), fr_obs2ref(2,:));
  k_valid = intersect(k_validobs, k_validobs2ref);

  fr_obs = fr_obs(:,k_valid);
  fr_obs2ref = fr_obs2ref(:,k_valid);  
  
  % Display correspondence scatterplot using the nsd toolbox  
  nsd.show.matching(im_ref, im_obs, fr_obs(:, :)', fr_obs2ref(:, :)', []);
  figure(4+s);  subplot(1,2,1); title('Right'); subplot(1,2,2); title('Left');
end


%% Plots
if do_plot
  for k=1:6
    export_fig(sprintf('pvc_demo_plot_%d.png',k) '-transparent');
  end
end

