function [] = eval_pvw()
%------------------------------------------------------------------------
%
% Photorealistic virtual world (PVW) correspondences
%
%------------------------------------------------------------------------
camera = 1:4;
location = {1:60, 1:61, 1:38, 1:41};
orientation = [3 1];
timeofday = 9:2:17;
indir = '/Volumes/JEBYRNE-BACKUP/datasets/iccv2011_data';
outdir = '/Volumes/JEBYRNE-BACKUP/datasets/pvw';
verbose = false;
desc = {'nsd', 'daisy', 'sift', 'orb', 'brisk', 'freak'};

addpath('/Users/jebyrne/software/mexopencv');
addpath(genpath('/Users/jebyrne/dev/pvw/deps'));


%% Translation evaluation
for i=camera
  for j=location{i}(2:end)
    for t=1:length(timeofday)
      fprintf('[pvw][%d/%d][%d/%d][%d/%d]: translation evaluation\n', i, length(camera), j, length(location{i}), t, length(timeofday));
      
      % Correspondences
      mat = load(fullfile(outdir, sprintf('asgn_%d_%dto%d_%d.mat', i, j, j-1, 2)));
      fr_obs = mat.fr_obs(:,1:100:end);
      fr_obs2ref = mat.fr_obs2ref(:,1:100:end);
      fr_ref = mat.fr_ref(:,1:100:end);
      
      % Images
      im_obs = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', i, j-1, 2, timeofday(t)))));
      im_ref = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', i, j, 2, timeofday(t)))));
      
      % Descriptors
      Z(i,j,t,:) = eval_descriptors(im_obs, im_ref, fr_obs, fr_obs2ref, fr_ref, desc, 10);
      fprintf('[pvw][%d/%d][%d/%d][%d/%d]: NSD=%f, DAISY=%f, SIFT=%f, ORB=%f, BRISK=%f, FREAK=%f\n', i, length(camera), j, length(location{i}), t, length(timeofday), Z(i,j,t,1), Z(i,j,t,2), Z(i,j,t,3), Z(i,j,t,4), Z(i,j,t,5), Z(i,j,t,6));
    end
  end
end
save(fullfile(outdir, sprintf('eval_translation.mat')), 'Z');



%% Orientation evaluation
Z = [];
for i=camera
  for j=location{i}
    for r=1:length(orientation)
      for t=1:length(timeofday)
        fprintf('[pvw][%d/%d][%d/%d][%d/%d][%d/%d]: orientation evaluation\n', i, length(camera), j, length(location{i}), r, length(orientation), t, length(timeofday));
        
        % Correspondences
        mat = load(fullfile(outdir, sprintf('asgn_%d_%d_2to%d.mat', i, j, orientation(r))));
        fr_obs = mat.fr_obs(:,1:100:end);
        fr_obs2ref = mat.fr_obs2ref(:,1:100:end);
        fr_ref = mat.fr_ref(:,1:100:end);
        
        % Images
        im_obs = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', i, j, orientation(r), timeofday(t)))));
        im_ref = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', i, j, 2, timeofday(t)))));
        
        % Descriptors
        Z(i,j,t,r,:) = eval_descriptors(im_obs, im_ref, fr_obs, fr_obs2ref, fr_ref, desc, 10);
        fprintf('[pvw][%d/%d][%d/%d][%d/%d]: NSD=%f, DAISY=%f, SIFT=%f, ORB=%f, BRISK=%f, FREAK=%f\n', i, length(camera), j, length(location{i}), t, length(timeofday), Z(i,j,t,r,1), Z(i,j,t,r,2), Z(i,j,t,r,3), Z(i,j,t,r,4), Z(i,j,t,r,5), Z(i,j,t,r,6));
      end
    end
  end
end
save(fullfile(outdir, sprintf('eval_orientation.mat')), 'Z');


%% Translation+orientation evaluation
Z = [];
for i=camera
  for j=location{i}(2:end)
    for t=1:length(timeofday)
      fprintf('[pvw][%d/%d][%d/%d][%d/%d]: translation+orientation evaluation\n', i, length(camera), j, length(location{i}), t, length(timeofday));
      
      % Correspondences
      mat = load(fullfile(outdir, sprintf('asgn_%d_%dto%d_1to3.mat', i, j, j-1)));
      fr_obs = mat.fr_obs(:,1:100:end);
      fr_obs2ref = mat.fr_obs2ref(:,1:100:end);
      fr_ref = mat.fr_ref(:,1:100:end);
      
      % Images
      im_obs = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', i, j-1, 3, timeofday(t)))));
      im_ref = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', i, j, 1, timeofday(t)))));
      
      % Descriptors
      Z(i,j,t,:) = eval_descriptors(im_obs, im_ref, fr_obs, fr_obs2ref, fr_ref, desc, 10);
      fprintf('[pvw][%d/%d][%d/%d][%d/%d]: NSD=%f, DAISY=%f, SIFT=%f, ORB=%f, BRISK=%f, FREAK=%f\n', i, length(camera), j, length(location{i}), t, length(timeofday), Z(i,j,t,1), Z(i,j,t,2), Z(i,j,t,3), Z(i,j,t,4), Z(i,j,t,5), Z(i,j,t,6));
    end
  end
end
save(fullfile(outdir, sprintf('eval_translation_and_orientation.mat')), 'Z');


%% Dense evaluation
Z = [];
stride = [2 4 8 16 32 64];
for i=camera
  for j=location{i}(2:end)   
    for t=1:length(timeofday)
      for s=1:length(stride)
        fprintf('[pvw][%d/%d][%d/%d][%d/%d]: dense translation evaluation\n', i, length(camera), j, length(location{i}), t, length(timeofday));
        
        % Images
        im_obs = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', i, j-1, 2, timeofday(t)))));
        im_ref = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', i, j, 2, timeofday(t)))));
        
        % Correspondences
        mat = load(fullfile(outdir, sprintf('asgn_%d_%dto%d_%d.mat', i, j, j-1, 2)));
        fr_obs = mat.fr_obs(:,1:100:end);
        fr_obs2ref = mat.fr_obs2ref(:,1:100:end);
        fr_ref = mat.fr_ref(:,1:100:end);
        
        % Strided correspondence
        fr_obs = stride(s)*round(fr_obs ./ stride(s));  fr_obs(3,:)=1; fr_obs(4,:)=0;
        fr_obs2ref = stride(s)*round(fr_obs2ref ./ stride(s));  fr_obs2ref(3,:)=1; fr_obs2ref(4,:)=0;
        fr_ref = stride(s)*round(fr_ref ./ stride(s));  fr_ref(3,:)=1; fr_ref(4,:)=0;
        k_validobs = nsd.util.inmat(size(im_obs), fr_obs(1,:), fr_obs(2,:));
        k_validobs2ref = nsd.util.inmat(size(im_obs), fr_obs2ref(1,:), fr_obs2ref(2,:));
        k_valid = intersect(k_validobs, k_validobs2ref);
        fr_obs = fr_obs(:,k_valid);
        fr_obs2ref = fr_obs2ref(:,k_valid);
        
        % Descriptors
        Z(i,j,t,s,:) = eval_descriptors(im_obs, im_ref, fr_obs, fr_obs2ref, fr_ref, desc, 10);
        fprintf('[pvw][%d/%d][%d/%d][%d/%d][%d/%d]: NSD=%f, DAISY=%f, SIFT=%f, ORB=%f, BRISK=%f, FREAK=%f\n', i, length(camera), j, length(location{i}), t, length(timeofday), s, length(stride), Z(i,j,t,s,1), Z(i,j,t,s,2), Z(i,j,t,s,3), Z(i,j,t,s,4), Z(i,j,t,s,5), Z(i,j,t,s,6));
      end
      %figure(100); plot(squeeze(Z(i,j,t,:,:)), 'MarkerSize', 10); legend(desc); grid on;  drawnow; 
    end
  end
end
save(fullfile(outdir, sprintf('eval_dense_translation.mat')), 'Z'); 



%% Plots
outdir = '/Volumes/JEBYRNE-BACKUP/datasets/pvw';

% Translation - mean per camera
mat = load(fullfile(outdir, sprintf('eval_translation.mat')));
for i=camera
  figure(10+i); hold on;  
  z = squeeze(mean(mat.Z(i,location{i}(2:end),:,:),3));  % remove zeros
  plot(z,'.-'); grid on; legend(upper(desc), 'Location','SouthEast');
  xlabel(sprintf('Position (Camera %d)', i));
  ylabel('Translation Matching Score');  
  ylim([0, 1]);
  export_fig(fullfile(outdir, sprintf('pvw_eval_translation_1_%d.png', i)), '-transparent');
end

% Translation - overall
mat = load(fullfile(outdir, sprintf('eval_translation.mat')));
for k=1:length(location)
  T(:,k,:,:) = mean(mat.Z(:,location{k},:,:),2);  % remove zeros
end
t = reshape(T, [size(T,1)*size(T,2)*size(T,3), size(T,4)]);
t = t(find(sum(t,2) > 0), :);
t = mean(t, 1);
figure(1); bar(t); set(gca,'xticklabel', upper(desc));
ylabel('Translation Matching Score');
ylim([0, 1]);
grid on;
export_fig(fullfile(outdir, 'pvw_eval_translation_2.png'), '-transparent');


% Matching Performance per byte
d_bytes = [512/8, (200*4), 128, 32, 64, 64];  % does not work, need another metric
d_bytes = exp(-(d_bytes-32).^2/128.^2);
figure(100); bar([t' (t.*d_bytes)']); set(gca,'xticklabel', upper(desc));
legend({'Matching','Storage Weighted Matching'}); grid on;
ylabel('Matching Score');
export_fig(fullfile(outdir, 'pvw_eval_swm.png'), '-transparent');


% Mean matching score vs. time 
t = squeeze(mean(reshape(T, [size(T,1)*size(T,2) size(T,3) size(T,4)]),1));
figure(2); plot(t,'.-'); grid on; legend(upper(desc), 'Location','SouthEast');
set(gca,'xticklabel', {'9am','10am','11am','12pm','1pm','2pm','3pm','4pm','5pm'});
ylabel('Translation Matching Score');
xlabel('Time of Day');
ylim([0, 1]);
grid on;
export_fig(fullfile(outdir, 'pvw_eval_translation_3.png'), '-transparent');

% Orientation
mat = load(fullfile(outdir, sprintf('eval_orientation.mat')));
for i=camera
  figure(20+i); hold on;  
  z = squeeze(mean(mean(mat.Z(i,location{i}(2:end),:,:,:),3),4));  % remove zeros
  plot(z,'.-'); grid on; legend(upper(desc), 'Location','SouthEast');
  xlabel(sprintf('Position (Camera %d)', i));
  ylabel('Rotation Matching Score');  
  ylim([0, 1]);
  export_fig(fullfile(outdir, sprintf('pvw_eval_rotation_1_%d.png', i)), '-transparent');
end


mat = load(fullfile(outdir, sprintf('eval_orientation.mat')));
for k=1:length(location)
  R(:,k,:,:,:) = mean(mat.Z(:,location{k},:,:,:),2);  % remove zeros
end
r = reshape(R, [size(R,1)*size(R,2)*size(R,3)*size(R,4), size(R,5)]);
r = r(find(sum(r,2) > 0), :);
r = mean(r, 1);
figure(3); bar(r); set(gca,'xticklabel', upper(desc));
ylabel('Rotation Matching Score');
grid on;
ylim([0,1]);
export_fig(fullfile(outdir, 'pvw_eval_rotation_2.png'), '-transparent');

% Mean matching score vs. time 
r = permute(R, [1 2 4 3 5]);
r = squeeze(mean(reshape(r, [size(r,1)*size(r,2)*size(r,3) size(r,4) size(r,5)]),1));
figure(4); plot(r,'.-'); grid on; legend(upper(desc),'Location','SouthEast');
set(gca,'xticklabel', {'9am','10am','11am','12pm','1pm','2pm','3pm','4pm','5pm'});
ylabel('Rotation Matching Score');
xlabel('Time of Day');
ylim([0, 1]);
grid on;
export_fig(fullfile(outdir, 'pvw_eval_rotation_3.png'), '-transparent');

% Translation+orientation
mat = load(fullfile(outdir, sprintf('eval_translation_and_orientation.mat')));
for i=camera
  figure(30+i); hold on;  
  z = squeeze(mean(mat.Z(i,location{i}(2:end),:,:),3));  % remove zeros
  plot(z,'.-'); grid on; legend(upper(desc), 'Location','SouthEast');
  xlabel(sprintf('Position (Camera %d)', i));
  ylabel('Translation+Rotation Matching Score');  
  ylim([0, 1]);
  export_fig(fullfile(outdir, sprintf('pvw_eval_translation_rotation_1_%d.png', i)), '-transparent');
end

% Translation - overall
mat = load(fullfile(outdir, sprintf('eval_translation_and_orientation.mat')));
T = [];
for k=1:length(location)
  T(:,k,:,:) = mean(mat.Z(:,location{k},:,:),2);  % remove zeros
end
t = reshape(T, [size(T,1)*size(T,2)*size(T,3), size(T,4)]);
t = t(find(sum(t,2) > 0), :);
t = mean(t, 1);
figure(5); bar(t); set(gca,'xticklabel', upper(desc));
ylabel('Translation+Rotation Matching Score');
ylim([0, 1]);
grid on;
export_fig(fullfile(outdir, 'pvw_eval_translation_rotation_2.png'), '-transparent');

% Mean matching score vs. time 
t = squeeze(mean(reshape(T, [size(T,1)*size(T,2) size(T,3) size(T,4)]),1));
figure(6); plot(t,'.-'); grid on; legend(upper(desc), 'Location','SouthEast');
set(gca,'xticklabel', {'9am','10am','11am','12pm','1pm','2pm','3pm','4pm','5pm'});
ylabel('Translation+Rotation Matching Score');
xlabel('Time of Day');
ylim([0, 1]);
grid on;
export_fig(fullfile(outdir, 'pvw_eval_translation_rotation_3.png'), '-transparent');


% Dense Translation - mean per stride
mat = load(fullfile(outdir, sprintf('eval_dense_translation.mat')));
for k=1:length(location)
  R(:,k,:,:,:) = mean(mat.Z(:,location{k},:,:,:),2);  % remove zeros
end
r = squeeze(mean(reshape(R, [size(R,1)*size(R,2)*size(R,3), size(R,4), size(R,5)]), 1));
figure(7); plot(r,'.-'); grid on; legend(upper(desc), 'Location','NorthEast');
xlabel(sprintf('Stride (log_2)'));
ylabel('Mean Matching Score');
ylim([0, 1]);
export_fig(fullfile(outdir, sprintf('pvw_eval_densestride.png', i)), '-transparent');




%% Helper functions
function Z = eval_descriptors(im_obs, im_ref, fr_obs, fr_obs2ref, fr_ref, desc, r_match)

opt = nsd.opts(); 
opt.descriptor.features.spyr.n_bands = 8;
opt.descriptor.features.spyr.n_scales = 8;
opt.descriptor.features.spyr.do_signed_orientation = true; 
verbose = false;

if isa(im_obs, 'uint16')
  im_obs = uint8(im_obs./255); % [pvw][1/4][11/60][1/2][1/5]
end
if isa(im_ref, 'uint16')
  im_ref = uint8(im_ref./255);
end
im_obs_padded = padarray(im_obs, [100 100], 'symmetric', 'both');
fr_obs_padded = fr_obs;  fr_obs_padded(1:2,:) = fr_obs_padded(1:2,:)+100;
im_ref_padded = padarray(im_ref, [100 100], 'symmetric', 'both');
fr_obs2ref_padded = fr_obs2ref;  fr_obs2ref_padded(1:2,:) = fr_obs2ref_padded(1:2,:)+100;

for k = 1:length(desc)
  switch desc{k}
    case 'nsd'
      opt.descriptor.mode = 'seedoflife';
      [d_obs, di_obs] = nsd.descriptor(im_obs_padded, opt.descriptor, fr_obs_padded);
      [d_ref, di_ref] = nsd.descriptor(im_ref_padded, opt.descriptor, fr_obs2ref_padded);
      
    case 'daisy'
      dzy_obs = compute_daisy(im_obs_padded);
      d_obs = dzy_obs.descs(sub2ind(size(im_obs_padded'), fr_obs_padded(2,:), fr_obs_padded(1,:)), :)';
      dzy_ref = compute_daisy(im_ref_padded);
      d_ref = dzy_ref.descs(sub2ind(size(im_ref_padded'), fr_obs2ref_padded(2,:), fr_obs2ref_padded(1,:)), :)';
      
    case 'sift'
      opt.descriptor.mode = 'sift';
      [d_obs, di_obs] = nsd.descriptor(single(im_obs_padded), opt.descriptor, fr_obs_padded);
      [d_ref, di_ref] = nsd.descriptor(single(im_ref_padded), opt.descriptor, fr_obs2ref_padded);
      
    case 'orb'
      extractor = cv.DescriptorExtractor('ORB');
      d_obs = to_binary(extractor.compute(im_obs_padded, frame2keypoint(fr_obs_padded))');
      d_ref = to_binary(extractor.compute(im_ref_padded, frame2keypoint(fr_obs2ref_padded))');            
      
    case 'brisk'
      extractor = cv.DescriptorExtractor('BRISK');
      d_obs = to_binary(extractor.compute(im_obs_padded, frame2keypoint(fr_obs_padded))');
      d_ref = to_binary(extractor.compute(im_ref_padded, frame2keypoint(fr_obs2ref_padded))');            

    case 'freak'
      extractor = cv.DescriptorExtractor('FREAK');
      d_obs = to_binary(extractor.compute(im_obs_padded, frame2keypoint(fr_obs_padded))');
      d_ref = to_binary(extractor.compute(im_ref_padded, frame2keypoint(fr_obs2ref_padded))');            
      
    case 'surf'
      extractor = cv.DescriptorExtractor('SURF');
      d_obs = extractor.compute(im_obs_padded, frame2keypoint(fr_obs_padded))';
      d_ref = extractor.compute(im_ref_padded, frame2keypoint(fr_obs2ref_padded))';            
      
    otherwise
      disp('unknown method');
  end
  
  % Descriptor distances
  D = sqdist(d_obs, d_ref);
  
  % Greedy bipartite matching
  [dists, perm] = sort(D(:),'ascend');
  numFramesA = size(D, 1);
  numFramesB = size(D, 2);
  [aIdx, bIdx] = ind2sub([numFramesA, numFramesB], perm(1:numel(dists)));
  edges = [aIdx bIdx];
  k_asgn = benchmarks.helpers.greedyBipartiteMatching(numFramesA, numFramesB, edges);  % seedoflife/deps/vlbenchmarks-1.0d
  
  % Matching
  k_match = find(sum((fr_obs2ref - fr_obs2ref(:,k_asgn)).^2,1) < r_match^2);
  y = zeros(length(k_asgn),1);
  y(k_match) = 1;
  Z(k) = mean(y);
  
  % Debugging
  if verbose
    nsd.show.matching(im_ref, im_obs, fr_obs', fr_obs2ref', fr_ref', figure(1));
    nsd.show.matching(im_ref, im_obs, fr_obs', fr_obs2ref(:,k_asgn)', fr_ref', figure(2));
  end  
end


%% Convert opencv byte packed binary to double binary
function B = to_binary(D)
n_bytes = size(D,1);  
n_desc = size(D,2);
B = [];
for k=1:n_bytes
  B = [B; (dec2bin(D(k,:), 8)-'0')'];
end

%% Convert nsd frame to opencv keypoint
function kp = frame2keypoint(fr)
n_desc = size(fr,2);
for i=1:n_desc
  kp(i).pt(1) = fr(2,i);
  kp(i).pt(2) = fr(1,i);  
  kp(i).size = 1;  
  kp(i).angle = -1;    
  kp(i).response = 0;      
  kp(i).octave = 0;        
  kp(i).class_id = -1;          
end
  
  

%% BUILD NOTES
% 
% >> mexopencv.make
% run matlab & from command line to fix tiff problem



