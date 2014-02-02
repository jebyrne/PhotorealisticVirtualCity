function [] = pvc_correspondence(indir, outdir)
%------------------------------------------------------------------------
%
% Photorealistic virtual world (PVW) correspondences
%
% http://people.csail.mit.edu/biliana/projects/iccv2011/
%
% Jeffrey Byrne <jebyrne@cis.upenn.edu>
% 
%------------------------------------------------------------------------

camera = 1:4;
location = {1:60, 1:61, 1:38, 1:41};
orientation = 1:3;
timeofday = 9:2:17;
%indir = '/Volumes/JEBYRNE-BACKUP/datasets/iccv2011_data';
%outdir = '/Volumes/JEBYRNE-BACKUP/datasets/pvw';
m = 480;  % image size
n = 640;
verbose = true;


%% Preprocessing
[V,U] = meshgrid(1:n, 1:m);
fr = [U(:) V(:) ones(480*640,1) zeros(480*640,1)]';  % dense pixel frame


%% Translation correspondence
for i=camera
  XYZ_obs = [];
  for j=location{i}   
    fprintf('[pvw][%d/%d][%d/%d]: translation correspondence\n', i, length(camera), j, length(location{i})); 
    X = csvread(fullfile(indir, sprintf('X_%d_%d_%d.txt', i, j, 2))); 
    Y = csvread(fullfile(indir, sprintf('Y_%d_%d_%d.txt', i, j, 2))); 
    Z = csvread(fullfile(indir, sprintf('Z_%d_%d_%d.txt', i, j, 2)));
    XYZ_ref = [X(:) Y(:) Z(:)]';

    if ~isempty(XYZ_obs)
      % Nearest neighbors in 3D world position
      k_validobs = find(~isnan(sum(XYZ_obs, 1)));
      k_validref = find(~isnan(sum(XYZ_ref, 1)));
      kdtree = vl_kdtreebuild(XYZ_ref(:, k_validref));  % current image is reference, (non NaN only)      
      [k_nn, dist] = vl_kdtreequery(kdtree, XYZ_ref(:,k_validref), XYZ_obs(:,k_validobs));
      k_matched = find(dist < 1);
            
      % Save me
      fr_obs = fr(:, k_validobs(k_matched));
      fr_ref = fr(:, k_validref);
      fr_obs2ref = fr(:, k_validref(k_nn(k_matched)));
      im_obs = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', i, j-1, 2, 9))));
      im_ref = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', i, j, 2, 9))));
      save(fullfile(outdir, sprintf('asgn_%d_%dto%d_%d.mat', i, j, j-1, 2)), 'fr_ref', 'fr_obs', 'fr_obs2ref', 'im_ref', 'im_obs');
      
      % Show me
      if verbose
        nsd.show.matching(im_ref, im_obs, fr_obs(:,1:100:end)', fr_obs2ref(:,1:100:end)', fr_ref(:,1:100:end)', figure(1));            
      end
    end
    XYZ_obs = XYZ_ref;    
  end
end


%% Orientation correspondence
XYZ_obs = [];
for i=camera
  for j=location{i}
    for k=[2 1 3]
      fprintf('[pvw][%d/%d][%d/%d]: orientation correspondence\n', i, length(camera), j, length(location{i}));
      X = csvread(fullfile(indir, sprintf('X_%d_%d_%d.txt', i, j, k)));
      Y = csvread(fullfile(indir, sprintf('Y_%d_%d_%d.txt', i, j, k)));
      Z = csvread(fullfile(indir, sprintf('Z_%d_%d_%d.txt', i, j, k)));
      
      % Nearest neighbors in 3D world position
      if k == 2
        XYZ_ref = [X(:) Y(:) Z(:)]';
      else
        XYZ_obs = [X(:) Y(:) Z(:)]';        
        
        % Nearest neighbors in 3D world position
        k_validobs = find(~isnan(sum(XYZ_obs, 1)));
        k_validref = find(~isnan(sum(XYZ_ref, 1)));
        kdtree = vl_kdtreebuild(XYZ_ref(:, k_validref));  % current image is reference, (non NaN only)
        [k_nn, dist] = vl_kdtreequery(kdtree, XYZ_ref(:,k_validref), XYZ_obs(:,k_validobs));
        k_matched = find(dist < 1);
        
        % Save me
        fr_obs = fr(:, k_validobs(k_matched));
        fr_ref = fr(:, k_validref);
        fr_obs2ref = fr(:, k_validref(k_nn(k_matched)));
        im_obs = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', i, j, k, 9))));
        im_ref = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', i, j, 2, 9))));
        save(fullfile(outdir,  sprintf('asgn_%d_%d_2to%d.mat', i, j, k)), 'fr_ref', 'fr_obs', 'fr_obs2ref', 'im_ref', 'im_obs');        
        
        % Show me
        if verbose
          nsd.show.matching(im_ref, im_obs, fr_obs(:,1:100:end)', fr_obs2ref(:,1:100:end)', fr_ref(:,1:100:end)', figure(1));
        end
      end      
    end
  end
end


%% Translation+orientation correspondence
for i=camera
  XYZ_obs = [];
  for j=location{i}(2:end)
    fprintf('[pvw][%d/%d][%d/%d]: translation+orientation correspondence\n', i, length(camera), j, length(location{i})); 
    X = csvread(fullfile(indir, sprintf('X_%d_%d_%d.txt', i, j, 1))); 
    Y = csvread(fullfile(indir, sprintf('Y_%d_%d_%d.txt', i, j, 1))); 
    Z = csvread(fullfile(indir, sprintf('Z_%d_%d_%d.txt', i, j, 1)));
    XYZ_ref = [X(:) Y(:) Z(:)]';

    X = csvread(fullfile(indir, sprintf('X_%d_%d_%d.txt', i, j-1, 3))); 
    Y = csvread(fullfile(indir, sprintf('Y_%d_%d_%d.txt', i, j-1, 3))); 
    Z = csvread(fullfile(indir, sprintf('Z_%d_%d_%d.txt', i, j-1, 3)));
    XYZ_obs = [X(:) Y(:) Z(:)]';
    
    % Nearest neighbors in 3D world position
    k_validobs = find(~isnan(sum(XYZ_obs, 1)));
    k_validref = find(~isnan(sum(XYZ_ref, 1)));
    kdtree = vl_kdtreebuild(XYZ_ref(:, k_validref));  % current image is reference, (non NaN only)
    [k_nn, dist] = vl_kdtreequery(kdtree, XYZ_ref(:,k_validref), XYZ_obs(:,k_validobs));
    k_matched = find(dist < 1);
    
    % Save me
    fr_obs = fr(:, k_validobs(k_matched));
    fr_ref = fr(:, k_validref);
    fr_obs2ref = fr(:, k_validref(k_nn(k_matched)));
    im_obs = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', i, j-1, 3, 9))));
    im_ref = rgb2gray(imread(fullfile(indir, sprintf('scene_%d_%d_%d_%d.png', i, j, 1, 9))));
    save(fullfile(outdir, sprintf('asgn_%d_%dto%d_1to3.mat', i, j, j-1)), 'fr_ref', 'fr_obs', 'fr_obs2ref', 'im_ref', 'im_obs');
    
    % Show me
    if verbose
      nsd.show.matching(im_ref, im_obs, fr_obs(:,1:100:end)', fr_obs2ref(:,1:100:end)', fr_ref(:,1:100:end)', figure(1));
    end
  end
end

