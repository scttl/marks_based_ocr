function imgs = create_cluster_training_data(Clust, Comps, num)
% CREATE_CLUSTER_TRAINING_DATA  Parse & create images clustered image files.
%
%   imgs = CREATE_CLUSTER_TRAINING_DATA(Clust, Comps, [num])
%
%   This function takes an already clustered (via cluster_comps) set of page
%   image data, and extracts individual lines as a cell array of images for
%   subsequent OCR processing.
%
%   Clust is an array of structs, each of which is assumed to contain several
%   fields of a particular format.  See cluster_comps for details.
%
%   Comps should be a cell array indexed by page, where each entry of each
%   entry represents a pixel, and each 'on' pixel is labelled with the component
%   number to which it belongs.
%
%   imgs is the corresponding cell array containing image versions of the 
%   sentences in vals (stored as logical arrays, with a value of 1 representing
%   'on' pixels, and 0 'off').
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: create_cluster_training_data.m,v 1.1 2006-06-12 20:53:27 scottl Exp $
%
% REVISION HISTORY
% $Log: create_cluster_training_data.m,v $
% Revision 1.1  2006-06-12 20:53:27  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%
max_cases = inf;

bg_val = 0;
fg_val = 1;

%these parameters control how noisy the images become.  See add_img_noise
%also, setting rem_pct to a number less than 0, disables adding noise to the 
%image
rem_pct = .01;
add_pct = .01;
thresh  = .5;
fspec_args = {'gaussian', [3 3], .5};


% CODE START %
%%%%%%%%%%%%%%
tic;
if nargin < 2 || nargin > 3
    error('incorrect number of arguments specified');
elseif nargin == 3
    max_cases = num;
end

%get line bounding boxes identifying each line
Pos = get_lines(Clust, Comps, max_cases);
fprintf('\n%.2fs: finished determining line boundaries\n', toc);

%now create images of each sentence
line = 1;
for p = 1:size(Pos,1)
    for l = 1:size(Pos{p},1)
        fprintf('  item: %d\r', l);
        x = Pos{p}(l,:);
        imgs{line,1} = zeros(x(4)-x(2)+1, x(3)-x(1)+1);
        imgs{line,1}(find(Comps{p}(x(2):x(4), x(1):x(3)) ~= bg_val)) = fg_val;
        line = line + 1;
    end
end
%we must ensure that all images are of the same height, so add zeros to the
%bottom of shorter rows, and chop off any leading rowspace
count = size(imgs,1);
first_row = inf;
for i = 1:count
    row = 1;
    while ~ any(imgs{i}(row,:))
        row = row + 1;
    end
    first_row = min(first_row, row);
end
for i = 1:count
    imgs{i} = imgs{i}(first_row:end,:);
end
max_height = 0;
for i = 1:count
    max_height = max(size(imgs{i},1), max_height);
end
for i = 1:count
    if size(imgs{i},1) ~= max_height
        imgs{i}(end+1:max_height,:) = 0;
    end
end
fprintf('\n%.2fs: finished creating images\n', toc);

%add noise to the images (using the GUI)
[rem_pct, add_pct, thresh, fspec_args] = noise_gui(imgs{1});
if rem_pct >= 0
    imgs = add_img_noise(imgs, rem_pct, add_pct, thresh, fspec_args{:});
    fprintf('%.2fs: finished adding noise to images\n', toc);
else
    fprintf('%.2fs: no noise being added\n', toc);
end



% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
