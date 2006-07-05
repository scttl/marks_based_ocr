function imgs = create_cluster_training_data(Comps, num)
% CREATE_CLUSTER_TRAINING_DATA  Create line images from component information
%
%   imgs = CREATE_CLUSTER_TRAINING_DATA(Comps, [num])
%
%   This function takes an already clustered (via cluster_comps) set of page
%   image data, and extracts individual lines as a cell array of images for
%   subsequent OCR processing.
%
%   Comps is a struct containing various fields including pg, position, and
%   associated neighbours of each connected component found while processing
%   page images.  This information will be used as the basis to extract lines
%   of an image to create each training case.  See cluster_comps for the format
%   of Comps.
%
%   imgs is the corresponding cell array containing image versions of the 
%   sentences in vals (stored as logical arrays, with a value of 1 representing
%   'on' pixels, and 0 'off').
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: create_cluster_training_data.m,v 1.3 2006-07-05 00:50:56 scottl Exp $
%
% REVISION HISTORY
% $Log: create_cluster_training_data.m,v $
% Revision 1.3  2006-07-05 00:50:56  scottl
% Changes based on new Cluster and Component structures.
%
% Revision 1.2  2006/06/19 21:50:03  scottl
% remove line equalization code, since we now handle lines of varying height
%
% Revision 1.1  2006/06/12 20:53:27  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%
max_cases = inf;

%these parameters control how noisy the images initially become (may be changed
%once the GUI is viewed).  See add_img_noise.  Also, setting rem_pct to a number
%less than 0, disables adding noise to the image
rem_pct = .01;
add_pct = .01;
thresh  = .5;
fspec_args = {'gaussian', [3 3], .5};


% CODE START %
%%%%%%%%%%%%%%
tic;
if nargin < 1 || nargin > 2
    error('incorrect number of arguments specified');
elseif nargin == 2
    max_cases = num;
end

%get line bounding boxes identifying each line
Pos = get_lines(Comps, max_cases);
fprintf('\n%.2fs: finished determining line boundaries\n', toc);

%now create images of each sentence
line = 1;
for pp = 1:size(Pos,1)
    M = ~imread(Comps.files{pp});  %flip pixel intensities
    for ll = 1:size(Pos{pp},1)
        fprintf('  item: %d\r', line);
        x = Pos{pp}(ll,:);
        imgs{line} = M(x(2):x(4), x(1):x(3));

        %trim leading and trailing background columns
        first_col = 1;
        while ~ any(imgs{line}(:,first_col))
            first_col = first_col + 1;
        end
        last_col = size(imgs{line},2);
        while ~ any(imgs{line}(:,last_col))
            last_col = last_col - 1;
        end
        imgs{line} = imgs{line}(:,first_col:last_col);

        line = line + 1;
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
