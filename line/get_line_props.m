function offs = get_baselines(imgs, thresh)
% get_baselines  Determine baseline offset in each image line passed.
%
%   offs = get_baselines(imgs, [threshold])
%
%   Given an individual or cell array of line images, this function attempts
%   to find the baseline (assuming the image contains a single line textual 
%   image).  The baseline is the bottom row at which most characters start
%   (exceptions are characters like j,g,q etc. which hang lower).
%
%   imgs should either be a single logical array representing the image, or it 
%   can be a cell array of logical images, each of which will be processed.
%
%   threshold is optional and if specified determines the minimum percentage
%   of pixels that must be 'on' for that row to be considered a baseline.  Since
%   there are fewer characters that hang below the baseline, this should be
%   possible to distinguish in typical lines.  The last such row found to exceed
%   this threshold is returned in the appropriate entry in the vector Offs.  If
%   not specified, threshold defaults to .20 (20% 'on' pixels)
%
%   If no rows meet the minimum threshold, the last row number is returned.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: get_line_props.m,v 1.2 2006-07-05 00:53:10 scottl Exp $
%
% REVISION HISTORY
% $Log: get_line_props.m,v $
% Revision 1.2  2006-07-05 00:53:10  scottl
% changed cluster and component structure, audited this file.
%
% Revision 1.1  2006/06/19 21:50:46  scottl
% initial revision.
%


% LOCAL VARIABLES %
%%%%%%%%%%%%%%%%%%%
base_thresh = .20;


% CODE START %
%%%%%%%%%%%%%%
if nargin < 1 || nargin > 2
    error('incorrect number of arguments specified!');
elseif nargin == 2
    base_thresh = thresh;
end

if iscell(imgs)
    num_imgs = length(imgs);
    for ii=1:num_imgs
        fprintf('processing image %d\n', ii);
        sz = size(imgs{ii});
        offs(ii) = sz(1);
        row_sums = sum(imgs{ii},2);
        idx = find((row_sums ./ sz(2)) >= base_thresh, 1, 'last');
        if idx
            offs(ii) = idx;
        end
    end
else
    %single image matrix
    fprintf('processing image 1\n');
    sz = size(imgs);
    offs = sz(1);
    row_sums = sum(imgs,2);
    idx = find((row_sums ./ sz(2)) >= base_thresh, 1, 'last');
    if idx
        offs = idx;
    end
end

