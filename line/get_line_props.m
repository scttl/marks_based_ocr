function [baselines, xheights] = get_line_props(imgs, varargin)
% GET_LINE_PROPS  Determine properties of line images (baseline, x-height)
%
%   [BASELINES, XHEIGHTS] = get_line_props(IMGS, [NAME1, VALUE1]...)
%
%   Given an individual or cell array of line images, this function attempts
%   to find the baseline (assuming the image contains a single line textual 
%   image), and x-height.  The baseline is the bottom row at which most 
%   characters start (exceptions are characters like j,g,q etc. which hang 
%   lower).  The x-height is the height of most lower-case letters (like x, a,
%   c, o, etc.), starting from the baseline.  Exceptions include capitals, and
%   large lower-case letters (like T, h, R, b, d, etc.)
%   NOTE: both these values are specified as the number of pixels offset from
%   the top region boundary of the line.
%
%   IMGS should either be a single logical array representing one line image, 
%   or it can be a cell array of logical images, each of which should
%   correspond to one line.  All will be processed.
%   NOTE: each image is assumed to have 0 represent background pixels, and
%   positive values representing foreground pixels.  Row sums are used to
%   calculate the properties.
%
%   The optional parameters allow one to override the default values defined
%   for the variables in LOCAL VARIABLES below.  These must be specified as
%   name, value pairs with the NAME being a string giving the name of the
%   variable to override, and the VALUE specifying its new value.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: get_line_props.m,v 1.3 2006-09-22 17:57:41 scottl Exp $
%
% REVISION HISTORY
% $Log: get_line_props.m,v $
% Revision 1.3  2006-09-22 17:57:41  scottl
% renamed from get_baselines, and moved to line directory.  Rewritten to
% include calculation of x-height of each line too.
%
% Revision 1.2  2006-07-05 00:53:10  scottl
% changed cluster and component structure, audited this file.
%
% Revision 1.1  2006/06/19 21:50:46  scottl
% initial revision.
%


% LOCAL VARIABLES %
%%%%%%%%%%%%%%%%%%%

%this parameter controls what portion of the pixels must be 'on' in a row to be
%considered a baseline.  The last such row found, is returned as a baseline for
%the image.
base_thresh = .20;

%similarly, this parameter controls what portion of the pixels must be 'on' in
%a row to be considered the x-height line.  The first such row found, is
%returned, with its value calculated as the number of pixels away from the
%baseline it resides (always positive)
xheight_thresh = .20;



% CODE START %
%%%%%%%%%%%%%%
if nargin < 1
    error('must specify a line image (at least one)!');
elseif nargin > 1
    process_optional_args(varargin{:});
end

if ~iscell(imgs)
    %single image passed
    imgs = {imgs};
end

num_imgs = length(imgs);

%preallocate the return parameters
baselines = NaN(num_imgs);
xheights = NaN(num_imgs);

for ii=1:num_imgs
    if num_imgs > 1 fprintf('processing image %d\n', ii); end
    sz = size(imgs{ii});
    baselines(ii) = sz(1) - 1; %use last row if all sums < thresh
    xheights(ii) = 0;  %use first row if all sums < thresh
    row_sums = sum(imgs{ii},2);
    idx = find((row_sums ./ sz(2)) >= base_thresh, 1, 'last');
    if idx
        baselines(ii) = idx - 1;
    end
    idx = find((row_sums ./ sz(2)) >= xheight_thresh, 1, 'first');
    if idx
        xheights(ii) = idx - 1;
    end
end
