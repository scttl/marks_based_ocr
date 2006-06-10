function imgs = add_img_noise(imgs, rp, ap, th, varargin)
% ADD_IMG_NOISE  Adds salt & pepper, and blurs the image(s) to make them noisy
%
%   imgs = add_img_noise(imgs, [rem_pct, add_pct, thresh, FSPECIAL_ARGS])
%
%   This function takes the binary image or images passed, turns on some
%   percentage of the 'off' pixels, turns off some percentage of the 'on' pixels
%   then passes a gaussian or other filter over the image to blur/sharpen it,
%   then rebinarizes it based on the threshold given.  The resulting noisy
%   binary image(s) are then returned.
%   
%   imgs should either be a matrix of binary values, or a cell array containing
%   matrices of binary values, each of which will have the same amount of noise
%   added.
%
%   rem_pct and add_pct are optional, and if specified should be between 0 and 1
%   give the percentage of on (respectively off) pixels to be switched.  This is
%   done randomly based on a uniform distribution.  If not passed, no pixels are
%   switched.
%
%   thresh is optional, and if specified it determines the resulting pixel 
%   values to be taken as a minimal to determine if the pixel should be 'on' 
%   once it is rebinarized.  Again it should be between 0 and 1, and defaults to
%   a value of 0.5.
%
%   the last set of optional parameters exactly follows the arguments to 
%   FSPECIAL, and allows one to select the size and type of filter to pass over
%   the image to blur/sharpen it.  For example, to create a spherical filter
%   5 pixels by 5 pixels, with a standard deviation of .5 you would tack:
%   ... 'gaussian', [3 3], .5)  on to the end of the argument list supplied
%   above.  See FSPECIAL for more information.  If not specified a 3x3 
%   'gaussian' filter with standard deviation of 0.5 is used.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: add_img_noise.m,v 1.1 2006-06-10 21:01:47 scottl Exp $
%
% REVISION HISTORY
% $Log: add_img_noise.m,v $
% Revision 1.1  2006-06-10 21:01:47  scottl
% Initial revision.
%


% LOCAL VARIABLES %
%%%%%%%%%%%%%%%%%%%

rem_pct = 0;
add_pct = 0;
thresh = 0.5;
fspec_args = {'gaussian', [3 3], 0.5};


% CODE START %
%%%%%%%%%%%%%%
if nargin < 1
    error('must specify the image to add noise to!');
elseif nargin >= 2
    rem_pct = rp;
    if nargin >= 3
        add_pct = ap;
        if nargin >= 4
            thresh = th;
            if nargin >= 5
                fspec_args = varargin;
            end
        end
    end
end

if ~ iscell(imgs)
    %single image to process
    imgs = process(imgs, rem_pct, add_pct, thresh, fspec_args);
else
    %multiple image cell array
    num_imgs = size(imgs,1);

    for i=1:num_imgs
        fprintf('  item %d      \r', i);
        imgs{i} = process(imgs{i}, rem_pct, add_pct, thresh, fspec_args);
    end
    fprintf('\n');
end

% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%process an individual image, adding noise to it
function img = process(img, rem_pct, add_pct, thresh, fspec_args)

%binarize the input if it isn't already binary
img = img > 0.5;
sz = size(img);

%add salt & pepper noise
all_idx = 1:(sz(1) * sz(2));
on_idx = find(img > 0);
len = length(on_idx);
switch_pos = randperm(len);
switch_pos = switch_pos(1:round(len * rem_pct));
img(on_idx(switch_pos)) = 0;

off_idx = setdiff(all_idx, on_idx);
len = length(off_idx);
switch_pos = randperm(len);
switch_pos = switch_pos(1:round(len * add_pct));
img(off_idx(switch_pos)) = 1;

%pass the filter over the image
H = fspecial(fspec_args{:});
img = imfilter(single(img), H);

%threshold the image to rebinarize and convert to logical array
img = img > thresh;
