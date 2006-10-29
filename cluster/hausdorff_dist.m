function d = hausdorff_dist(img1, other_imgs, pct)
%  HAUSDORFF_DIST  Calculate the Hausdorff distance from 1 image to others.
%
%    distance = HAUSDORFF_DIST(img1, other_imgs, [percent])
%
% img1 is an image intensity matrix, typically the average of
% all items in a cluster.
%
% other_imgs is either an individual image intensity
% matrix, or it is a cell array of image intensity matrices.
%
% percent is optional and if specified, determines the top percentage of pixel
% distances to be counted in the distance calculations (in ascending order).  
% If not specified, all pixels are included.  Lowering the value results in a 
% softer distance and can be useful for noisy documents.
%
% dist returned is a scalar (if other_imgs is a single image), or a vector (if
% other_imgs is a cell array), and is normalized by the area of the smaller of 
% the two images passed.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: hausdorff_dist.m,v 1.5 2006-10-29 17:23:17 scottl Exp $
%
% REVISION HISTORY
% $Log: hausdorff_dist.m,v $
% Revision 1.5  2006-10-29 17:23:17  scottl
% fix for handling blank images.
%
% Revision 1.4  2006/10/18 15:42:05  scottl
% implement ability to specify whether averaging between the two images will
% be used (or the maximal distance will be used)
%
% Revision 1.3  2006/08/14 01:38:06  scottl
% fix bug dealing with completely blank images.
%
% Revision 1.2  2006/07/21 21:37:39  scottl
% rewritten and vectorized (now compares multiple clusters at a time).  Also
% implemented soft matching.
%
% Revision 1.1  2006/06/03 20:55:48  scottl
% Initial check-in.


% LOCAL VARS %
%%%%%%%%%%%%%%
on_thresh = 0.5;  %when is a fractional pixel considered on?
percent = 1; %what portion of the closest pixel distances should be considered?

%should we average the "distances" between the images, or take the maximum?
use_avg = true;

% CODE START %
%%%%%%%%%%%%%%
if nargin < 2 || nargin > 3
    error('incorrect number of args!');
elseif nargin == 3
    percent = pct;
end

S1 = size(img1);

if iscell(other_imgs)
    num = length(other_imgs);
    if(size(other_imgs,1) ~= num)
        other_imgs = other_imgs';
    end
    S2 = zeros(num,2);
    for ii=1:num
        S2(ii,:) = size(other_imgs{ii});
    end
else
    num = 1;
    S2 = size(other_imgs);
end

%binarize img1
bin_img1 = img1 >= on_thresh;

%ensure there are at least some 'on' pixels in img1, otherwise there will be
%nothing to match!
img1_sum = sum(bin_img1(:));
if img1_sum == 0
    warning('MBOCR:blankImg', ...
            'trying to calculate the distance from a blank image!\n');
    d = inf(num,1);
    for ii=1:num
        if sum(other_imgs{ii}(:)) == 0
            d(ii) = 0;
        end
    end
    return;
end

%first get the maximum row and col size
max_row = max([S1(1); S2(:,1)]);
max_col = max([S1(2); S2(:,2)]);

%glue together the matrices in other_imgs with enough space between them to
%be able to calculate the Euclidean distance transform correctly.
other_mat = zeros(2*S1(1) + max_row, (2*S1(2) + max_col)*num, 'single');
img1_mat = logical(other_mat);
row_off = S1(1)+1;
start_col = S1(2)+1;
for ii = 1:num
    col_off = start_col;
    other_mat(row_off:row_off+S2(ii,1)-1, col_off:col_off+S2(ii,2)-1) = ...
              other_imgs{ii} >= on_thresh;
    start_col = start_col + 2*S1(2)+max_col;
end

%pass the first img over this matrix as a filter to get points of maximum
%overlap.  These will provide the centers for overlaying the img.
[row_scores, row_idx] = min(abs(filter2(img1, ...
                        other_mat(row_off:row_off+max_row-1,:))-img1_sum));
[col_idx, col_idx] = min(reshape(row_scores, 2*S1(2)+max_col, num));

%construct the glued together, maximally overlapping versions of the first image
half_h = ceil(S1(1)/2);
half_w = ceil(S1(2)/2);
if rem(S1(1),2)
    %odd height
    end_h = half_h - 2;
else
    %even height
    end_h = half_h - 1;
end
if rem(S1(2),2)
    %odd width
    end_w = half_w - 2;
else
    %even width
    end_w = half_w - 1;
end
start_col = 1;
for ii= 1:num
    tsc = start_col + col_idx(ii);
    tsr = row_off + row_idx(tsc);
    img1_mat(tsr-half_h:tsr+end_h,tsc-half_w:tsc+end_w) = bin_img1;
    start_col = start_col + 2*S1(2) + max_col;
end

%perform the Euclidean distance transform over the images
img1_dists = bwdist(img1_mat);
other_dists = bwdist(other_mat);

%read off the distances at the 'on' pixel locations in img1, centered about
%the maximally overlapping parts of the other img.  Sort them and read the
%smallest value that is within the appropriate percentile of pixels.  In this
%case, we can speed calculations up, since there are the same number of 
%points in each img1.
img1_to_other_d = sort(reshape(other_dists(img1_mat), img1_sum, num));
%calculate the soft distance by taking the largest value within threshold pixels
img1_to_other_d = img1_to_other_d(ceil(img1_sum * percent), :);

%calculate the distances at the 'on' pixel locations of the other images
%centered about the glued-together versions of img1.  Because the other images
%vary in size, this can't be calculated as efficiently as above.
other_to_img1_d = zeros(1,num);
other_sums = accumarray(reshape(repmat(1:num,2*S1(2)+max_col,1), num * ...
                        (2*S1(2) + max_col), 1), sum(other_mat, 'double'));
dist_vals = img1_dists(other_mat == 1);
start_idx = 1;
for ii=1:num
    if other_sums(ii) == 0
        %there are no 'on' pixels in the image, assume an infinite distance
        other_to_img1_d(ii) = inf;
    else
        this_d = sort(dist_vals(start_idx:start_idx+other_sums(ii)-1));
        other_to_img1_d(ii) = this_d(ceil(other_sums(ii) * percent));
    end
    start_idx = start_idx + other_sums(ii);
end

%now take the mean/maximum of the two distances to get the overall distance.
if use_avg
    d = mean([img1_to_other_d; other_to_img1_d])';
else
    d = max([img1_to_other_d; other_to_img1_d])';
end
