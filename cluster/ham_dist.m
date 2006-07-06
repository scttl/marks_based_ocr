function d = ham_dist(img1, other_imgs)
%  HAM_DIST  Calculate the normalized Hamming distance from 1 image to others.
%
% img1 is an image intensity matrices, typically the average of
% all items in a cluster.  other_imgs is either an individual image intensity
% matrix, or it is a cell array of image intensity matrices.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: ham_dist.m,v 1.2 2006-07-06 17:53:14 scottl Exp $
%
% REVISION HISTORY
% $Log: ham_dist.m,v $
% Revision 1.2  2006-07-06 17:53:14  scottl
% rewritten to compute distances (after rename).  Vectorized to compute
% multiple distances in a single call.
%
% Revision 1.1  2006/06/03 20:55:48  scottl
% Initial check-in.
%

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

%glue together copies of the max. sized entry, take the hamming distance 
%over the entire thing, then determine the distance for each entry
max_rows = max([S1(1); S2(:,1)]);
max_cols = max([S1(2); S2(:,2)]);
areas = min(repmat(prod(S1), num,1), prod(S2,2));

img1mat = zeros(max_rows, max_cols);
img1mat(1:S1(1), 1:S1(2)) = img1;
img1mat = repmat(img1mat, num, 1);

othermat = zeros(num*max_rows, max_cols);
start_row = 1;
for ii=1:num
    othermat(start_row:start_row+S2(ii,1)-1, 1:S2(ii,2)) = other_imgs{ii};
    start_row = start_row + max_rows;
end
if num > 1
    d = accumarray(reshape(repmat(1:num,max_rows,1), num * max_rows,1), ...
        sum(xor(img1mat, othermat),2)) ./ areas;
else
    d = sum(sum(xor(img1mat, othermat))) ./ areas;
end
