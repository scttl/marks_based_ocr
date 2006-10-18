function d = conv_euc_dist(img1, other_imgs, norm_sq1, other_norms)
%  CONV_EUC_DIST  Calculate Euclidian distance between images after convolving
%
%    distance = CONV_EUC_DIST(img1, other_imgs, norm_sq1, other_norms)
%
% img1 is an image intensity matrices, typically the average of
% all items in a cluster.  other_imgs is either an individual image intensity
% matrix, or it is a cell array of image intensity matrices.
%
% norm_sq1 is a scalar representing the sqaured values of the L2-norm of img1.
% other_norms is either a scalar or a vector representing the squared L2-norms
% of the other_imgs.
%
% dist returned is a scalar (if other_imgs is a single image), or a vector (if
% other_imgs is a cell array), and is normalized by the area of the smaller of 
% the two images passed.

% CVS INFO %
%%%%%%%%%%%%
% $Id: conv_euc_dist.m,v 1.3 2006-10-18 15:45:38 scottl Exp $
%
% REVISION HISTORY
% $Log: conv_euc_dist.m,v $
% Revision 1.3  2006-10-18 15:45:38  scottl
% small efficiency fix.
%
% Revision 1.2  2006/07/22 04:09:38  scottl
% rewritten based on new Clust and Comps structures.
%
% Revision 1.1  2006/06/03 20:55:47  scottl
% Initial check-in.

% LOCAL VARS %
%%%%%%%%%%%%%%

% CODE START %
%%%%%%%%%%%%%%
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

eq_idx = (S1(1) == S2(:,1) & S1(2) == S2(:,2));
d = inf(num,1);
areas = min(repmat(prod(S1), num, 1), prod(S2,2));
num_eq = sum(eq_idx);

if num_eq > 0
    %no need to convolve these images
    d(eq_idx) = euc_dist(img1,other_imgs(eq_idx),norm_sq1,other_norms(eq_idx));
end
neq_idx = ~ eq_idx;
num_neq = sum(neq_idx);

if num_neq > 0
    max_row = max([S1(1); S2(neq_idx,1)]);
    max_col = max([S1(2); S2(neq_idx,2)]);
    
    %glue together the remaining matricies in other_imgs with enough space
    %between them to be able to calculate the maximal overlapping position
    %correctly.
    other_mat = zeros(2*S1(1) + max_row, (2*S1(2) + max_col)*num_neq);
    img1_mat = other_mat;
    row_off = S1(1)+1;
    start_col = S1(2)+1;
    for ii = find(neq_idx == 1)'
        col_off = start_col;
        other_mat(row_off:row_off+S2(ii,1)-1, col_off:col_off+S2(ii,2)-1) = ...
                  other_imgs{ii};
        start_col = start_col + 2*S1(2)+max_col;
    end
    
    %pass the first img over this matrix as a filter to get points of maximum
    %overlap.  These will provide the centers for overlaying img.
    img1_sum = sum(img1(:));
    [row_scores, row_idx] = min(abs(filter2(img1, ...
                            other_mat(row_off:row_off+max_row-1,:))-img1_sum));
    [col_idx, col_idx] = min(reshape(row_scores, 2*S1(2)+max_col, num_neq));

    %glue together copies of the first img on the maximally overlapping points
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
    for ii= 1:num_neq
        tsc = start_col + col_idx(ii);
        tsr = row_off + row_idx(tsc);
        img1_mat(tsr-half_h:tsr+end_h,tsc-half_w:tsc+end_w) = img1;
        start_col = start_col + 2*S1(2) + max_col;
    end

    %calculate the distance between these overlapping matrices
    d(neq_idx) = sqrt(norm_sq1 + other_norms(neq_idx) - 2 * accumarray(...
                 reshape(repmat(1:num_neq,2*S1(2)+max_col,1), num_neq * ...
                 (2*S1(2)+max_col), 1), sum(other_mat .* img1_mat))) ...
                 ./ areas(neq_idx);
end
