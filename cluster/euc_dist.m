function d = euc_dist(img1, other_imgs, norm_sq1, other_norms)
%  EUC_DIST  Calculate the normalized Euclidean distance from 1 image to others.
%
%    distance = EUC_DIST(img1, other_imgs, norm_sq1, other_norms)
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
    S2 = size(other_imgs);
end

eq_idx = (S1(1) == S2(:,1) & S1(2) == S2(:,2));
num_eq = sum(eq_idx);
if num_eq > 0
    d(eq_idx) = sqrt(norm_sq1 + other_norms(eq_idx) - 2 * accumarray(...
                reshape(repmat(1:num_eq,S1(1),1), num_eq * S1(1),1), ...
                sum(repmat(img1,num_eq,1) .* cell2mat(other_imgs(eq_idx)),2)));
end
neq_idx = ~ eq_idx;
num_neq = sum(neq_idx);
if num_neq > 0
    rows = min(repmat(S1(1),num_neq,1), S2(neq_idx,1));
    cols = min(repmat(S1(2),num_neq,1), S2(neq_idx,2));
    areas = min(repmat(prod(S1), num_neq,1), prod(S2(neq_idx),2));
    idx = find(neq_idx);
    diffs = zeros(num_neq,1);

    for ii=1:num_neq
        MM1 = img1(1:rows(ii), 1:cols(ii));
        MM2 = other_imgs{idx(ii)}(1:rows(ii), 1:cols(ii));
        diffs(ii) = sum(sum(MM1 .* MM2));
    end
    d(neq_idx) = sqrt(norm_sq1 + other_norms(neq_idx) - 2* diffs) ./ areas;
end
