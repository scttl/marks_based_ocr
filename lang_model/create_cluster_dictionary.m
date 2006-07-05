function [Clust, Comps] = create_cluster_dictionary(Clust, Comps, ...
                                   num, P)
% CREATE_CLUSTER_DICTIONARY Create a list and counts of blobs (chars)
%
%   [Clust, Comps] = CREATE_CLUSTER_DICTIONARY(Clust, Comps, [num], ...
%                             [Pos])
%
%   This function uses already clustered page data to come up with baseline
%   offsets and counts of cluster transitions for use in constructing a 
%   dictionary.
%
%   num is optional, and if specified determines the number of clusters (i.e.
%   characters) to use.  This is done based on a descending order of number of 
%   elements belonging to that cluster.  If not given, it defaults to 50
%
%   Pos is optional, and if specified it should be a cell array of nx4 matrices
%   giving line bounding boxes (indexed by page).  It should have the same
%   number of entries as Comps, and should be like Pos returned in get_lines.
%   If not passed, it will be calculated in this function (which takes longer).
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: create_cluster_dictionary.m,v 1.4 2006-07-05 01:00:04 scottl Exp $
%
% REVISION HISTORY
% $Log: create_cluster_dictionary.m,v $
% Revision 1.4  2006-07-05 01:00:04  scottl
% re-written after changing cluster and component structures.  Character counts,
% and bigram model is now stored in Clust.
%
% Revision 1.3  2006/06/21 21:47:12  scottl
% use if test and iterate over all items instead of using i_offs
%
% Revision 1.2  2006/06/19 21:37:57  scottl
% implemented baseline offset calculation, addition of a 'space' cluster.
%
% Revision 1.1  2006/06/12 20:57:50  scottl
% initial revision
%


% LOCAL VARS %
%%%%%%%%%%%%%%
num_clusts = 50;
baseline_thresh = 0.2;
intensity_thresh = 0.5;  %min value for cluster pixel be considered 'on'
bg_val = 0;

add_space = true;  %should we add a space character bitmap and bigram count?
space_width = 12;  %typical value between words (the space width)
space_height = 12; %how tall should the space character be from the baseline 

smoothing_counts = 1;  %add plus-one smoothing to ensure all transitions are
                       %possible.
D = {};
Pos = {};


% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 3 || nargin > 4
    error('incorrect number of arguments passed');
elseif nargin >= 3
    num_clusts = num;
    if nargin == 4
        Pos = P;
    end
end

%sort the clusters to ensure they are ordered by number of elements
[Clust, Comps] = sort_clusters(Clust, Comps);

num_pgs = size(Comps.pg_size,1);
if isempty(Pos) || length(Pos) < num_pgs
    %must run get_lines to determine all line boundaries
    fprintf('%.2fs: determining line boundaries of all pages\n', toc);
    Pos = get_lines(Comps);
    fprintf('%.2fs: finished determining line boundaries\n', toc);
end

if ~ Clust.found_offsets
    fprintf('%.2fs: calculating baseline offsets for each cluster\n', toc);
    if ~ Comps.found_offsets
        fprintf('%.2fs: calculating baseline offsets for each component\n',toc);
        for pp=1:num_pgs
            idx = find(Comps.pg == pp);
            M = ~imread(Comps.files{pp});
            for line=1:size(Pos{pp},1)
                lpos = Pos{pp}(line,:);
                lidx = find(Comps.pos(idx,2) >= lpos(2) & ...
                            Comps.pos(idx,4) <= lpos(4));
                lidx = idx(lidx);
                offs = get_baselines(M(lpos(2):lpos(4),lpos(1):lpos(3)) ~= ...
                       bg_val, baseline_thresh);
                Comps.offset(lidx) = (Comps.pos(lidx,4) - lpos(2)) - offs;
            end
        end
        Comps.found_offsets = true;
        fprintf('%.2fs: finished calculating component offsets\n',toc);
    end
    %take the mode of each Component offset belonging to that cluster
    for cc=1:Clust.num
        offs = single(Comps.offset(Clust.comps{cc}));
        off_vals = unique(offs);
        if length(off_vals) == 1
            Clust.offset(cc) = offs(1);
        elseif length(off_vals) > 1
            [Dummy, mf_bin] = max(hist(offs,off_vals));
            Clust.offset(cc) = off_vals(mf_bin);
        end
    end
    Clust.found_offsets = true;
    fprintf('%.2fs: finished calculating baseline offsets\n', toc);
end

%add the 'space character' Cluster to the end if specified
if add_space
    Clust.num = Clust.num + 1;
    Clust.num_comps(Clust.num) = 0;
    Clust.comps{Clust.num} = [];
    Clust.avg{Clust.num} = bg_val + zeros(space_height, space_width);
    Clust.norm_sq(Clust.num) = 0;
    Clust.offset(Clust.num) = 0;
end
fprintf('%.2fs: finished counting chars and creating bitmaps\n', toc);

%now create the bigram counts
idx = find(Comps.nb(:,3) ~= 0);
Trans = double([Comps.clust(idx), Comps.clust(Comps.nb(idx,3))]);
trans_dist = Comps.pos(Trans(:,2),1) - Comps.pos(Trans(:,1),3);
if add_space
    %add space transitions (based on distance)
    blank_clust = Clust.num;
    blank_trans = floor(trans_dist / space_width);
    Clust.num_comps(Clust.num) = sum(blank_trans);
    blank_idx = find(blank_trans > 0);
    Trans = [Trans; blank_clust+zeros(length(blank_idx),1), Trans(blank_idx,2)];
    Trans(blank_idx,2) = blank_clust;
    blank_trans(blank_idx) = blank_trans(blank_idx) - 1;
    num_b_trans = sum(blank_trans(blank_idx));
    Trans = [Trans; blank_clust+zeros(num_b_trans,2)];
end
Clust.num_trans = size(Trans,1);
Clust.bigram = zeros(Clust.num);
for ii=1:size(Trans,1)
    fr = Trans(ii,1); to = Trans(ii,2);
    Clust.bigram(fr,to) = Clust.bigram(fr,to) + 1;
end

%add smothing counts to the totals;
Clust.bigram = Clust.bigram + smoothing_counts;
Z = sum(Clust.bigram,2);
if any(Z == 0)
    %this can happen if a component never transitions to another component
    %i.e. only appeaars on the far right edge of pages.  Smoothing should
    %normally take care of this, but to prevent dividing by 0, we augment
    %the sum
    zero_rows = find(Z == 0);
    warning('No transitions seen from cluster %d\n', zero_rows);
    Z(zero_rows) = 1;
end
Clust.bigram = Clust.bigram ./ repmat(Z,1,Clust.num);

fprintf('%.2fs: finished creating character bigram matrix\n', toc);
