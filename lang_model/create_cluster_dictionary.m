function [D, bitmaps, base_offs] = create_cluster_dictionary(Clust, Comps, ...
                                   num, P)
% CREATE_CLUSTER_DICTIONARY Create a list and counts of blobs (chars)
%
%   [D, bitmaps, base_offs] = CREATE_CLUSTER_DICTIONARY(Clust, Comps, [num], ...
%                             [Pos])
%   Files should be a cell array (or single string), and give the full 
%   path and name of the file(s) to be used in creating the dictionary.
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
%   D is a struct with fields: char, char_count, char_bigram
%
%   bitmaps is a cell array of logical arrays showing the individual clusters
%   (i.e. characters).
%
%   base_offs is a vector giving positive or negative vertical offsets of the
%   individual clusters (i.e. characters) from the baseline of each line of
%   text.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: create_cluster_dictionary.m,v 1.3 2006-06-21 21:47:12 scottl Exp $
%
% REVISION HISTORY
% $Log: create_cluster_dictionary.m,v $
% Revision 1.3  2006-06-21 21:47:12  scottl
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
Clust = sort_clusters(Clust);

num_pgs = length(Comps);
if isempty(Pos) || length(Pos) < num_pgs
    %must run get_lines to determine all line boundaries
    fprintf('%.2fs: determining line boundaries of all pages\n', toc);
    Pos = get_lines(Clust, Comps);
    fprintf('%.2fs: finished determining line boundaries\n', toc);
end

fprintf('%.2fs: calculating baseline offsets for each cluster\n', toc);
base_offs = cell(1,num_clusts);
for p=1:num_pgs
    for line=1:size(Pos{p},1)
        l = Pos{p}(line,1); t = Pos{p}(line,2); r = Pos{p}(line,3);
        b = Pos{p}(line,4);
        baseline = get_baselines(Comps{p}(t:b,l:r) ~= bg_val, baseline_thresh);

        comp_ids = unique(Comps{p}(t:b, l:r));
        comp_ids = comp_ids(comp_ids ~= bg_val);
        for i=1:length(comp_ids)
            [cl, off] = get_cl_off(Clust, comp_ids(i));
            if cl <= num_clusts
                base_offs{cl} = [base_offs{cl}, Clust(cl).pos(off,4) - ...
                                (t-1+baseline)];
            end
        end
    end
end

%now take the mode of each array of offsets to determine the best overall
%offset
for i=1:length(base_offs)
    max_bin = max(base_offs{i});
    min_bin = min(base_offs{i});
    [Dummy, idx] = max(hist(base_offs{i}, max_bin - min_bin + 1));
    base_offs{i} = min_bin + idx - 1;
end
base_offs = cell2mat(base_offs);
fprintf('%.2fs: finished calculating baseline offsets\n', toc);

for i=1:num_clusts
    D.char_id(i) = i;
    D.char_count(i) = Clust(i).num;
    bitmaps{i} = zeros(size(Clust(i).avg));
    bitmaps{i}(find(Clust(i).avg >= intensity_thresh)) = 1;
end

%add the 'space character' bitmap to the end if specified
if add_space
    bitmaps{num_clusts+1} = bg_val + zeros(space_height, space_width);
    base_offs(num_clusts+1) = 0;
end
fprintf('%.2fs: finished counting chars and creating bitmaps\n', toc);

%now create the bigram counts
if add_space
    D.char_bigram = zeros(num_clusts+1);
    D.char_count(num_clusts+1) = 0;
else
    D.char_bigram = zeros(num_clusts);
end
for i = 1:num_clusts
    right_nbs = Clust(i).nb(:,3);
    for j=1:length(right_nbs)
        if right_nbs(j) == 0
            continue;
        end
        [cl, off] = get_cl_off(Clust, right_nbs(j));
        trans_clust = i;
        if add_space
            %check to see if 1 or more spaces exist between this character and
            %its right neighbour
            dist = Clust(cl).pos(off,1) - Clust(i).pos(j,3);
            while dist >= space_width
                D.char_bigram(trans_clust,num_clusts+1) = ...
                              D.char_bigram(trans_clust,num_clusts+1) + 1;
                D.char_count(num_clusts+1) = D.char_count(num_clusts+1) +1;
                dist = dist - space_width;
                trans_clust = num_clusts + 1;
            end
        end
        if cl <= num_clusts
            D.char_bigram(trans_clust, cl) = D.char_bigram(trans_clust,cl) + 1;
        end
    end
end
fprintf('%.2fs: finished creating character bigram matrix\n', toc);
