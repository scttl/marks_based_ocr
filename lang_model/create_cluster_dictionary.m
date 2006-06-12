function [D, bitmaps] = create_cluster_dictionary_dictionary(Clust, Comps, num)
% CREATE_CLUSTER_DICTIONARY Create a list and counts of blobs (chars)
%
%   [D, bitmaps] = CREATE_CLUSTER_DICTIONARY(Clust, Comps, [num])
%   Files should be a cell array (or single string), and give the full 
%   path and name of the file(s) to be used in creating the dictionary.
%
%   num is optional, and if specified determines the number of clusters (i.e.
%   characters) to use.  This is done based on a descending order of number of 
%   elements belonging to that cluster.  If not given, it defaults to 50
%
%   D is a struct with fields: char, char_count, char_bigram
%
%   bitmaps is a cell array of logical arrays showing the individual clusters
%   (i.e. characters).
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: create_cluster_dictionary.m,v 1.1 2006-06-12 20:57:50 scottl Exp $
%
% REVISION HISTORY
% $Log: create_cluster_dictionary.m,v $
% Revision 1.1  2006-06-12 20:57:50  scottl
% initial revision
%


% LOCAL VARS %
%%%%%%%%%%%%%%
num_clusts = 50;
avg_thresh = 0.5;  %min value for cluster pixel intensity be considered 'on'
D = {};


% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 2 || nargin > 3
    error('incorrect number of arguments passed');
elseif nargin == 3
    num_clusts = num;
end

%sort the clusters to ensure they are ordered by number of elements
Clust = sort_clusters(Clust);

for i=1:num_clusts
    D.char_id(i) = i;
    D.char_count(i) = Clust(i).num;
    bitmaps{i} = zeros(size(Clust(i).avg));
    bitmaps{i}(find(Clust(i).avg >= avg_thresh)) = 1;
end
fprintf('%.2fs: finished counting chars and creating bitmaps\n', toc);

D.char_bigram = zeros(num_clusts);
for i = 1:num_clusts
    right_nbs = Clust(i).nb(:,3);
    right_nbs = right_nbs(right_nbs > 0);
    for j=1:length(right_nbs)
        [cl, off] = get_cl_off(Clust, right_nbs(j));
        if cl <= num_clusts
            D.char_bigram(i, cl) = D.char_bigram(i,cl) + 1;
        end
    end
end
fprintf('%.2fs: finished creating character bigram matrix\n', toc);
