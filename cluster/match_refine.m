function [Clust, chg_list] = match_refine(Clust, rl, dm, thr)
% match_refine  Attempt to match near-identical Clusters
%
%   [Clust, chg_list] = match_refine(Clust, refine_list, dist_metric, thresh)
%
%   Clust should be an array of structs, each of which is assumed to contain 
%   several fields of a particular format.  See cluster_comps for details.
%
%   refine_list is optional and if specified, determines which items are
%   to be refined (ie searched for a match).  It should be a vector of cluster
%   indices and will be processed in the order given.  If not specified, this 
%   method will attempt to refine all clusters, starting from the highest 
%   numbered one.
%
%   dist_metric is optional and if specified determines the type of distance
%   metric used for matching.  Valid options for this parameter are 'euc'
%   (straight Euclidian distance -- the default), 'conv_euc' (Euclidian distance
%   after convolving the matrices to find maximal overlapping point), 
%   'hausdorff' to use Hausdorff distance, or 'ham' for Hamming distance
%
%   thresh is optional and if specified, determines the normalized maximal 
%   Euclidian distance allowable for two cluster averages to be considered 
%   matching.  If not specified, it defaults to .009
%
%   the refined list of clusters is returned in Clust, as is the indices of 
%   those clusters that changed in chg_list

% CVS INFO %
%%%%%%%%%%%%
% $Id: match_refine.m,v 1.1 2006-06-03 20:55:48 scottl Exp $
%
% REVISION HISTORY
% $Log: match_refine.m,v $
% Revision 1.1  2006-06-03 20:55:48  scottl
% Initial check-in.
%
%

% LOCAL VARS %
%%%%%%%%%%%%%%
dist_metric = 'euc';              %default distance metric
distance_thresh = 0.009;          %default distance theshold

display_images = false;  %set this to true to display matches as they are found


% CODE START %
%%%%%%%%%%%%%%
if nargin < 1 || nargin > 4
    error('incorrect number of arguments specified!');
elseif nargin >= 2
    refine_list = rl;
    if nargin >= 3
        dist_metric = dm;
        if nargin == 4
            distance_thresh = thr;
        end
    end
end

num_clusts = size(Clust, 1);
chg_list = [];

%go through each cluster to see if it matches another cluster
if nargin < 2
    refine_list = num_clusts:-1:1;
end
keep_list = 1:num_clusts;

while ~isempty(refine_list)
    r = refine_list(1);
    fprintf('                                                         \r');
    fprintf('cluster: %d  -- ', r);
    match_found = false;

    for i=keep_list
        if i == r
            continue;
        end

        if strcmp(dist_metric,'euc')
            [match_found, Mr, Mi, d] = euc_match(Clust(r).avg, Clust(i).avg, ...
                                   distance_thresh);
        elseif strcmp(dist_metric, 'conv_euc')
            [match_found, Mr, Mi, d] = conv_euc_match(Clust(r).avg, ...
                                   Clust(i).avg, distance_thresh);
        elseif strcmp(dist_metric, 'hausdorff')
            [match_found, Mr, Mi, d] = hausdorff_match(Clust(r).avg, ...
                                   Clust(i).avg, distance_thresh);
        elseif strcmp(dist_metric, 'ham')
            [match_found, Mr, Mi, d] = ham_match(Clust(r).avg, ...
                                   Clust(i).avg, distance_thresh);
        else
            error('incorrect distance metric specified!');
        end

        if match_found
            fprintf('found a match between %d and %d\r', r, i);
            if display_images
                clf;
                subplot(1,2,1), imshow(Clust(r).avg), xlabel('r'), ...
                title(sprintf('distance: %f', d));
                subplot(1,2,2), imshow(Clust(i).avg), xlabel('i');
                drawnow;
                pause(.5);
            end

            %merge the clusters together
            Clust(i) = add_and_reaverage(Clust(i), Clust(r), Mi, Mr);
            pos = find(keep_list == r,1);
            keep_list = [keep_list(1:pos-1), keep_list(pos+1:end)];
            chg_list = [chg_list, i];
            break;
        end
    end
    if ~ match_found
        fprintf('no match\r');
    end
    refine_list = refine_list(2:end);
end

%refinement complete, now delete those items not on the keep list
[Dummy, Dummy, chg_list] = intersect(unique(chg_list), keep_list);
Clust = Clust(keep_list);
