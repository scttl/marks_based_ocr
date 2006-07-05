function [Clust, Comps] = match_refine(Clust, Comps, dm, thr)
% match_refine  Attempt to match near-identical Clusters
%
%   [Clust, Comps, chg_list] = match_refine(Clust, Comps, refine_list, ...
%                                           dist_metric, thresh)
%
%   Clust should be a struct containing several fields specifying which 
%   components belong to each cluster, and their averages etc.  See 
%   cluster_comps for details.
%
%   Comps should be a struct containing component information.  See 
%   cluster_comps for details.
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

% CVS INFO %
%%%%%%%%%%%%
% $Id: match_refine.m,v 1.2 2006-07-05 01:16:58 scottl Exp $
%
% REVISION HISTORY
% $Log: match_refine.m,v $
% Revision 1.2  2006-07-05 01:16:58  scottl
% rewritten based on new Cluster and Component structures.
%
% Revision 1.1  2006/06/03 20:55:48  scottl
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
if nargin < 2 || nargin > 4
    error('incorrect number of arguments specified!');
elseif nargin >= 3
    dist_metric = dm;
    if nargin >= 4
        distance_thresh = thr;
    end
end

%go through each unrefined cluster and attempt to group it with other clusters
rr = find(Clust.refined == false, 1, 'first');
while ~isempty(rr)
    fprintf('                                                         \r');
    fprintf('cluster: %d  -- ', rr);

    D = euc_dist(Clust.avg{rr}, Clust.avg, Clust.norm_sq(rr), Clust.norm_sq);
%        elseif strcmp(dist_metric, 'conv_euc')
%            [match_found, Mr, Mi, d] = conv_euc_match(Clust(r).avg, ...
%                                   Clust(i).avg, distance_thresh);
%        elseif strcmp(dist_metric, 'hausdorff')
%            [match_found, Mr, Mi, d] = hausdorff_match(Clust(r).avg, ...
%                                   Clust(i).avg, distance_thresh);
%        elseif strcmp(dist_metric, 'ham')
%            [match_found, Mr, Mi, d] = ham_match(Clust(r).avg, ...
%                                   Clust(i).avg, distance_thresh);
%        else
%            error('incorrect distance metric specified!');
%        end
    match_idcs = find(D <= distance_thresh);
    match_idcs = match_idcs(match_idcs ~= rr);
    num_match = length(match_idcs);

    if num_match > 0
        fprintf('found %d matches for cluster %d\r', num_match, rr);
        if display_images
            clf;
            num_match = min(num_match, 5);  %only show the first 5 matches
            subplot(1,num_match+1,1), imshow(Clust.avg{rr});
            xlabel('r'); title('Main Cluster');
            for ii=1:num_match
                subplot(1,1+num_match,1+ii), imshow(Clust.avg{match_idcs(ii)});
                xlabel(match_idcs(ii)); 
                title(sprintf('dist: %.3f', D(match_idcs(ii))));
            end
            drawnow;
            pause(.5);
        end
        %merge the clusters together, note that this will re-order the clusters
        [Clust, Comps] = add_and_reaverage(Clust, Comps, rr, match_idcs);
    else
        Clust.refined(rr) = true;
        fprintf('no match\r');
    end
    rr = find(Clust.refined == false, 1, 'first');
end
