function [Clust, Comps] = match_refine(Clust, Comps, dm, thr)
% match_refine  Attempt to match near-identical Clusters
%
%   [Clust, Comps, chg_list] = match_refine(Clust, Comps, [dist_metric, thresh])
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
%   thresh is optional and if specified, determines the maximal 
%   distance allowable for two cluster averages to be considered 
%   matching.  If not specified, it defaults to .009 (suitable for use with
%   'euc' as a distance metric).
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: match_refine.m,v 1.5 2006-08-24 21:40:04 scottl Exp $
%
% REVISION HISTORY
% $Log: match_refine.m,v $
% Revision 1.5  2006-08-24 21:40:04  scottl
% added ability to use the mode instead of taking the average of cluster
% intensities while refining.
%
% Revision 1.4  2006/08/14 01:34:34  scottl
% Updates based on new Clust.changed field, and changes to add_and_reaverage.
%
% Revision 1.3  2006/07/22 04:10:30  scottl
% cleaned up comments, implemented alternate distance metrics, modified
% add_and_reaverage
%
% Revision 1.2  2006/07/05 01:16:58  scottl
% rewritten based on new Cluster and Component structures.
%
% Revision 1.1  2006/06/03 20:55:48  scottl
% Initial check-in.


% LOCAL VARS %
%%%%%%%%%%%%%%
dist_metric = 'euc';              %default distance metric
distance_thresh = 0.009;          %default distance threshold

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

    if strcmp(dist_metric, 'euc')
        D = euc_dist(Clust.avg{rr},Clust.avg,Clust.norm_sq(rr),Clust.norm_sq);
    elseif strcmp(dist_metric, 'hausdorff')
        D = hausdorff_dist(Clust.avg{rr}, Clust.avg);
    elseif strcmp(dist_metric, 'ham')
        D = ham_dist(Clust.avg{rr}, Clust.avg);
    elseif strcmp(dist_metric, 'conv_euc')
        D = conv_euc_dist(Clust.avg{rr},Clust.avg,Clust.norm_sq(rr), ...
            Clust.norm_sq);
    else
        error('incorrect distance metric specified!');
    end
    match_idcs = find(D <= distance_thresh);
    match_idcs = match_idcs(match_idcs ~= rr);
    %sort matches in increasing distance order
    [srt_idx,srt_idx] = sort(D(match_idcs));
    match_idcs = match_idcs(srt_idx);
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
        if strcmp(dist_metric, 'hausdorff')
            %take the mode instead of averaging pixel intensenties
            [Clust, Comps,idx] = add_and_reaverage(Clust, Comps, rr, ...
                                 match_idcs, true);
        else
            [Clust, Comps,idx] = add_and_reaverage(Clust,Comps,rr,match_idcs);
        end
        %mark the new cluster as changed, so it will be scanned on subsequent
        %checks
        Clust.changed(idx) = true;
        Clust.refined(idx) = true;
    else
        Clust.refined(rr) = true;
        fprintf('no match\r');
    end
    rr = find(Clust.refined == false, 1, 'first');
end
