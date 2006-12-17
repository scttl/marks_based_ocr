function [Clust, Comps] = match_refine(Clust, Comps, varargin)
% match_refine  Attempt to match near-identical Clusters
%
%   [Clust, Comps, chg_list] = match_refine(Clust, Comps, [VAR1, VAL1]...)
%
%   Clust should be a struct containing several fields specifying which 
%   components belong to each cluster, and their averages etc.  See 
%   cluster_comps for details.
%
%   Comps should be a struct containing component information.  See 
%   get_comps for details.
%
%   optional parameters defined in LOCAL VARS below can be overridden by
%   passing in a string representation of the variable name in VAR1, and its
%   new value in VAL1.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: match_refine.m,v 1.9 2006-12-17 20:13:52 scottl Exp $
%
% REVISION HISTORY
% $Log: match_refine.m,v $
% Revision 1.9  2006-12-17 20:13:52  scottl
% show the 5 farthest matches that are still within the threshold
% (intead of the 5 closest), when displaying matches.
%
% Revision 1.8  2006-11-13 17:56:47  scottl
% small spacing improvements.
%
% Revision 1.7  2006/10/18 15:39:22  scottl
% optional argument reorganization.
%
% Revision 1.6  2006/08/30 17:38:23  scottl
% implement batches for Hausdorff matching to prevent memory issues for large
% cluster sizes.
%
% Revision 1.5  2006/08/24 21:40:04  scottl
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

%process at most this many clusters at a time to prevent memory issues
haus_batchsize = 2000; 

%by default take the average when combining clusters (instead of the mode).
use_avg = true;

display_images = false;  %set this to true to display matches as they are found


% CODE START %
%%%%%%%%%%%%%%
if nargin < 2
    error('incorrect number of arguments specified!');
elseif nargin > 2
    process_optional_args(varargin{:});
end

%go through each unrefined cluster and attempt to group it with other clusters
rr = find(Clust.refined == false, 1, 'first');
while ~isempty(rr)
    fprintf('                                                         \r');
    fprintf('cluster: %d  -- ', rr);

    if strcmp(dist_metric, 'euc')
        D = euc_dist(Clust.avg{rr},Clust.avg,Clust.norm_sq(rr),Clust.norm_sq);
    elseif strcmp(dist_metric, 'hausdorff')
        recs_rem = Clust.num;
        D = [];
        while recs_rem > haus_batchsize
            D = [hausdorff_dist(Clust.avg{rr}, ...
                 Clust.avg(recs_rem-haus_batchsize+1:recs_rem)); D];
            recs_rem = recs_rem - haus_batchsize;
        end
        D = [hausdorff_dist(Clust.avg{rr}, Clust.avg(1:recs_rem)); D];
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
            num_match = min(num_match, 5);  %only show the last 5 matches
            subplot(1,num_match+1,1), imshow(Clust.avg{rr});
            xlabel('r'); title('Main Cluster');
            for ii=1:num_match
                pos = length(match_idcs) - ii + 1;
                subplot(1,1+num_match,1+ii), imshow(Clust.avg{match_idcs(pos)});
                xlabel(match_idcs(pos)); 
                title(sprintf('dist: %.3f', D(match_idcs(pos))));
            end
            drawnow;
            pause(.5);
        end
        %merge the clusters together, note that this will re-order the clusters
        [Clust, Comps,idx] = add_and_reaverage(Clust, Comps, rr, ...
                             match_idcs, 'use_avg', use_avg);
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
fprintf('\n');
