function [Clust, Comps] = merge_refine(Clust,Comps,varargin)
% MERGE_REFINE  Merge together inappropriately separated clusters
%
%   [CLUST, COMPS] = merge_refine(CLUST, COMPS, [VAR1, VAL1]...) 
%   CLUST should be a struct containing cluster information.  See cluster_comps
%   for details on the format of each field.
%
%   COMPS should be a struct containing information for each component.  See 
%   get_comps for details.
%
%   Default values for variables defined in LOCAL VARS can be overriden by 
%   specifying the variable name as a string in VAR1 and its new value in 
%   VAL1 (and repeating these pairs for any other variables whose values are 
%   to be overridden).
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: merge_refine.m,v 1.7 2006-10-29 17:24:54 scottl Exp $
%
% REVISION HISTORY
% $Log: merge_refine.m,v $
% Revision 1.7  2006-10-29 17:24:54  scottl
% change to cluster struct, to use descender and ascender offsets, instead
% of a single offset field.
%
% Revision 1.6  2006/10/18 15:47:17  scottl
% changes to introduce scale invariant matches, better parameter
% processing
%
% Revision 1.5  2006/08/24 21:40:04  scottl
% added ability to use the mode instead of taking the average of cluster
% intensities while refining.
%
% Revision 1.4  2006/08/14 01:36:43  scottl
% finished re-implementation based on new Clust and Comps structs.
%
% Revision 1.3  2006/07/05 01:21:23  scottl
% started rewriting based on new Cluster and Component structures.  Not
% complete yet!
%
% Revision 1.2  2006/06/12 20:56:01  scottl
% changed comps to a cell array (indexed by page) to allow different sized pages
%
% Revision 1.1  2006/06/03 20:55:48  scottl
% Initial check-in.


% LOCAL VARS %
%%%%%%%%%%%%%%
max_dist = 3;
valid_pct = .85;
min_match_comp = 3;

%should we display matches onscreen (wastes resources)
display_matches = false;

%should we thin merged clusters when created?
use_thinned_imgs = true;

%should we resize merged clusters when created?
resize_imgs = false;
resize_method = 'nearest';

%these will contain the list of extraneous components and clusters to be
%deleted at the end of processing
del_comps = [];
del_clusts = [];


% CODE START %
%%%%%%%%%%%%%%

if nargin < 2
    error('incorrect number of arguments specified!');
elseif nargin > 2
    process_optional_args(varargin{:});
end

%initially put all elements in the refine list if not passed
rr = find(Clust.refined == false,1,'first');
while ~ isempty(rr)

    %determine if this item can be refined by merging it with a cluster
    %that almost always appears adjacent to it
    fprintf('                                                        \r');
    fprintf('cluster: %d -- ', rr);
    if Clust.num_comps(rr) < min_match_comp
        %too few components in this cluster for a valid merge
        Clust.refined(rr) = true;
        rr = find(Clust.refined == false,1,'first');
        continue;
    end
    rr_comps = Clust.comps{rr};

    %first look for a valid left neighbour cluster
    r_comps = rr_comps;
    lnb_comps = Comps.nb(r_comps,1);
    valid_idx = lnb_comps ~= 0;
    r_comps = r_comps(valid_idx);
    lnb_comps = lnb_comps(valid_idx);
    l_dist = Comps.pos(r_comps,1) - Comps.pos(lnb_comps,3);
    valid_idx = l_dist <= max_dist;
    r_comps = r_comps(valid_idx);
    lnb_comps = lnb_comps(valid_idx);
    if length(lnb_comps) >= min_match_comp
        %determine the most frequent cluster of the remaining valid
        %left neighbour components
        lnb_clusts = single(Comps.clust(lnb_comps));
        [mf_clust, mf_num] = mode(lnb_clusts);

        %ensure that the most frequent lnb cluster has the required minimum
        %number of components, makes up at least valid_pct of the total number 
        %of that clusters components as well as valid_pct of r's total number of
        %components.
        if mf_num >= min_match_comp && ...
           (mf_num / Clust.num_comps(mf_clust)) >= valid_pct && ...
           (mf_num / Clust.num_comps(rr)) >= valid_pct
            %valid merge found, update neighbours, cluster averages etc.
            fprintf('valid left merge found\r');
            if display_matches
                clf;
                subplot(1,2,1), imshow(Clust.avg{mf_clust});
                xlabel(Clust.num_comps(mf_clust));
                subplot(1,2,2), imshow(Clust.avg{rr});
                xlabel(Clust.num_comps(rr));
                drawnow;
                pause(.5);
            end

            merge_idx = Comps.clust(lnb_comps) == mf_clust;
            r_merge_comps = r_comps(merge_idx);
            lnb_merge_comps = lnb_comps(merge_idx);

            %update the merged components.  Note that the right pos, neighbour
            %and neighbour distance don't change
            %positions
            Comps.pos(r_merge_comps,1) = Comps.pos(lnb_merge_comps,1);
            Comps.pos(r_merge_comps,2) = min(...
               Comps.pos(r_merge_comps,2), Comps.pos(lnb_merge_comps,2));
            Comps.pos(r_merge_comps,4) = max(...
               Comps.pos(r_merge_comps,4), Comps.pos(lnb_merge_comps,4));

            %neighbours
            Comps.nb(r_merge_comps,1) = Comps.nb(lnb_merge_comps,1);
            Comps.nb_dist(r_merge_comps,1) = Comps.nb_dist(lnb_merge_comps,1);
            %for top and bottom, take the closer of the two possible neighbours
            %provided it isn't one of the components being merged.
            [dist,idx] = min([Comps.nb_dist(r_merge_comps,2), ...
                             Comps.nb_dist(lnb_merge_comps,2)],[],2);
            top_nbs = [Comps.nb(r_merge_comps,2), Comps.nb(lnb_merge_comps,2)];
            top_r_rejects = Comps.nb(r_merge_comps,2) == lnb_merge_comps;
            top_lnb_rejects = Comps.nb(lnb_merge_comps,2) == r_merge_comps;
            flip_r = top_r_rejects & (idx == 1);
            dist(flip_r) = Comps.nb_dist(lnb_merge_comps(flip_r),2);
            flip_lnb = top_lnb_rejects & (idx == 2);
            dist(flip_lnb) = Comps.nb_dist(r_merge_comps(flip_lnb),2);
            idx(flip_r) = 2;
            idx(flip_lnb) = 1;
            Comps.nb_dist(r_merge_comps,2) = dist;
            Comps.nb(r_merge_comps,2) = top_nbs(sub2ind(size(top_nbs), ...
                                        (1:length(idx))', idx));
            [dist,idx] = min([Comps.nb_dist(r_merge_comps,4), ...
                             Comps.nb_dist(lnb_merge_comps,4)],[],2);
            bot_nbs = [Comps.nb(r_merge_comps,4), Comps.nb(lnb_merge_comps,4)];
            bot_r_rejects = Comps.nb(r_merge_comps,4) == lnb_merge_comps;
            bot_lnb_rejects = Comps.nb(lnb_merge_comps,4) == r_merge_comps;
            flip_r = bot_r_rejects & (idx == 1);
            dist(flip_r) = Comps.nb_dist(lnb_merge_comps(flip_r),4);
            flip_lnb = bot_lnb_rejects & (idx == 2);
            dist(flip_lnb) = Comps.nb_dist(r_merge_comps(flip_lnb),4);
            idx(flip_r) = 2;
            idx(flip_lnb) = 1;
            Comps.nb_dist(r_merge_comps,4) = dist;
            Comps.nb(r_merge_comps,4) = bot_nbs(sub2ind(size(bot_nbs), ...
                                        (1:length(idx))', idx));
            %must also update any component that lists the left half of the now
            %merged components as a neighbour
            idx = (1:Comps.max_comp)';
            idx(lnb_merge_comps) = r_merge_comps;
            valid = Comps.nb ~= 0;
            Comps.nb(valid) = idx(Comps.nb(valid));

            %offsets
            if Comps.found_lines
                Comps.descender_off(r_merge_comps) = max(...
                                        Comps.descender_off(r_merge_comps), ...
                                        Comps.descender_off(lnb_merge_comps));
                Comps.ascender_off(r_merge_comps) = min(...
                                        Comps.ascender_off(r_merge_comps), ...
                                        Comps.ascender_off(lnb_merge_comps));
                %scale factor
                Comps.scale_factor(r_merge_comps) = Comps.modal_height ./ ...
                                   double(Comps.pos(r_merge_comps,4) - ...
                                   Comps.pos(r_merge_comps,2) + 1);
            end

            %ground truth assignment (doesn't change)

            %clusters (we'll require a new cluster for the merged components
            %if at least one of the components from the original right half 
            %cluster is not being merged.
            merge_clust = rr;
            if length(r_merge_comps) ~= Clust.num_comps(rr)
                Clust.num = Clust.num + 1;
                merge_clust = Clust.num;
                Comps.clust(r_merge_comps) = merge_clust;
                r_unmerged_comps = setdiff(Clust.comps{rr}, r_merge_comps);
                Clust.num_comps(rr) = length(r_unmerged_comps);
                Clust.comps{rr} = r_unmerged_comps;
                Clust.refined(rr) = false;
                Clust.changed(rr) = true;
                Clust.num_comps(merge_clust) = length(r_merge_comps);
                Clust.comps{merge_clust} = r_merge_comps;
                if ~isempty(Clust.bigram)
                    Clust.bigram(Clust.num,:) = NaN;
                    Clust.bigram(:,Clust.num) = NaN;
                end
                %recalculate average for rr
                Clust.mode_num(rr) = length(r_unmerged_comps);
                Clust.avg{rr} = recalc_avg(Comps, r_unmerged_comps, ...
                                use_thinned_imgs, resize_imgs, resize_method);
                Clust.norm_sq(rr) = sum(Clust.avg{rr}(:).^2);
                if Comps.found_lines
                    Clust.descender_off(rr) = int16(mode(single(...
                                   Comps.descender_off(r_unmerged_comps))));
                    Clust.ascender_off(rr) = int16(mode(single(...
                                   Comps.ascender_off(r_unmerged_comps))));
                end
            end
            if length(lnb_merge_comps) ~= Clust.num_comps(mf_clust)
                %some unmerged left-halfs remain
                lnb_unmerged_comps = setdiff(Clust.comps{mf_clust}, ...
                                     lnb_merge_comps);
                Clust.num_comps(mf_clust) = length(lnb_unmerged_comps);
                Clust.comps{mf_clust} = lnb_unmerged_comps;
                Clust.refined(mf_clust) = false;
                Clust.changed(mf_clust) = true;
                %recalculate average for mf_clust
                Clust.mode_num(mf_clust) = length(lnb_unmerged_comps);
                Clust.avg{mf_clust} = recalc_avg(Comps, lnb_unmerged_comps, ...
                                 use_thinned_imgs, resize_imgs, resize_method);
                Clust.norm_sq(mf_clust) = sum(Clust.avg{mf_clust}(:).^2);
                if Comps.found_lines
                    Clust.descender_off(mf_clust) = int16(mode(single(...
                                     Comps.descender_off(lnb_unmerged_comps))));
                    Clust.ascender_off(mf_clust) = int16(mode(single(...
                                     Comps.ascender_off(lnb_unmerged_comps))));
                end
            else
                %flag the left-half cluster for removal
                Clust.refined(mf_clust) = true;
                del_clusts = [del_clusts; mf_clust];
            end
            %now update the merged cluster average and offsets based on the
            %components.  Unfortunately this is computationally expensive
            %if there are a lot of components
            Clust.mode_num(merge_clust) = length(r_merge_comps);
            Clust.avg{merge_clust} = recalc_avg(Comps, r_merge_comps, ...
                                 use_thinned_imgs, resize_imgs, resize_method);
            Clust.norm_sq(merge_clust) = sum(Clust.avg{merge_clust}(:).^2);
            if Comps.found_lines
                Clust.descender_off(merge_clust) = int16(mode(single(...
                                        Comps.descender_off(r_merge_comps))));
                Clust.ascender_off(merge_clust) = int16(mode(single(...
                                        Comps.ascender_off(r_merge_comps))));
            end
            Clust.refined(merge_clust) = true;
            Clust.changed(merge_clust) = true;

            %flag the left-half components for removal
            del_comps = [del_comps; lnb_merge_comps];

            rr = find(Clust.refined == false,1,'first');
            continue;  %don't try and find a top match too
        else
            %no valid left neighbour match found for this cluster.
            Clust.refined(rr) = true;
        end
    end

    %now look for a valid top neighbour cluster
    r_comps = rr_comps;
    tnb_comps = Comps.nb(r_comps,2);
    valid_idx = tnb_comps ~= 0;
    r_comps = r_comps(valid_idx);
    tnb_comps = tnb_comps(valid_idx);
    t_dist = Comps.pos(r_comps,2) - Comps.pos(tnb_comps,4);
    valid_idx = t_dist <= max_dist;
    r_comps = r_comps(valid_idx);
    tnb_comps = tnb_comps(valid_idx);
    if length(tnb_comps) >= min_match_comp
        %determine the most frequent cluster of the remaining valid
        %top neighbour components
        tnb_clusts = single(Comps.clust(tnb_comps));
        [mf_clust, mf_num] = mode(tnb_clusts);

        %ensure that the most frequent tnb cluster has the required minimum
        %number of components, makes up at least valid_pct of the total number 
        %of that clusters components as well as valid_pct of r's total number of
        %components.
        if mf_num >= min_match_comp && ...
           (mf_num / Clust.num_comps(mf_clust)) >= valid_pct && ...
           (mf_num / Clust.num_comps(rr)) >= valid_pct
            %valid merge found, update neighbours, cluster averages etc.
            fprintf('valid top merge found\r');
            if display_matches
                clf;
                subplot(2,1,1), imshow(Clust.avg{mf_clust});
                xlabel(Clust.num_comps(mf_clust));
                subplot(2,1,2), imshow(Clust.avg{rr});
                xlabel(Clust.num_comps(rr));
                drawnow;
                pause(.5);
            end

            merge_idx = Comps.clust(tnb_comps) == mf_clust;
            r_merge_comps = r_comps(merge_idx);
            tnb_merge_comps = tnb_comps(merge_idx);

            %update the merged components.  Note that the bottom pos, neighbour
            %and neighbour distance don't change
            %positions
            Comps.pos(r_merge_comps,2) = Comps.pos(tnb_merge_comps,2);
            Comps.pos(r_merge_comps,1) = min(...
               Comps.pos(r_merge_comps,1), Comps.pos(tnb_merge_comps,1));
            Comps.pos(r_merge_comps,3) = max(...
               Comps.pos(r_merge_comps,3), Comps.pos(tnb_merge_comps,3));

            %neighbours
            Comps.nb(r_merge_comps,2) = Comps.nb(tnb_merge_comps,2);
            Comps.nb_dist(r_merge_comps,2) = Comps.nb_dist(tnb_merge_comps,2);
            %for left and right, take the closer of the two possible neighbours.
            [dist,idx] = min([Comps.nb_dist(r_merge_comps,1), ...
                             Comps.nb_dist(tnb_merge_comps,1)],[],2);
            lft_nbs = [Comps.nb(r_merge_comps,1), Comps.nb(tnb_merge_comps,1)];
            lft_r_rejects = Comps.nb(r_merge_comps,1) == tnb_merge_comps;
            lft_tnb_rejects = Comps.nb(tnb_merge_comps,1) == r_merge_comps;
            flip_r = lft_r_rejects & (idx == 1);
            dist(flip_r) = Comps.nb_dist(tnb_merge_comps(flip_r),1);
            flip_tnb = lft_tnb_rejects & (idx == 2);
            dist(flip_tnb) = Comps.nb_dist(r_merge_comps(flip_tnb),1);
            idx(flip_r) = 2;
            idx(flip_tnb) = 1;
            Comps.nb_dist(r_merge_comps,1) = dist;
            Comps.nb(r_merge_comps,1) = lft_nbs(sub2ind(size(lft_nbs), ...
                                        (1:length(idx))', idx));
            [dist,idx] = min([Comps.nb_dist(r_merge_comps,3), ...
                             Comps.nb_dist(tnb_merge_comps,3)],[],2);
            rgt_nbs = [Comps.nb(r_merge_comps,3), Comps.nb(tnb_merge_comps,3)];
            rgt_r_rejects = Comps.nb(r_merge_comps,3) == tnb_merge_comps;
            rgt_tnb_rejects = Comps.nb(tnb_merge_comps,3) == r_merge_comps;
            flip_r = rgt_r_rejects & (idx == 1);
            dist(flip_r) = Comps.nb_dist(tnb_merge_comps(flip_r),3);
            flip_tnb = rgt_tnb_rejects & (idx == 2);
            dist(flip_tnb) = Comps.nb_dist(r_merge_comps(flip_tnb),3);
            idx(flip_r) = 2;
            idx(flip_tnb) = 1;
            Comps.nb_dist(r_merge_comps,3) = dist;
            Comps.nb(r_merge_comps,3) = rgt_nbs(sub2ind(size(rgt_nbs), ...
                                        (1:length(idx))', idx));
            %must also update any component that lists the top half of the now
            %merged components as a neighbour
            idx = (1:Comps.max_comp)';
            idx(tnb_merge_comps) = r_merge_comps;
            valid = Comps.nb ~= 0;
            Comps.nb(valid) = idx(Comps.nb(valid));

            %offsets don't change since we'll just use the bottom components,
            %but scale factors do change
            if Comps.found_lines
                %scale factor
                Comps.scale_factor(r_merge_comps) = Comps.modal_height ./ ...
                                   double(Comps.pos(r_merge_comps,4) - ...
                                   Comps.pos(r_merge_comps,2) + 1);
            end

            %clusters (we'll require a new cluster for the merged components
            %if at least one of the components from the original bottom half 
            %cluster is not being merged.
            merge_clust = rr;
            if length(r_merge_comps) ~= Clust.num_comps(rr)
                Clust.num = Clust.num + 1;
                merge_clust = Clust.num;
                Comps.clust(r_merge_comps) = merge_clust;
                r_unmerged_comps = setdiff(Clust.comps{rr}, r_merge_comps);
                Clust.num_comps(rr) = length(r_unmerged_comps);
                Clust.comps{rr} = r_unmerged_comps;
                Clust.refined(rr) = false;
                Clust.changed(rr) = true;
                Clust.num_comps(merge_clust) = length(r_merge_comps);
                Clust.comps{merge_clust} = r_merge_comps;
                if ~isempty(Clust.bigram)
                    Clust.bigram(Clust.num,:) = NaN;
                    Clust.bigram(:,Clust.num) = NaN;
                end
                %recalculate average for rr
                Clust.mode_num(rr) = length(r_unmerged_comps);
                Clust.avg{rr} = recalc_avg(Comps, r_unmerged_comps, ...
                                use_thinned_imgs, resize_imgs, resize_method);
                Clust.norm_sq(rr) = sum(Clust.avg{rr}(:).^2);
                if Comps.found_lines
                    Clust.descender_off(rr) = int16(mode(single(...
                                   Comps.descender_off(r_unmerged_comps))));
                    Clust.ascender_off(rr) = int16(mode(single(...
                                   Comps.ascender_off(r_unmerged_comps))));
                end
            end
            if length(tnb_merge_comps) ~= Clust.num_comps(mf_clust)
                %some unmerged top-halfs remain
                tnb_unmerged_comps = setdiff(Clust.comps{mf_clust}, ...
                                     tnb_merge_comps);
                Clust.num_comps(mf_clust) = length(tnb_unmerged_comps);
                Clust.comps{mf_clust} = tnb_unmerged_comps;
                Clust.refined(mf_clust) = false;
                Clust.changed(mf_clust) = true;
                %recalculate average for mf_clust
                Clust.mode_num(mf_clust) = length(tnb_unmerged_comps);
                Clust.avg{mf_clust} = recalc_avg(Comps, tnb_unmerged_comps, ...
                                use_thinned_imgs, resize_imgs, resize_method);
                Clust.norm_sq(merge_clust) = sum(Clust.avg{mf_clust}(:).^2);
                if Comps.found_lines
                    Clust.descender_off(merge_clust) = int16(mode(single(...
                                    Comps.descender_off(tnb_unmerged_comps))));
                    Clust.ascender_off(merge_clust) = int16(mode(single(...
                                    Comps.ascender_off(tnb_unmerged_comps))));
                end
            else
                %flag the top-half cluster for removal
                Clust.refined(mf_clust) = true;
                del_clusts = [del_clusts; mf_clust];
            end
            %now update the Cluster average and offsets based on these
            %components.  Unfortunately this is computationally expensive
            %if there are a lot of components
            Clust.mode_num(merge_clust) = length(r_merge_comps);
            Clust.avg{merge_clust} = recalc_avg(Comps, r_merge_comps, ...
                                use_thinned_imgs, resize_imgs, resize_method);
            Clust.norm_sq(merge_clust) = sum(Clust.avg{merge_clust}(:).^2);
            if Comps.found_lines
                Clust.descender_off(merge_clust) = int16(mode(single(...
                                        Comps.descender_off(r_merge_comps))));
                Clust.ascender_off(merge_clust) = int16(mode(single(...
                                        Comps.ascender_off(r_merge_comps))));
            end
            Clust.refined(merge_clust) = true;
            Clust.changed(merge_clust) = true;

            %flag the now merged top-neighbour components for removal
            del_comps = [del_comps; tnb_merge_comps];

        else
            %no valid top neighbour match found for this cluster.
            Clust.refined(rr) = true;
        end
    else
        %no valid top or left neighbour match found for this cluster.
        Clust.refined(rr) = true;
    end
    rr = find(Clust.refined == false,1,'first');
end

%now cleanup deleted components and clusters
keep_comps = setdiff(1:Comps.max_comp, del_comps);
%first renumber the neighbours and cluster component refereces based on the 
%deleted components
idx = zeros(Comps.max_comp,1);
idx(keep_comps) = (1:length(keep_comps))';
valid = Comps.nb ~= 0;
Comps.nb(valid) = idx(Comps.nb(valid));
for ii=1:Clust.num
    Clust.comps{ii} = idx(Clust.comps{ii});
end
%now remove the unneccessary components.
Comps.max_comp = length(keep_comps);
Comps.clust = Comps.clust(keep_comps);
Comps.pos = Comps.pos(keep_comps,:);
Comps.pg = Comps.pg(keep_comps);
Comps.nb = Comps.nb(keep_comps,:);
Comps.nb_dist = Comps.nb_dist(keep_comps,:);
if Comps.found_lines
    Comps.descender_off = Comps.descender_off(keep_comps);
    Comps.ascender_off = Comps.ascender_off(keep_comps);
    Comps.scale_factor = Comps.scale_factor(keep_comps);
    Comps.line = Comps.line(keep_comps);
end
if Comps.found_true_labels
    Comps.truth_label = Comps.truth_label(keep_comps);
end

%now remove the unnessary clusters, updating component cluster references based
%on the deleted clusters
keep_clust = setdiff(1:Clust.num, del_clusts);
Clust.num = Clust.num-length(del_clusts);
Clust.num_comps = Clust.num_comps(keep_clust);
Clust.mode_num = Clust.mode_num(keep_clust);
Clust.comps = Clust.comps(keep_clust);
Clust.avg = Clust.avg(keep_clust);
Clust.norm_sq = Clust.norm_sq(keep_clust);
Clust.refined = Clust.refined(keep_clust);
Clust.changed = Clust.changed(keep_clust);
if Clust.found_offsets
    Clust.descender_off = Clust.descender_off(keep_clust);
    Clust.ascender_off = Clust.ascender_off(keep_clust);
end
if ~isempty(Clust.bigram)
    Clust.bigram = Clust.bigram(keep_clust,keep_clust);
end
for ii = 1:Clust.num
    Comps.clust(Clust.comps{ii}) = ii;
end


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%recalculate average: Use the listing of components passed to determine the
%average intensity image of them.
function new_avg = recalc_avg(Comps, idcs, use_thin, use_resize, rsz_method)
    imgs = get_comp_imgs(Comps, idcs);
    if use_resize && Comps.found_lines
        %normalize the images
        for ii=1:length(imgs)
            imgs{ii} = imresize(imgs{ii}, Comps.scale_factor(idcs(ii)), ...
                       rsz_method);
        end
    elseif use_resize
        warning('MBOCR:NoScaleFactor', ...
                'unable to renormalize images because scale_factor not set');
    end

    new_avg = imgs{1};
    new_tot = 1;
    for ii=2:length(imgs)
        [ht,wd] = size(imgs{ii});
        if all(size(new_avg) == [ht,wd])
            new_avg = (new_tot/(new_tot+1) .* new_avg) + ...
                      (1/(new_tot+1) .* imgs{ii});
        else
            new_avg = (new_tot/(new_tot+1) .* new_avg) + ...
                      (1/(new_tot+1) .* imresize(imgs{ii}, size(new_avg)));
        end
        new_tot = new_tot + 1;
    end

    if use_thin
        %thin the new image too
        new_avg = double(bwmorph(new_avg, 'thin', Inf));
    end
