function [Clust, Comps] = merge_refine(Clust,Comps,dm,md,vp,mn)
% MERGE_REFINE  Merge together inappropriately separated clusters
%
%   [Clust, Comps, chg_list] = merge_refine(Clust, Comps, refine_list, ...
%                    dist_metric, max_dist, valid_pct, min_num)
%
%   Clust should be a struct containing cluster information.  See cluster_comps
%   for details on the format of each field.
%
%   Comps should be a struct containing information for each component.  See 
%   cluster_comps for details.
%
%   dist_metric is optional and if specified determines the type of distance
%   metric used for matching.  Valid options for this parameter are 'euc'
%   (straight Euclidian distance -- the default), 'conv' (Euclidian distance
%   after convolving the matrices to find maximal overlapping point), or
%   'hausdorff' to use Hausdorff distance.
%
%   max_dist is optional and if specified, determines the amount of tolerance 
%   (maximum number of pixels separating the closest parts of 2 components to
%   be merged).  If not specified, it defaults to 3 pixels.
%
%   valid_pct is optional and if speecified gives the percentage of neighbouring
%   components that must belong to the same cluster for the merge to be 
%   considered valid.  Note that this percentage must hold for both clusters, 
%   and is based on the total number of components in that cluster.  If not 
%   specified, it defaults to .85
%
%   min_num is optional and if specified determines the minimum number of
%   matching components that must be present to be considered for a merging. 
%   This alleviates merging clusters with single components that would otherwise
%   get matched.  If not specified, the defulat value is 3

% CVS INFO %
%%%%%%%%%%%%
% $Id: merge_refine.m,v 1.3 2006-07-05 01:21:23 scottl Exp $
%
% REVISION HISTORY
% $Log: merge_refine.m,v $
% Revision 1.3  2006-07-05 01:21:23  scottl
% started rewriting based on new Cluster and Component structures.  Not
% complete yet!
%
% Revision 1.2  2006/06/12 20:56:01  scottl
% changed comps to a cell array (indexed by page) to allow different sized pages
%
% Revision 1.1  2006/06/03 20:55:48  scottl
% Initial check-in.
%
%

% LOCAL VARS %
%%%%%%%%%%%%%%
dist_metric = 'euc';
max_dist = 3;
valid_pct = .85;
min_match_comp = 3;

bg_val = 0;   %pixel value for background colours

%should we display matches onscreen (wastes resources)
display_matches = false;


% CODE START %
%%%%%%%%%%%%%%

if nargin < 2 || nargin > 6
    error('incorrect number of arguments specified!');
elseif nargin >= 3
    dist_metric = dm;
    if nargin >= 4
        max_dist = md;
        if nargin >= 5
            valid_pct = vp;
            if nargin == 6
                min_match_comp = mn;
            end
        end
    end
end

num_clusts = length(Clust);

%initially put all elements in the refine list if not passed
rr = find(Clust.refined == false,1,'first');
while ~ isempty(rr)

    %determine if this item can be refined by merging it with a cluster
    %that almost always appears adjacent to it
    fprintf('                                                        \r');
    fprintf('cluster: %d -- ', r);
    lcl = []; tcl = [];
    if Clust.num_comps(rr) < min_match_comp
        %too few components in this cluster for a valid merge
        Clust.refined(rr) = true;
        rr = find(Clust.refined == false,1,'first');
        continue;
    end
    r_comps = Clust.comps(rr);

    %first look for a valid left neighbour cluster
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
        lnb_clusts = Comps.clust(lnb_comps);
        lnb_clust_list = unique(lnb_clusts);
        if length(lnb_clust_list) == 1
            mf_num = length(lnb_clusts);
            mf_clust = lnb_clust_list;
        else
            [mf_num, mf_clust] = max(hist(lnb_clusts, lnb_club_list));
            mf_clust = lnb_clust_list(mf_clst);
        end

        %ensure that the most frequent lnb cluster has the required minimum
        %number of components, makes up at least valid_pct of the total number 
        %of that clusters components as well as valid_pct of r's total number of
        %components.
        if mf_num >= min_match_comp && ...
           (mf_num / Clust.num_comps(mf_clust)) >= valid_pct && ...
           (mf_num / Clust.num_comps(rr)) >= valid_pct
            %valid merge found, update neighbours, component averages etc.
            fprintf('valid left merge found\r');
            if display_matches
                clf;
                subplot(1,2,1), imshow(Clust.avg{mf_clust});
                xlabel(Clust.num_comps(mf_clst));
                subplot(1,2,2), imshow(Clustavg{rr});
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
            %for top and bottom, take the closer of the two possible neighbours.
            %If there is a tie, take the maximally overlapping neighbour
            [Dummy, idx] = 

            %clusters -- we'll create a new cluster below
            Comps.clust(r_merge_comps) = num_clusts+1;

            %offsets

            %remove the now merged left-neighbour components


            num_clusts = num_clusts + 1;
            keep_list = [keep_list, num_clusts];
            chg_list = [chg_list, num_clusts];

                        %flag this item for removal from r
                        pos = find(r_keep == i,1);
                        r_keep = [r_keep(1:pos-1), r_keep(pos+1:end)];

                        %remove the item from the mfc
                        pos = find(mf_keep == idx,1);
                        mf_keep = [mf_keep(1:pos-1), mf_keep(pos+1:end)];
                    end
                end
            end


            %delete the necessary components from r and mfc
            if r == mf_clst
                %if we are merging items from the same cluster, we have to
                %update our keep_lists based on what is common to both only
                r_keep = intersect(r_keep, mf_keep);
            end
            if length(r_keep) == 0
                r_pos = find(keep_list == r,1);
                keep_list = [keep_list(1:r_pos-1), keep_list(r_pos+1:end)];
            else
                Clust(r).comp = Clust(r).comp(r_keep);
                Clust(r).pos = Clust(r).pos(r_keep,:);
                Clust(r).nb = Clust(r).nb(r_keep,:);
                Clust(r).pg = Clust(r).pg(r_keep);
                Clust(r).num = length(r_keep);
                chg_list = [chg_list, r];
            end
            if r ~= mf_clst
                if length(mf_keep) == 0
                    mf_pos = find(keep_list == mf_clst,1);
                    keep_list = [keep_list(1:mf_pos-1),keep_list(mf_pos+1:end)];
                else
                    Clust(mf_clst).comp = Clust(mf_clst).comp(mf_keep);
                    Clust(mf_clst).pos = Clust(mf_clst).pos(mf_keep,:);
                    Clust(mf_clst).nb = Clust(mf_clst).nb(mf_keep,:);
                    Clust(mf_clst).pg = Clust(mf_clst).pg(mf_keep);
                    Clust(mf_clst).num = length(mf_keep);
                    chg_list = [chg_list, mf_clst];
                end
            end

            refine_list = refine_list(2:end);
            continue;  %don't try and find a top match too
        end
    end

    %look for a top match
    for i=1:Clust(r).num
        tnb = Clust(r).nb(i,2);
        if tnb ~= 0
            [cl,off] = get_cl_off(Clust, Clust(r).nb(i,2));
            if (Clust(r).pos(i,2) - Clust(cl).pos(off,4)) > pixel_thresh
                continue;  %match distance too large
             else
                %ensure that there are 'on' pixels within the distance
                top = min(Clust(r).pos(i,2), Clust(cl).pos(off,4)) - ...
                       pixel_thresh + 1;
                bot = max(Clust(r).pos(i,2), Clust(cl).pos(off,4)) + ...
                       pixel_thresh - 1;
                left = min(Clust(cl).pos(off,1), Clust(r).pos(i,1));
                right = max(Clust(cl).pos(off,3), Clust(r).pos(i,3));
                if top < 1
                    top = 1;
                end
                if bot > size(Comps{Clust(r).pg(i)},1)
                    bot = size(Comps{Clust(r).pg(i)},1);
                end
                [tr, tc] = find(Comps{Clust(cl).pg(off)}(top:bot,left:right) ...
                           == Clust(cl).comp(off));
                [br, bc] = find(Comps{Clust(r).pg(i)}(top:bot, left:right) ...
                           == Clust(r).comp(i));
                match_found = false;
                for j=1:length(tr)
                    for k=1:length(br)
                        if abs(br(k) - tr(j)) + abs(bc(k) - tc(j)) <= ...
                           pixel_thresh
                            tcl = [tcl, [cl; off]];
                            match_found = true;
                            break;
                        end
                    end
                    if match_found
                        break;
                    end
                end
            end
        end
    end
    if ~isempty(tcl)
        %determine the most frequently occuring cluster in tcl
        [mf_cnt, mf_clst] = max(hist(tcl(1,:), max(tcl(1,:)) - ...
                            min(tcl(1,:)) + 1));
        mf_clst = mf_clst + min(tcl(1,:)) - 1;
        if mf_cnt >= min_match_comp && (mf_cnt / Clust(r).num) >= valid_pct ...
           && (mf_cnt / Clust(mf_clst).num) >= valid_pct
            fprintf('valid top merge found\r');
            if display_matches
                clf;
                subplot(2,1,1), imshow(Clust(mf_clst).avg), ...
                    xlabel(Clust(mf_clst).num);
                subplot(2,1,2), imshow(Clust(r).avg), ...
                    xlabel(Clust(r).num);
                drawnow;
                pause(2);
            end

            %add the new cluster to hold the merged items
            num_clusts = num_clusts + 1;
            keep_list = [keep_list, num_clusts];
            chg_list = [chg_list, num_clusts];
            Clust(num_clusts).num = 0;

            r_keep = 1:Clust(r).num;
            mf_keep = 1:Clust(mf_clst).num;
            for i=1:Clust(r).num
                tnb = Clust(r).nb(i,2);
                if tnb ~= 0
                    idx = find(Clust(mf_clst).comp == tnb, 1);
                    if ~isempty(idx)
                        this.comp = Comps.max_comp + 1;
                        Comps.max_comp = Comps.max_comp + 1;
                        mfc_pos = Clust(mf_clst).pos(idx,:);
                        r_pos = Clust(r).pos(i,:);
                        ll = min(mfc_pos(1), r_pos(1)); tt = mfc_pos(2);
                        rr = max(mfc_pos(3), r_pos(3)); bb = r_pos(4);
                        this.pos = [ll, tt, rr, bb];
                        
                        %update the new merge components neighbours
                        mfc_nb = Clust(mf_clst).nb(idx,:);
                        r_nb = Clust(r).nb(i,:);
                        if mfc_nb(1) == 0 || mfc_nb(1) == r_nb(1) || ...
                           mfc_nb(1) == Clust(r).comp(i)
                            lnb = r_nb(1);
                        elseif r_nb(1) == 0 || ... 
                           r_nb(1) == Clust(mf_clst).comp(idx)
                            lnb = mfc_nb(1);
                        else
                            %take the maximally overlapping match of the two
                            [cl,off] = get_cl_off(Clust, mfc_nb(1));
                            pos = Clust(cl).pos(off,:);
                            mfc_ovrlp = min(bb, pos(4)) - max(tt, pos(2));
                            [cl,off] = get_cl_off(Clust, r_nb(1));
                            pos = Clust(cl).pos(off,:);
                            r_ovrlp = min(bb, pos(4)) - max(tt, pos(2));
                            if r_ovrlp >= mfc_ovrlp
                                lnb = r_nb(1);
                            else
                                lnb = mfc_nb(1);
                            end
                        end
                        if mfc_nb(3) == 0 || mfc_nb(3) == r_nb(3) || ...
                           mfc_nb(3) == Clust(r).comp(i)
                            rnb = r_nb(3);
                        elseif r_nb(4) == 0 || ...
                           r_nb(3) == Clust(mf_clst).comp(idx)
                            rnb = mfc_nb(3);
                        else
                            [cl,off] = get_cl_off(Clust, mfc_nb(3));
                            pos = Clust(cl).pos(off,:);
                            mfc_ovrlp = min(bb, pos(4)) - max(tt, pos(2));
                            [cl,off] = get_cl_off(Clust, r_nb(3));
                            pos = Clust(cl).pos(off,:);
                            r_ovrlp = min(bb, pos(4)) - max(tt, pos(2));
                            if r_ovrlp >= mfc_ovrlp
                                rnb = r_nb(3);
                            else
                                rnb = mfc_nb(3);
                            end
                        end
                        this.nb = [lnb, mfc_nb(2), rnb, r_nb(4)];

                        this.pg = Clust(r).pg(i);
                        this.num = 1;

                        %update the cluster average
                        this.avg = zeros(bb-tt+1, rr-ll+1);
                        this.avg(find(Comps{this.pg}(tt:bb,ll:rr) ~= ...
                               bg_val)) = 1;

                        if Clust(num_clusts).num == 0
                            Clust(num_clusts) = this;
                        else
                            if any(size(this.avg)~=size(Clust(num_clusts).avg))
                                %resize this.avg to be the same size
                                this.avg = imresize(this.avg, ...
                                           size(Clust(num_clusts).avg));
                            end
                            Clust(num_clusts) = add_and_reaverage(...
                            Clust(num_clusts), this, Clust(num_clusts).avg, ...
                            this.avg);
                        end

                        %update Comps appropriately
                        region = Comps{this.pg}(tt:bb, ll:rr);
                        chg_idx = find(region == Clust(r).comp(i) | ...
                                       region == Clust(mf_clst).comp(idx));
                        region(chg_idx) = this.comp;
                        Comps{this.pg}(tt:bb, ll:rr) = region;

                        %update neighbours of other elements that point at
                        %either of these components
                        for j=keep_list
                            chg_idx = find(Clust(j).nb == Clust(r).comp(i) | ...
                                      Clust(j).nb == Clust(mf_clst).comp(idx));
                            if ~isempty(chg_idx)
                                Clust(j).nb(chg_idx) = this.comp;
                            end
                        end

                        %flag this item for removal from r
                        pos = find(r_keep == i,1);
                        r_keep = [r_keep(1:pos-1), r_keep(pos+1:end)];

                        %remove the item from the mfc
                        pos = find(mf_keep == idx,1);
                        mf_keep = [mf_keep(1:pos-1), mf_keep(pos+1:end)];
                    end
                end
            end

            %delete the necessary components from r and mfc
            if r == mf_clst
                r_keep = intersect(r_keep, mf_keep);
            end
            if length(r_keep) == 0
                r_pos = find(keep_list == r,1);
                keep_list = [keep_list(1:r_pos-1), keep_list(r_pos+1:end)];
            else
                Clust(r).comp = Clust(r).comp(r_keep);
                Clust(r).pos = Clust(r).pos(r_keep,:);
                Clust(r).nb = Clust(r).nb(r_keep,:);
                Clust(r).pg = Clust(r).pg(r_keep);
                Clust(r).num = length(r_keep);
                chg_list = [chg_list, r];
            end
            if r ~= mf_clst
                if length(mf_keep) == 0
                    mf_pos = find(keep_list == mf_clst,1);
                    keep_list = [keep_list(1:mf_pos-1),keep_list(mf_pos+1:end)];
                else
                    Clust(mf_clst).comp = Clust(mf_clst).comp(mf_keep);
                    Clust(mf_clst).pos = Clust(mf_clst).pos(mf_keep,:);
                    Clust(mf_clst).nb = Clust(mf_clst).nb(mf_keep,:);
                    Clust(mf_clst).pg = Clust(mf_clst).pg(mf_keep);
                    Clust(mf_clst).num = length(mf_keep);
                    chg_list = [chg_list, mf_clst];
                end
            end
        end
    end
    refine_list = refine_list(2:end);
end
fprintf('\n');

[Dummy, Dummy, chg_list] = intersect(unique(chg_list), keep_list);
Clust = Clust(keep_list);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
