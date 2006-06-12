function [Clust, Comps, chg_list] = merge_refine(Clust,Comps,rl,dm,md,vp,mn)
% MERGE_REFINE  Merge together inappropriately separated clusters
%
%   [Clust, Comps, chg_list] = merge_refine(Clust, Comps, refine_list, ...
%                    dist_metric, max_dist, valid_pct, min_num)
%
%   Clust should be an array of structs, each of which is assumed to contain 
%   several fields of a particular format.  See cluster_comps for details.
%
%   Comps should be a cell array indexed by page, where each entry of each entry
%   represents a pixel, and each 'on' pixel is labelled with the component 
%   number to which it belongs.
%
%   refine_list is optional and if specified, determines which items are
%   to be refined (ie searched for a match).  It should be a vector of cluster
%   indices and will be processed in the order given.  If not specified, this 
%   method will attempt to refine all clusters, starting from the highest 
%   numbered one.
%
%   dist_metric is optional and if specified determines the type of distance
%   metric used for matching.  Valid options for this parameter are 'euc'
%   (straight Euclidian distance -- the default), 'conv' (Euclidian distance
%   after convolving the matrices to find maximal overlapping point), or
%   'hausdorff' to use Hausdorff distance.  NOTE: for mergers, the distance
%   metric specified is ignored (only used for consistency with other refine
%   calls).
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
% $Id: merge_refine.m,v 1.2 2006-06-12 20:56:01 scottl Exp $
%
% REVISION HISTORY
% $Log: merge_refine.m,v $
% Revision 1.2  2006-06-12 20:56:01  scottl
% changed comps to a cell array (indexed by page) to allow different sized pages
%
% Revision 1.1  2006/06/03 20:55:48  scottl
% Initial check-in.
%
%

% LOCAL VARS %
%%%%%%%%%%%%%%
dist_metric = 'euc';
pixel_thresh = 3;
valid_pct = .85;
min_match_comp = 3;

bg_val = 0;   %pixel value for background colours

%should we display matches onscreen (wastes resources)
display_matches = false;


% CODE START %
%%%%%%%%%%%%%%

if nargin < 2 || nargin > 7
    error('incorrect number of arguments specified!');
elseif nargin >= 3
    refine_list = rl;
    if nargin >= 4
        dist_metric = dm;
        if nargin >= 5
            pixel_thresh = md;
            if nargin >= 6
                valid_pct = vp;
                if nargin == 7
                    min_match_comp = mn;
                end
            end
        end
    end
end

chg_list = [];
num_clusts = size(Clust, 1);
max_comp  = max(max(Comps{end}));

%initially put all elements in the refine list if not passed
if nargin < 3
    refine_list = num_clusts:-1:1;
end
keep_list = 1:num_clusts;

while ~ isempty(refine_list)

    %determine if the first item in the list can be refined by merging it
    %with a neighbouring component that is always the same cluster
    r = refine_list(1);
    fprintf('                                                        \r');
    fprintf('cluster: %d -- ', r);
    lcl = []; tcl = [];
    if Clust(r).num < min_match_comp
        %too few components in this cluster for a valid merge
        refine_list = refine_list(2:end);
        continue;
    end

    for i=1:Clust(r).num
        lnb = Clust(r).nb(i,1);
        if lnb ~= 0
            [cl,off] = get_cl_off(Clust, Clust(r).nb(i,1));
            if (Clust(r).pos(i,1) - Clust(cl).pos(off,3)) > pixel_thresh
                continue;  %match distance too large
            else
                %ensure that there are 'on' pixels within the distance
                top = min(Clust(r).pos(i,2), Clust(cl).pos(off,2));
                bot = max(Clust(r).pos(i,4), Clust(cl).pos(off,4));
                left = min(Clust(cl).pos(off,3), Clust(r).pos(i,1)) - ...
                       pixel_thresh + 1;
                right = max(Clust(cl).pos(off,3), Clust(r).pos(i,1)) + ...
                       pixel_thresh - 1;
                if left < 1
                    left = 1;
                end
                if right > size(Comps{Clust(r).pg(i)},2)
                    right = size(Comps{Clust(r).pg(i)},2);
                end
                [lr, lc] = find(Comps{Clust(cl).pg(off)}(top:bot,left:right) ...
                           == Clust(cl).comp(off));
                [rr, rc] = find(Comps{Clust(r).pg(i)}(top:bot,left:right) ...
                           == Clust(r).comp(i));
                match_found = false;
                for j=1:length(lr)
                    for k=1:length(rr)
                        if abs(rr(k) - lr(j)) + abs(rc(k) - lc(j)) <= ...
                           pixel_thresh
                            lcl = [lcl, [cl; off]];
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
    if ~isempty(lcl)
        % determine the most frequently occuring cluster in lcl
        [mf_cnt, mf_clst] = max(hist(lcl(1,:), max(lcl(1,:)) - ...
                            min(lcl(1,:)) + 1));
        mf_clst = mf_clst + min(lcl(1,:)) - 1;
        if mf_cnt >= min_match_comp && (mf_cnt / Clust(r).num) >= valid_pct ...
           && (mf_cnt / Clust(mf_clst).num) >= valid_pct
            %valid merge found, update neighbours, component averages etc.
            fprintf('valid left merge found\r');
            if display_matches
                clf;
                subplot(1,2,1), imshow(Clust(mf_clst).avg), ...
                    xlabel(Clust(mf_clst).num);
                subplot(1,2,2), imshow(Clust(r).avg), ...
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
                lnb = Clust(r).nb(i,1);
                if lnb ~= 0
                    idx = find(Clust(mf_clst).comp == lnb, 1);
                    if ~isempty(idx)
                        this.comp = max_comp + 1;
                        max_comp = max_comp + 1;
                        mfc_pos = Clust(mf_clst).pos(idx,:);
                        r_pos = Clust(r).pos(i,:);
                        ll = mfc_pos(1); tt = min(mfc_pos(2), r_pos(2));
                        rr = r_pos(3); bb = max(mfc_pos(4), r_pos(4));
                        this.pos = [ll, tt, rr, bb];
                        
                        %update the new merge components neighbours
                        mfc_nb = Clust(mf_clst).nb(idx,:);
                        r_nb = Clust(r).nb(i,:);
                        if mfc_nb(2) == 0 || mfc_nb(2) == r_nb(2) || ...
                           mfc_nb(2) == Clust(r).comp(i)
                            tnb = r_nb(2);
                        elseif r_nb(2) == 0 || ...
                           r_nb(2) == Clust(mf_clst).comp(idx)
                            tnb = mfc_nb(2);
                        else
                            %take the maximally overlapping match of the two
                            %@@ this may not be optimal: --- -- ---
                            %                             llllrrrr
                            [cl,off] = get_cl_off(Clust, mfc_nb(2));
                            pos = Clust(cl).pos(off,:);
                            mfc_ovrlp = min(rr, pos(3)) - max(ll, pos(1));
                            [cl,off] = get_cl_off(Clust, r_nb(2));
                            pos = Clust(cl).pos(off,:);
                            r_ovrlp = min(rr, pos(3)) - max(ll, pos(1));
                            if r_ovrlp >= mfc_ovrlp
                                tnb = r_nb(2);
                            else
                                tnb = mfc_nb(2);
                            end
                        end
                        if mfc_nb(4) == 0 || mfc_nb(4) == r_nb(4) || ...
                           mfc_nb(4) == Clust(r).comp(i)
                            bnb = r_nb(4);
                        elseif r_nb(4) == 0 || ...
                           r_nb(4) == Clust(mf_clst).comp(idx)
                            bnb = mfc_nb(4);
                        else
                            [cl,off] = get_cl_off(Clust, mfc_nb(4));
                            pos = Clust(cl).pos(off,:);
                            mfc_ovrlp = min(rr, pos(3)) - max(ll, pos(1));
                            [cl,off] = get_cl_off(Clust, r_nb(4));
                            pos = Clust(cl).pos(off,:);
                            r_ovrlp = min(rr, pos(3)) - max(ll, pos(1));
                            if r_ovrlp >= mfc_ovrlp
                                bnb = r_nb(4);
                            else
                                bnb = mfc_nb(4);
                            end
                        end
                        this.nb = [mfc_nb(1), tnb, r_nb(3), bnb];

                        this.pg = Clust(r).pg(i);
                        this.num = 1;

                        %update the cluster average
                        this.avg = zeros(bb-tt+1, rr-ll+1);
                        this.avg(find(Comps{this.pg}(tt:bb,ll:rr) ~= ...
                               bg_val)) = 1;

                        %update the cluster
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
                        this.comp = max_comp + 1;
                        max_comp = max_comp + 1;
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
