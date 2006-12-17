function [Clust, Comps] = split_refine(Clust, Comps, varargin)
% SPLIT_REFINE   Regroup inappropriately mergred clusters as separate clusters
%
%   [CLUST, COMPS] = split_refine(CLUST, COMPS, [VAR1, VAL1]...)
%
%   CLUST should be struct containing several fields of a particular format.  
%   See cluster_comps for details.
%
%   COMPS should be a struct containing component information.  Again see 
%   get_comps for details.
%
%   default values for parameters defined in LOCAL VARS below can be overridden
%   by specifying the name as a string in VAR1, and the new value in VAL1.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: split_refine.m,v 1.10 2006-12-17 20:17:08 scottl Exp $
%
% REVISION HISTORY
% $Log: split_refine.m,v $
% Revision 1.10  2006-12-17 20:17:08  scottl
% bugfix in checking when value is NaN
%
% Revision 1.9  2006/10/18 15:49:32  scottl
% changes to allow scale-invariant splits, better processing of optional
% arguments, some efficiency improvements.
%
% Revision 1.8  2006/08/24 21:40:04  scottl
% added ability to use the mode instead of taking the average of cluster
% intensities while refining.
%
% Revision 1.7  2006/08/14 01:35:14  scottl
% Updates based on new Clust.changed field.
%
% Revision 1.6  2006/08/07 21:21:01  scottl
% fixup small display bug, as well as a larger incorrect neighbour updating bug.
%
% Revision 1.5  2006/08/06 22:08:24  scottl
% bugfix to handle split matches of different size when re-averaging.
%
% Revision 1.4  2006/07/22 04:12:58  scottl
% implemented other distance metrics, cleaned up display bug
%
% Revision 1.3  2006/07/05 01:20:02  scottl
% rewritten based on new Cluster and Component structures.
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
dist_thresh = .013;
max_splits  = 2; %refers to the number of cuts --> max_splits + 1 pieces

bg_val = 0;  %the background pixel value in Comps
min_width = inf;

%should we use the average or the mode when combining clusters?
use_avg = true;

%by default we don't rescale components, unless the boolean below is set to
%true.  This also requires that Comps.scale_factor exists (and thus that
%get_lines() has been run)
resize_imgs = false;
%which method should be used to resize components (see imresize() for choices)
resize_method = 'nearest';  

%by default we don't work with thinned representations of the images
use_thinned_imgs = false;

%should we display matches onscreen (wastes resources)
display_matches = false;


% CODE START %
%%%%%%%%%%%%%%
if nargin < 2
    error('incorrect number of arguments specified!');
elseif nargin > 2
    process_optional_args(varargin{:});
end

%start by determining the pixel width of the smallest element in all the
%clusters (only split if the cluster width is at least twice this).
min_width = min(Comps.pos(:,3) - Comps.pos(:,1) + 1);

modal_height = NaN;
if resize_imgs && Comps.found_lines
    modal_height = Comps.modal_height;
end

rr = find(Clust.refined == false, 1, 'first');
while ~isempty(rr)
    fprintf('                                                         \r');
    fprintf('cluster: %d  -- ', rr);

    %determine if this cluster can be refined by spliting it into up to 
    %split+1 pieces and matching each part with another cluster

    %first ensure that it is large enough to be split
    c_width = size(Clust.avg{rr}, 2);
    if resize_imgs && Comps.found_lines
        %since averages have been rescaled, look at typical component width
        idx = Clust.comps{rr};
        c_width = mode(double(Comps.pos(idx,3) - Comps.pos(idx,1) + 1));
    end
    if (c_width >= 2 * min_width)
        idcs = [1:rr-1,rr+1:Clust.num];
        mres = find_match(max_splits, Clust.avg{rr}, Clust, idcs, min_width, ...
                          dist_metric, dist_thresh, modal_height, ...
                          resize_method, use_thinned_imgs);
        if ~isempty(mres)
            %appropriate match found, reaverage items as appropriate
            fprintf('match found!\r');
            num_matches = length(mres.clust);
            if display_matches
                clf;
                subplot(2, num_matches+1, 1), imshow(Clust.avg{rr}), ...
                        xlabel('original');
            end

            num_new_comps = Clust.num_comps(rr);
            old_idcs = Clust.comps{rr};
            %these matrices list the indices of the neighbour components, their
            %components, and the median L-R position of those components that
            %list one of the split comoponents as a top, right,or bot. neighbour
            r_nbs = [];
            t_nbs = [];
            b_nbs = [];
            for ii= 1:length(old_idcs)
                this_rnbs = find(Comps.nb(:,3) == old_idcs(ii));
                r_nbs = [r_nbs; [repmat(ii,length(this_rnbs),1), this_rnbs]];
                this_tnbs = find(Comps.nb(:,2) == old_idcs(ii));
                t_nbs = [t_nbs; [repmat(ii,length(this_tnbs),1), this_tnbs, ...
                            Comps.pos(this_tnbs,1) + ...
                            (floor(diff(Comps.pos(this_tnbs,[1,3]),1,2) ./2))]];
                this_bnbs = find(Comps.nb(:,4) == old_idcs(ii));
                b_nbs = [b_nbs; [repmat(ii,length(this_bnbs),1), this_bnbs, ...
                            Comps.pos(this_bnbs,1) + ...
                            (floor(diff(Comps.pos(this_bnbs,[1,3]),1,2) ./2))]];
            end

            while ~isempty(mres.clust)
                mc = mres.clust(1);
                if isempty(mres.sep_pos)
                    %right-most matching piece.  Use refined existing components
                    [Clust, Comps, mc] = add_and_reaverage(Clust, Comps, ...
                                         mc, rr, 'use_avg', use_avg);
                    Clust.changed(mc) = true;
                else
                    %not right-most piece.  Create new components for the 
                    %split part and update positions, neighbours etc. of the 
                    %right part.
                    new_idcs = Comps.max_comp + (1:num_new_comps)';
                    Comps.max_comp = Comps.max_comp + num_new_comps;
                    Comps.clust(new_idcs) = mc;
                    Comps.pos(new_idcs,:) = Comps.pos(old_idcs,:);
                    rem_pos = Comps.pos(old_idcs,1) + mres.sep_pos(1) - 1;
                    Comps.pos(new_idcs,3) = rem_pos;
                    Comps.pos(old_idcs,1) = rem_pos + 1;
                    %potentially update the top and bottom position too
                    Comps.pos(new_idcs,2) = Comps.pos(new_idcs,2)+mres.top(1)-1;
                    Comps.pos(new_idcs,4) = Comps.pos(new_idcs,4) - ...
                                           (size(mres.avg{1},1) - mres.bot(1));
                    Comps.pg(new_idcs) = Comps.pg(old_idcs);
                    Comps.nb(new_idcs,:) = Comps.nb(old_idcs,:);
                    Comps.nb(new_idcs,3) = old_idcs;
                    Comps.nb(old_idcs,1) = new_idcs;
                    Comps.nb_dist(new_idcs,:) = Comps.nb_dist(old_idcs,:);
                    Comps.nb_dist(new_idcs,3) = 1;
                    Comps.nb_dist(old_idcs,1) = 1;

                    %potentially update any components that list the old
                    %component as a top, right, or bottom neighbour.
                    %Neighbours listing the split component as a right neighbour
                    %should be updated to point at the left-most i.e first
                    %piece of the split component
                    if length(mres.clust) == num_matches
                        %left-most piece
                        Comps.nb(r_nbs(:,2),3) = new_idcs(r_nbs(:,1));
                    end
                    %update top and bottom neighbours whos median L-R position
                    %is to the left of the right position of the split
                    %component piece
                    top_idx = find(t_nbs(:,3)<=Comps.pos(new_idcs(t_nbs(:,1))));
                    Comps.nb(t_nbs(top_idx,2),2) = new_idcs(t_nbs(top_idx,1));
                    t_nbs = t_nbs(setdiff((1:size(t_nbs,1))',top_idx),:);
                    bot_idx = find(b_nbs(:,3)<=Comps.pos(new_idcs(b_nbs(:,1))));
                    Comps.nb(b_nbs(bot_idx,2),4) = new_idcs(b_nbs(bot_idx,1));
                    b_nbs = b_nbs(setdiff((1:size(b_nbs,1))',bot_idx),:);

                    if Comps.found_lines
                        Comps.line(new_idcs) = Comps.line(old_idcs);
                        Comps.descender_off(new_idcs) =  ...
                               Comps.descender_off(old_idcs) + ...
                               int16(Comps.pos(new_idcs,4)) - ...
                               int16(Comps.pos(old_idcs,4));
                        Comps.ascender_off(new_idcs) =  ...
                               Comps.ascender_off(old_idcs) + ...
                               int16(Comps.pos(new_idcs,2)) - ...
                               int16(Comps.pos(old_idcs,2));
                        Comps.scale_factor(new_idcs) = ...
                               Comps.scale_factor(old_idcs);
                    end

                    if Comps.found_true_labels
                        Comps.truth_label(new_idcs) = ...
                               Comps.truth_label(old_idcs);
                    end

                    %add these new components to the appropriate cluster
                    num_prev_comps = Clust.num_comps(mc);
                    Clust.num_comps(mc) = Clust.num_comps(mc) + num_new_comps;
                    Clust.comps{mc} = [Clust.comps{mc}; new_idcs];

                    if use_avg
                        Clust.mode_num(mc) = Clust.mode_num(mc) + num_new_comps;
                        %note that the size of the matching averages may differ
                        %slightly, so we resize the match if required
                        if any(size(Clust.avg{mc}) ~= size(mres.avg{1}))
                            mres.avg{1} = imresize(mres.avg{1}, ...
                                          size(Clust.avg{mc}));
                        end
                        Clust.avg{mc} = (num_prev_comps/Clust.num_comps(mc) ...
                                        .* Clust.avg{mc}) + ...
                                        (num_new_comps/Clust.num_comps(mc) ...
                                        .* mres.avg{1});
                    else
                        %take the mode
                        idcs = [mc; rr];
                        [md_val,md_idx] = max(Clust.mode_num(idcs));
                        Clust.avg{mc} = Clust.avg{idcs(md_idx)};
                        Clust.mode_num(mc) = md_val;
                    end

                    Clust.norm_sq(mc) = sum(Clust.avg{mc}(:).^2);
                    Clust.changed(mc) = true;

                    if length(mres.sep_pos) == 1
                        %slice off the now removed part of the average from the 
                        %old cluster, and trim any blankspace
                        Clust.avg{rr} = Clust.avg{rr}(:,mres.sep_pos(1)+1:end);
                        [row, col] = find(Clust.avg{rr} ~= bg_val);
                        row = sort(row); col = sort(col);
                        Clust.avg{rr} = Clust.avg{rr}(row(1):row(end),...
                                                      col(1):col(end));
                        Clust.norm_sq(rr) = sum(Clust.avg{rr}(:).^2);
                    end
                end

                if display_matches
                    %add the new piece to the display
                    pos = num_matches+1 - (length(mres.clust)-1);
                    subplot(2, num_matches+1, pos), imshow(mres.avg{1}), ...
                        xlabel(sprintf('piece %d\n', pos - 1));
                    subplot(2, num_matches+1, num_matches+pos+1), ...
                        imshow(Clust.avg{mc}), xlabel(...
                        sprintf('matching cluster %d\n', mc));
                end
                mres.clust = mres.clust(2:end);
                mres.sep_pos = mres.sep_pos(2:end);
                mres.avg = mres.avg(2:end);
                mres.top = mres.top(2:end);
                mres.bot = mres.bot(2:end);
            end

            if display_matches
                drawnow;
                pause(2);
            end
        else
            Clust.refined(rr) = true;
            fprintf('no match\r');
        end
    else
        Clust.refined(rr) = true;
        fprintf('too narrow\r');
    end
    rr = find(Clust.refined == false, 1, 'first');
end



% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function match = find_match(num_splits, C, Clust, idcs, min_width, ...
                 dist_metric, thresh, modal_height, resize_method, thin_img)
%this function attempts to recursively match halves using at most num_splits
%splits.  Match will be empty if no successful match could be found otherwise, 
%it will be a struct with fields: match.clust, match.sep_pos, and match.avg
%The first field will contain at least 1 scalar item specifying the 
%matching cluster index, the second will contain size(match.clust)-1 items 
%listing the positions of the separators used to find the matchings, and the
%final will be a cell array listing the intensity images were used in the
%matching of each piece

bg_val = 0;  %value used for background pixels

%attempt first to find a single valid match of all of C in Clust
if strcmp(dist_metric, 'euc')
    D = euc_dist(C, Clust.avg(idcs), sum(C(:).^2), Clust.norm_sq(idcs));
elseif strcmp(dist_metric, 'hausdorff')
    D = hausdorff_dist(C, Clust.avg(idcs));
elseif strcmp(dist_metric, 'ham')
    D = ham_dist(C, Clust.avg(idcs));
elseif strcmp(dist_metric, 'conv_euc')
    D = conv_euc_dist(C, Clust.avg(idcs), sum(C(:).^2), Clust.norm_sq(idcs));
else
    error('incorrect distance metric specified!');
end

[val,idx] = min(D);
if val <= thresh
    %match found
    match.clust = idcs(idx);
    match.sep_pos = [];
    match.avg = {C};
    match.top = 1;
    match.bot = size(C,1);
    return;
end

if num_splits >= 1
    %recursive case, first look for a left-side match
    c_width = size(C,2);
    cut_pos = min_width;
    while(cut_pos <= (c_width - min_width))
        L = C(:,1:cut_pos);
        R = C(:,cut_pos+1:end);

        %first trim any blank space from L to create a tight bounding box
        [lrow, lcol] = find(L ~= bg_val);
        if isempty(lrow) && isempty(lcol)
            %left side is completely empty, just move on to next cut position
            warning('MBOCR:LeftEmpty', 'left side completely empty');
            cut_pos = cut_pos + 1;
            continue;
        end
        lrow = sort(lrow);
        lcol = sort(lcol);
        L = L(lrow(1):lrow(end),lcol(1):lcol(end));

        if ~isnan(modal_height)
            L = imresize(L, modal_height ./ size(L,1), resize_method);
        end
        if thin_img
            L = double(bwmorph(L, 'thin', Inf));
        end

        if strcmp(dist_metric, 'euc')
            D = euc_dist(L, Clust.avg(idcs), sum(L(:).^2), Clust.norm_sq(idcs));
        elseif strcmp(dist_metric, 'hausdorff')
            D = hausdorff_dist(L, Clust.avg(idcs));
        elseif strcmp(dist_metric, 'ham')
            D = ham_dist(L, Clust.avg(idcs));
        elseif strcmp(dist_metric, 'conv_euc')
            D = conv_euc_dist(L, Clust.avg(idcs), sum(L(:).^2), ...
                Clust.norm_sq(idcs));
        else
            error('incorrect distance metric specified!');
        end
        [val, idx] = min(D);
        if val <= thresh
            %ensure we can find a right half match in up to one fewer splits
            %first trim blank space form R
            [row, col] = find(R ~= bg_val);
            row = sort(row);
            col = sort(col);
            R = R(row(1):row(end), col(1):col(end));
            mres = find_match(num_splits-1, R, Clust, idcs, min_width, ...
                   dist_metric, thresh, modal_height, resize_method, thin_img);
            if ~isempty(mres)
                %right match found!
                match.clust = [idcs(idx), mres.clust];
                match.sep_pos = [cut_pos, mres.sep_pos];
                match.avg = {L, mres.avg{:}};
                if length(mres.top) == 1
                    match.top = [lrow(1), row(1)];
                    match.bot = [lrow(end), row(end)];
                else
                    match.top = [lrow(1), mres.top];
                    match.bot = [lrow(end), mres.bot];
                end
                return;
            end
        end
        cut_pos = cut_pos + 1;
    end
end
%if we get down here, a match hasn't been found, set match to an empty array
match = [];
