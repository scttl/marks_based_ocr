function [Clust, Comps] = split_refine(Clust,Comps, dm, ms, thr)
% SPLIT_REFINE   Regroup inappropriately mergred clusters as separate clusters
%
%   [Clust, Comps] = split_refine(Clust, Comps, dist_metric, max_splits, thresh)
%
%   Clust should be struct containing several fields of a particular format.  
%   See cluster_comps for details.
%
%   Comps should be a struct containing component information.  Again see 
%   cluster_comps for details.
%
%   dist_metric is optional and if specified determines the type of distance
%   metric used for matching.  Valid options for this parameter are 'euc'
%   (straight Euclidian distance -- the default), 'conv_euc' (Euclidian distance
%   after convolving the matrices to find maximal overlapping point),
%   'hausdorff' to use Hausdorff distance, or 'ham' to use Hamming distance.
%
%   max_splits is optional and if specified determines the maximum number of
%   splits to perform on each cluster to check for max_splits + 1 piece matches.
%   This number should be low to reduce running-time, and defaults to 2 if not 
%   specified.
%
%   thresh is optional and if specified, determines the normalized maximal
%   Euclidian distance allowable for two cluster averages to be considered 
%   matching.  If not specified it defaults to a value of .01

% CVS INFO %
%%%%%%%%%%%%
% $Id: split_refine.m,v 1.3 2006-07-05 01:20:02 scottl Exp $
%
% REVISION HISTORY
% $Log: split_refine.m,v $
% Revision 1.3  2006-07-05 01:20:02  scottl
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
dist_thresh = .01;
max_splits  = 2;

bg_val = 0;  %the background pixel value in Comps
non_nb_val = 0;  %the value to use in Clust if a neighbour doesn't exist
min_width = inf;

%should we display matches onscreen (wastes resources)
display_matches = true;


% CODE START %
%%%%%%%%%%%%%%
if nargin < 2 || nargin > 5
    error('incorrect number of arguments specified!');
elseif nargin >= 3
    dist_metric = dm;
    if nargin >= 4
        max_splits = ms;
        if nargin == 5
            dist_thresh = thr;
        end
    end
end

%start by determining the pixel width of the smallest element in all the
%clusters (only split if the cluster width is at least twice this).
min_width = min(Comps.pos(:,3) - Comps.pos(:,1));

rr = find(Clust.refined == false, 1, 'first');
while ~isempty(rr)
    fprintf('                                                         \r');
    fprintf('cluster: %d  -- ', rr);

    %determine if this cluster can be refined by spliting it into up to 
    %split+1 pieces and matching each part with another cluster

    % first ensure that it is large enough to be split
    c_width = size(Clust.avg{rr}, 2);
    if (c_width >= 2 * min_width)
        idcs = [1:rr-1,rr+1:Clust.num];
        mres = find_match(max_splits, Clust.avg{rr}, Clust, idcs, min_width, ...
                          dist_metric, dist_thresh);
        if ~isempty(mres)
            %appropriate match found, reaverage items as appropriate
            fprintf('match found!\r');
            if display_matches
                num_matches = length(mres.clust) + 1;
                clf;
                subplot(2, num_matches, 1), imshow(Clust.avg{rr}), ...
                        xlabel('original');
            end

            while ~isempty(mres.clust)
                mc = mres.clust(1);
                if isempty(mres.sep_pos)
                    %right-most matching piece.  Use refined existing components
                    [Clust, Comps] = add_and_reaverage(Clust, Comps, ...
                                     mc, rr);
                else
                    %not right-most piece.  Create new components for the 
                    %split part and update positions, neighbours etc. of the 
                    %right part.
                    num_new_comps = Clust.num_comps(rr);
                    new_idcs = Comps.max_comp + [1:num_new_comps]';
                    old_idcs = Clust.comps{rr};
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
                    %component as a top, right, or bottom neighbour
                    %@@@to do!
                    %find(Comps.nb(:,3) == any of old_idcs)

                    Comps.offset(new_idcs) = Comps.offset(old_idcs) + ...
                           int16(Comps.pos(new_idcs,4)) - ...
                           int16(Comps.pos(old_idcs,4));

                    %add these new components to the appropriate cluster
                    num_prev_comps = Clust.num_comps(mc);
                    Clust.num_comps(mc) = Clust.num_comps(mc) + num_new_comps;
                    Clust.comps{mc} = [Clust.comps{mc}; new_idcs];
                    Clust.avg{mc} = (num_prev_comps/Clust.num_comps(mc) .* ...
                                    Clust.avg{mc}) + ...
                                    (num_new_comps/Clust.num_comps(mc) .* ...
                                    mres.avg{1});
                    Clust.norm_sq(mc) = sum(sum(Clust.avg{mc}.^2));
                end

                if display_matches
                    %add the new piece to the display
                    pos = num_matches - length(mres.clust) + 1;
                    subplot(2, num_matches, pos), imshow(mres.avg{1}), ...
                        xlabel(sprintf('piece %d\n', pos - 1));
                    subplot(2, num_matches, num_matches + pos), ...
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
                 dist_metric, thresh)
%this function attempts to recursively match halves using at most num_splits
%splits.  Match will be empty if no successful match could be found otherwise, 
%it will be a struct with fields: match.clust, match.sep_pos, and match.avg
%The first field will contain at least 1 scalar item specifying the 
%matching cluster index, the second will contain size(match.clust)-1 items 
%listing the positions of the separators used to find the matchings, and the
%final will be a cell array listing the intensity images were used in the
%matching of each piece

%@@@ implement calls to other distance metrics

bg_val = 0;  %value used for backgroundpixels

%attempt first to find a single valid match of all of C in Clust
D = euc_dist(C, Clust.avg(idcs), sum(sum(C.^2)), Clust.norm_sq(idcs));
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
        lrow = sort(lrow);
        lcol = sort(lcol);
        L = L(lrow(1):lrow(end),lcol(1):lcol(end));

        D = euc_dist(L, Clust.avg(idcs), sum(sum(L.^2)), Clust.norm_sq(idcs));
        [val, idx] = min(D);
        if val <= thresh
            %ensure we can find a right half match in up to one fewer splits
            %first trim blank space form R
            [row, col] = find(R ~= bg_val);
            row = sort(row);
            col = sort(col);
            R = R(row(1):row(end), col(1):col(end));
            mres = find_match(num_splits-1, R, Clust, idcs, min_width, ...
                              dist_metric, thresh);
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
