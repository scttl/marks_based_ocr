function [Clust, Comps, chg_list] = split_refine(Clust,Comps,rl,dm,ms,thr)
% SPLIT_REFINE   Regroup inappropriately mergred clusters as separate clusters
%
%   [Clust, Comps, chg_list] = split_refine(Clust, Comps, refine_list, ...
%                              dist_metric, max_splits, thresh)
%
%   Clust should be an array of structs, each of which is assumed to contain 
%   several fields of a particular format.  See cluster_comps for details.
%
%   Comps should be a 2 or 3 dimensional image matrix, where each entry 
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
% $Id: split_refine.m,v 1.1 2006-06-03 20:55:48 scottl Exp $
%
% REVISION HISTORY
% $Log: split_refine.m,v $
% Revision 1.1  2006-06-03 20:55:48  scottl
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
if nargin < 2 || nargin > 6
    error('incorrect number of arguments specified!');
elseif nargin >= 3
    refine_list = rl;
    if nargin >= 4
        dist_metric = dm;
        if nargin >= 5
            max_splits = ms;
            if nargin == 6
                dist_thresh = thr;
            end
        end
    end
end

chg_list = [];
num_clusts = size(Clust, 1);
max_comp  = max(max(max(Comps)));

%start by determining the pixel width of the smallest element in all the
%clusters (only split if the cluster width is at least twice this).
for i=1:num_clusts
    clust_min = min(Clust(i).pos(:,3) - Clust(i).pos(:,1));
    if clust_min < min_width
        min_width = clust_min;
    end
end

%initially put all elements in the refine list if not passed
if nargin < 3
    refine_list = num_clusts:-1:1;
end
keep_list = 1:num_clusts;

while ~ isempty(refine_list)

    %determine if the first item in the list can be refined by spliting it
    %into up to split+1 pieces and matching each part with another element

    % first ensure that it is large enough to be split
    r = refine_list(1);
    fprintf('                                                              \r');
    fprintf('cluster: %d -- ', r);
    c_width = size(Clust(r).avg, 2);
    if (c_width >= 2 * min_width)
        r_pos = find(keep_list == r,1);
        ind = [keep_list(1:r_pos-1), keep_list(r_pos+1:end)];
        mres = find_match(max_splits, Clust(r).avg, Clust, ind, min_width, ...
                          dist_metric, dist_thresh);
        if ~isempty(mres)
            %appropriate match found, reaverage items as appropriate
            fprintf('match found!\r');
            chg_list = [chg_list, mres.clust];
            if display_matches
                num_matches = length(mres.clust) + 1;
                clf;
                subplot(2, num_matches, 1), imshow(Clust(r).avg), ...
                        xlabel('original');
            end

            %first find all components that list one of r's components as a
            %neighbour
            NB = cell(Clust(r).num);
            for i=1:num_clusts
                if i==r
                    continue;
                end
                for j=1:Clust(r).num
                    idx = find(Clust(i).nb == Clust(r).comp(j));
                    if ~isempty(idx)
                        if size(idx,2) > 1
                            %this can happen if Clust(i).nb contains a single 
                            %row, but has multiple matches within that row 
                            %(ex left and top neighbour).  Rare but it does 
                            %appear
                            idx = idx';
                        end
                        NB{j} = [NB{j}, [repmat(i,1,length(idx)); idx']];
                    end
                end
            end

            while ~isempty(mres.clust)
                if ~isempty(mres.sep_pos)
                    %must update the R positions of each item currently in this
                    %cluster (and possibly crop the T or B positions), as well 
                    %as updating the R neighbour (and maybe T and B)
                    newcl.pg   = Clust(r).pg;
                    newcl.num  = Clust(r).num;
                    newcl.comp = max_comp + (1:Clust(r).num)';
                    max_comp   = max_comp + newcl.num;
                    newcl.pos  = Clust(r).pos;
                    newcl.pos(:,3) = Clust(r).pos(:,1) + mres.sep_pos(1) - 1;
                    newcl.avg  = Clust(r).avg(:,1:mres.sep_pos(1));
                    newcl.nb   = Clust(r).nb;
                    newcl.nb(:,3) = Clust(r).comp;
                    [row,Dummy] = find(newcl.avg ~= bg_val);
                    row = sort(row);
                    if row(end) ~= size(newcl.avg,1)
                        %decrease the bottom positions
                        newcl.pos(:,4) = newcl.pos(:,4) - ...
                              (size(newcl.avg,1) - row(end));
                        newcl.avg = newcl.avg(1:row(end),:);
                    end
                    if row(1) ~= 1
                        %push top positions up
                        newcl.pos(:,2) = newcl.pos(:,2) + row(1) - 1;
                        newcl.avg = newcl.avg(row(1):end,:);
                    end

                    %update Comps to reflect the new component numbers
                    for i=1:newcl.num
                        x = newcl.pos(i,:);
                        pg = newcl.pg(i);
                        region = Comps(x(2):x(4), x(1):x(3), pg);
                        idx = find(region == Clust(r).comp(i));
                        region(idx) = newcl.comp(i);
                        Comps(x(2):x(4), x(1):x(3), pg) = region;
                    end

                    %loop through each component, checking whether the top and 
                    %bottom neighbours overlap the new component by at least 
                    %half the pixels, if they don't we may have to look for a 
                    %new top and/or bottom neighbour that matches better
                    %
                    %we also potentially update any neighbours pointing to this
                    %unsplit component to the appropriate overlapping split part
                    for i=1:newcl.num
                        tl = newcl.pos(i,1); tt = newcl.pos(i,2);
                        tr = newcl.pos(i,3); tb = newcl.pos(i,4);

                        %check the overlap with the current top neighbour
                        if newcl.nb(i,2) ~= non_nb_val
                            [tnbcl, tnboff] = get_cl_off(Clust, newcl.nb(i,2));
                            dif = min(tr, Clust(tnbcl).pos(tnboff,3)) - ...
                                  max(tl, Clust(tnbcl).pos(tnboff,1));
                        else 
                            dif = -Inf;
                        end
                        if dif < ((tr - tl) / 2)
                            vals = [];
                            col = tl;
                            while col ~= tr;
                                row = tt -1;
                                while row >= 1 && ...
                                      Comps(row, col, newcl.pg(i)) == bg_val
                                    row = row - 1;
                                end
                                if row ~= non_nb_val
                                    vals = [vals, Comps(row, col, newcl.pg(i))];
                                end
                                col = col + 1;
                            end
                            %the top neighbour is the maximum occuring value in 
                            %vals
                            if isempty(vals)
                                newcl.nb(i,2) = non_nb_val;  %no top neighbour
                            else
                                [Dummy, idx] = max(hist(vals, max(vals) - ...
                                               min(vals) + 1));
                                newcl.nb(i,2) = idx + min(vals) - 1;
                            end
                        end

                        %check the overlap with current bottom neighbour
                        if newcl.nb(i,4) ~= non_nb_val
                            [bnbcl, bnboff] = get_cl_off(Clust, newcl.nb(i,4));
                            dif = min(tr, Clust(bnbcl).pos(bnboff,3)) - ...
                                  max(tl, Clust(bnbcl).pos(bnboff,1));
                        else
                            dif = -Inf;
                        end
                        if dif < ((tr - tl) / 2)
                            vals = [];
                            col = tl;
                            while col ~= tr;
                                row = tb +1;
                                while row <= size(Comps,1) && ...
                                      Comps(row, col, newcl.pg(i)) == bg_val
                                    row = row + 1;
                                end
                                if row <= size(Comps,1)
                                    vals = [vals, Comps(row, col, newcl.pg(i))];
                                end
                                col = col + 1;
                            end
                            %the bottom neighbour is the maximum occuring value 
                            %in vals
                            if isempty(vals)
                                newcl.nb(i,4) = non_nb_val; %no bottom neighbour
                            else
                                [Dummy, idx] = max(hist(vals, max(vals) - ...
                                               min(vals) + 1));
                                newcl.nb(i,4) = idx + min(vals) - 1;
                            end
                        end

                        %now see if the component numbers of any of those 
                        %components which list this (unsplit) component as a 
                        %neighbour, should be updated to the now split component
                        keep_cols = [];
                        for k=1:size(NB{i},2)
                            [row,col] = ind2sub([Clust(NB{i}(1,k)).num, 4], ...
                                                NB{i}(2,k));
                            %any left neighbours should be associated with the
                            %first component, any right neighbours should be 
                            %associated with the last component
                            if col == 3 || ... %left neighbour
                               (col ~= 1 && ... % not right neighbour AND
                               ((newcl.pos(i,3) - ...
                               Clust(NB{i}(1,k)).pos(row,1)) >= ...
                               ((Clust(NB{i}(1,k)).pos(row,3) - ...
                               Clust(NB{i}(1,k)).pos(row,1))/2))) %max overlap
                                                     %top or bottom neighbour
                                %in either of these cases point this component 
                                %at the new split component
                                Clust(NB{i}(1,k)).nb(row,col) = newcl.comp(i);
                            else
                                keep_cols = [keep_cols,k];
                            end
                        end
                        NB{i} = NB{i}(:,keep_cols);
                    end

                    % update r's L position and neighbour of each item so it 
                    % is set correctly next time through the loop
                    Clust(r).pos(:,1) = Clust(r).pos(:,1) + mres.sep_pos(1);
                    Clust(r).avg = Clust(r).avg(:,mres.sep_pos(1)+1:end);
                    Clust(r).nb(:,1) = newcl.comp;
                    mres.sep_pos = mres.sep_pos(2:end);
                else
                    %last piece, everything should already be updated in terms
                    %of neighbours and component numbers etc.
                    newcl = Clust(r);

                    [row,Dummy] = find(newcl.avg ~= bg_val);
                    row = sort(row);
                    if row(end) ~= size(newcl.avg,1)
                        %decrease the bottom positions
                        newcl.pos(:,4) = newcl.pos(:,4) - ...
                              (size(newcl.avg,1) - row(end));
                        newcl.avg = newcl.avg(1:row(end),:);
                    end
                    if row(1) ~= 1
                        %push top positions up
                        newcl.pos(:,2) = newcl.pos(:,2) + row(1) - 1;
                        newcl.avg = newcl.avg(row(1):end,:);
                    end
                end
                if display_matches
                    %add newcl.avg to the display
                    pos = num_matches - length(mres.clust) + 1;
                    subplot(2, num_matches, pos), imshow(newcl.avg), ...
                        xlabel(sprintf('piece %d\n', pos - 1));
                    subplot(2, num_matches, num_matches + pos), ...
                        imshow(Clust(mres.clust(1)).avg), xlabel(...
                        sprintf('matching cluster %d\n', mres.clust(1)));
                end
                %now re-average the items and update mres
                Clust(mres.clust(1)) = add_and_reaverage(...
                           Clust(mres.clust(1)), newcl, mres.m{1}, mres.m2{1});
                mres.clust = mres.clust(2:end);
                mres.m = mres.m(2:end);
                mres.m2 = mres.m2(2:end);
            end
            % remove r from the cluster list
            keep_list = ind;

            if display_matches
                drawnow;
                pause(2);
            end
        else
            fprintf('no match\r');
        end
    end
    refine_list = refine_list(2:end);
end
fprintf('\n');

% remove those items not on the keep list, and sort the results
[Dummy, Dummy, chg_list] = intersect(unique(chg_list), keep_list);
%Clust = sort_clusters(Clust(keep_list));
Clust = Clust(keep_list);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function match = find_match(num_splits, C, Clust, ind, min_width, ...
                 dist_metric, thresh)
%this function attempts to recursively match halves using at most num_splits
%splits.
%match will be empty if no successful match could be found otherwise, it will
%be a struct with fields: match.clust, match.sep_pos, match.m, and match.m2.  
%The first field will contain at least 1 scalar item specifying the matching 
%cluster index, the second will contain size(match.clust)-1 items listing the
%positions of the separators used to find the matchings, and the last 2 fields
%are cell arrays of resized maximally overlapping versions of Clust(i) (in m) 
%and C (in m2) that were used to find the matching.

bg_val = 0;  %value used for backgroundpixels

%always attempt to find a single match over all of C
for i=1:size(Clust,1)
    if isempty(find(ind == i,1))
        continue;
    end
    if strcmp(dist_metric,'euc')
        [match_l, Mi, Mc] = euc_match(Clust(i).avg, C, thresh);
    elseif strcmp(dist_metric, 'conv_euc')
        [match_l, Mi, Mc] = conv_euc_match(Clust(i).avg, C, thresh);
    elseif strcmp(dist_metric, 'hausdorff')
        [match_l, Mi, Mc] = hausdorff_match(Clust(i).avg, C, thresh);
    elseif strcmp(dist_metric, 'ham')
        [match_l, Mi, Mc] = ham_match(Clust(i).avg, C, thresh);
    else
        error('incorrect distance metric specified!');
    end

    if match_l
        %match found
        match.clust = i;
        match.sep_pos = [];
        match.m = {Mi};
        match.m2 = {Mc};
        return;
    end
end

if num_splits >= 1
    %recursive case, first look for a left-side match
    c_width = size(C,2);
    for i=1:size(Clust,1)
        if isempty(find(ind == i,1))
            continue;
        end
        cut_pos = min_width;
        while (cut_pos <= (c_width - min_width))
            L = C(:,1:cut_pos);
            [row, col] = find(L ~= bg_val);
            row = sort(row);
            col = sort(col);
            L = L(row(1):row(end),col(1):col(end));
            R = C(:,cut_pos+1:end);
            if strcmp(dist_metric,'euc')
                [match_l, Mi, Mc] = euc_match(Clust(i).avg, L, thresh);
            elseif strcmp(dist_metric, 'conv_euc')
                [match_l, Mi, Mc] = conv_euc_match(Clust(i).avg, L, thresh);
            elseif strcmp(dist_metric, 'hausdorff')
                [match_l, Mi, Mc] = hausdorff_match(Clust(i).avg, L, thresh);
            elseif strcmp(dist_metric, 'ham')
                [match_l, Mi, Mc] = ham_match(Clust(i).avg, L, thresh);
            else
                error('incorrect distance metric specified!');
            end
            if match_l
                %ensure the right half matches in up to one fewer splits
                [row, col] = find(R ~= bg_val);
                row = sort(row);
                R = R(row(1):row(end),:);
                mres = find_match(num_splits-1, R, Clust, ind, min_width, ...
                                  dist_metric, thresh);
                if ~isempty(mres)
                    %right match found!
                    match.clust = [i, mres.clust];
                    match.sep_pos = [cut_pos, mres.sep_pos];
                    match.m = {Mi, mres.m{:}};
                    match.m2 = {Mc, mres.m2{:}};
                    return;
                end
            end
            cut_pos = cut_pos + 1;
        end
    end
end
%if we get down here, a match hasn't been found, set match to an empty array
match = [];

