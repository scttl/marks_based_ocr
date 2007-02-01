function [Clust, Comps] = subimage_refine(Clust, Comps, varargin)
% SUBIMAGE_REFINE   Refine cluster via repeated Hausdorff subimage matching
%
%   [CLUST, COMPS] = SUBIMAGE_REFINE(CLUST, COMPS, [VAR1, VAL1]...)
%
%   CLUST should be a struct like that returned from cluster_comps
%
%   COMPS should be a struct like that returned from get_comps
%
%   VAR1 and VAL1 are optional, and can be used to override the default values
%   for the LOCAL VARS defined below.  Each VAR1 should be a string giving the
%   exact name of the variable to override.  Each VAL1 should be the updated
%   value to use for that variable.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: subimage_refine.m,v 1.1 2007-02-01 18:09:45 scottl Exp $
%
% REVISION HISTORY
% $Log: subimage_refine.m,v $
% Revision 1.1  2007-02-01 18:09:45  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%how close must the Hausdorff distance be to be considered a match?
haus_match_thresh = 0; %1;

%at what point are cluster pixels considered on?
on_thresh = 0.5;

%how close to each edge should a match be to prevent being split
edge_tol = 3;

%how wide does a cluster have to be to be considered for refinement?
min_width = 5;

%how tall does a cluster have to be to be considered for refinement?
min_height = 5;

%how much padding should we put between our glued together cluster averages
pad = 1;

%should we display matches on screen as we find them?
display_matches = true;


% CODE START %
%%%%%%%%%%%%%%

%process input arguments
if nargin < 2
    error('Clust and Comps struct must be specified');
elseif nargin > 2
    process_optional_args(varargin{:});
end

img = cell(1,Clust.num);
dist = cell(1,Clust.num);
sz = zeros(Clust.num,2);
top_pad = ceil(1+haus_match_thresh);
for ii=1:Clust.num
    img{ii} = Clust.avg{ii} >= on_thresh;  
    %pad the top and bottom of the image to ensure good vertical matches
    pad_amt = zeros(top_pad,size(img{ii},2));
    img{ii} = [pad_amt; img{ii}; pad_amt];
    sz(ii,:) = size(img{ii});
    dist{ii} = bwdist(img{ii});
end

rr = find(Clust.refined == false, 1, 'first');
while ~isempty(rr)
    fprintf('                                                         \r');
    fprintf('cluster: %d  -- ', rr);

    %ensure its wide enough to consider matching
    if size(Clust.avg{rr},2) < min_width || size(Clust.avg{rr},1) < min_height
        Clust.refined(rr) = true;
        rr = find(Clust.refined == false, 1, 'first');
        continue;
    end

    rem_idx = [1:rr-1, rr+1:Clust.num];

    sub_img = img{rr};
    sub_dist = dist{rr};
    rem_img = imgcell2mat(img(rem_idx), 'elem_padding', pad);
    rem_dist = imgcell2mat(dist(rem_idx), 'elem_padding', pad, 'pad_val', inf);
    %[t,l,s] = haus_subimage_match(sub_img, sub_dist, rem_img, rem_dist, ...
    %                              'thresh', haus_match_thresh);
    [t,l,s] = subimage_match(sub_img, rem_img, 'thresh', haus_match_thresh);

    accum = cumsum(sz(rem_idx,2)+pad);
    accum = [1; accum];
    accum = accum(1:end-1);
    pos = 1;
    %determine which clusters match (and at what point)
    matches = [];
    while ~isempty(l)
        while pos < length(accum) && accum(pos+1) <= l(1)
            pos = pos + 1;
        end
        if pos == length(accum)
            %done
            break;
        else
            %found a match
            matches = [matches; rem_idx(pos), l(1) - accum(pos)];
        end
        l = l(2:end);
    end

    %merge the matching clusters
    matches = unique(matches, 'rows');
    num_matches = size(matches,1);
    fprintf('found %d matches\n', num_matches);
    match_idx = [];
    idx = 1;
    if num_matches > 0 && display_matches
        subplot(1,2,1), imshow(Clust.avg{rr});
        subplot(1,2,2), imshow(imgcell2mat(Clust.avg(unique(matches(:,1)))));
        pause(1);
    end
    while idx <= num_matches
        cl = matches(idx,1);
        pos = matches(idx,2);
        end_pos = pos;
        while idx < num_matches && matches(idx+1,1) == cl && ...
              matches(idx+1,2) == end_pos + 1;
            %multple adjacent matches for this cluster
            idx = idx + 1;
            end_pos = end_pos + 1;
        end
        pos = floor((end_pos+pos)/2);
        if (pos <= edge_tol && sz(cl,2) - (pos + sz(rr,2)) <= edge_tol) || ...
            pos >= sz(cl,2)
            %the match covers the cluster entirely
            fprintf('cluster %d, wholly matches\n', cl);
        else
            %need to split the matching cluster
            if pos <= edge_tol
                %matches left-half
                fprintf('cluster %d, left matches\n', cl);
                [Clust, Comps] = split_cluster(Clust, Comps, cl, pos+sz(rr,2));
                img{Clust.num} = img{cl}(:,pos+sz(rr,2):end);
                dist{Clust.num} = bwdist(img{Clust.num});
                sz(Clust.num,:) = size(img{Clust.num});
                img{cl} = img{cl}(:,1:pos-1);
                dist{cl} = bwdist(img{cl});
                sz(cl,:) = size(img{cl});
                n_idx = idx+1;
                while n_idx <= num_matches && matches(n_idx,1) == cl
                    %multiple (separate) matches in this cluster
                    matches(n_idx,1) = Clust.num;
                    matches(n_idx,2) = matches(n_idx,2) - sz(rr,2);
                    n_idx = n_idx + 1;
                end
            elseif sz(cl,2) - (pos + sz(rr,2)) <= edge_tol
                %matches right-half
                fprintf('cluster %d, right matches\n', cl);
                [Clust, Comps] = split_cluster(Clust, Comps, cl, pos, ...
                                 'keep_right', true);
                img{Clust.num} = img{cl}(:,1:pos-1);
                dist{Clust.num} = bwdist(img{Clust.num});
                sz(Clust.num,:) = size(img{Clust.num});
                img{cl} = img{cl}(:,pos+sz(rr,2):end);
                dist{cl} = bwdist(img{cl});
                sz(cl,:) = size(img{cl});
            else
                %matches in the middle
                fprintf('cluster %d, middle matches\n', cl);
                [Clust, Comps] = split_cluster(Clust, Comps, cl, pos, ...
                                 'keep_right', true);
                img{Clust.num} = img{cl}(:,1:pos-1);
                dist{Clust.num} = bwdist(img{Clust.num});
                sz(Clust.num,:) = size(img{Clust.num});
                img{cl} = img{cl}(:,pos:end);
                dist{cl} = bwdist(img{cl});
                sz(cl,:) = size(img{cl});
                [Clust, Comps] = split_cluster(Clust, Comps, cl, sz(rr,2));
                img{Clust.num} = img{cl}(:,sz(rr,2):end);
                dist{Clust.num} = bwdist(img{Clust.num});
                sz(Clust.num,:) = size(img{Clust.num});
                img{cl} = img{cl}(:,1:sz(rr,2)-1);
                dist{cl} = bwdist(img{cl});
                sz(cl,:) = size(img{cl});
                n_idx = idx+1;
                while n_idx <= num_matches && matches(n_idx,1) == cl
                    %multiple (separate) matches in this cluster
                    matches(n_idx,1) = Clust.num;
                    matches(n_idx,2) = matches(n_idx,2) - ...
                                       (pos + sz(rr,2));
                    n_idx = n_idx + 1;
                end
            end
        end
        match_idx = [match_idx; cl];
        idx = idx + 1;
    end

    keep_idx = setdiff(1:Clust.num, match_idx);
    [Clust, Comps] = add_and_reaverage(Clust, Comps, rr, match_idx);
    Clust.refined(rr) = true;
    sz = sz(keep_idx,:);
    img = img(keep_idx);
    dist = dist(keep_idx);


    rr = find(Clust.refined == false, 1, 'first');
end
fprintf('\n');



% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
