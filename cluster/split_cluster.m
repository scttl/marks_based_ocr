function [Clust,Comps] = split_cluster(Clust,Comps,idx,pos,varargin)
%  SPLIT_CLUSTER  Split a cluster into two clusters
%
%  [CLUST,COMPS] = split_cluster(CLUST, COMPS, IDX, POS, [VA1, VL1]...)
%
%  CLUST should be a struct containing cluster information.  See cluster_comps
%
%  COMPS should be a struct containing component information.  See get_comps
%
%  IDX should be a scalar indexing the Clust array specifying which cluster
%  to split into two.
%
%  POS should specify the column (in the cluster average image) at which to
%  split the cluster into two (it should specify the first column of the right
%  half of the split
%
%  One can override the default values for any of the LOCAL VARS defined below
%  by passing in name and value pairs.  VA1 should be a string representing the
%  name of the variable to override, and VL1 should be its new value.
%
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: split_cluster.m,v 1.1 2007-02-01 18:09:45 scottl Exp $
%
% REVISION HISTORY
% $Log: split_cluster.m,v $
% Revision 1.1  2007-02-01 18:09:45  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%by default we make a new cluster for the right half of the split.  Setting
%this to false, means we make a new cluster for the left half.
keep_right = false;


% CODE START %
%%%%%%%%%%%%%%

if nargin < 4
    error('incorrect number of arguments passed');
elseif nargin > 4
    process_optional_args(varargin{:});
end

if idx < 1 || idx > Clust.num
    error('first index passed is invalid: %d, clust size: %d', idx, Clust.num);
end

num_cols = size(Clust.avg{idx},2);
num_rows = size(Clust.avg{idx},1);

if pos < 1 || pos > num_cols
    error('invalid pos: %d, cluster size %d', pos, num_cols);
end

%update the cluster fields
new_idx = Clust.num + 1;
if keep_right
    [Clust.avg{new_idx},l,t,r,b] = crop_image(Clust.avg{idx}(:,1:pos-1));
    [Clust.avg{idx},rl,rt,rr,rb] = crop_image(Clust.avg{idx}(:,pos:num_cols));
else
    [Clust.avg{new_idx},rl,rt,rr,rb]=crop_image(Clust.avg{idx}(:,pos:num_cols));
    [Clust.avg{idx},l,t,r,b] = crop_image(Clust.avg{idx}(:,1:pos-1));
end
left_adj = [l-1, t-1, -(num_cols-r), -(num_rows-b)];
right_adj = [rl+pos-1, rt-1, -(num_cols - (rr+pos)), -(num_rows - rb)];

Clust.num = Clust.num+1;
Clust.num_comps(new_idx) = Clust.num_comps(idx);
Clust.mode_num(new_idx) = Clust.mode_num(idx);
Clust.norm_sq(idx) = sum(Clust.avg{idx}(:).^2);
Clust.norm_sq(new_idx) = sum(Clust.avg{new_idx}(:).^2);
Clust.refined(new_idx) = false;
Clust.changed(new_idx) = Clust.changed(idx);

if ~isempty(Clust.bigram)
    %@@@ this should be fixed (counts re-normalized)
    Clust.bigram(new_idx,:) = Clust.bigram(idx,:);
    Clust.bigram(:,new_idx) = Clust.bigram(:,idx);
end
if ~isempty(Clust.pos_count)
    for ii=1:length(Clust.pos_count)
        Clust.pos_count{ii}(new_idx,:) = Clust.pos_count{ii}(idx,:);
    end
end
if ~isempty(Clust.descender_off)
    Clust.descender_off(new_idx) = Clust.descender_off(idx);
end
if ~isempty(Clust.ascender_off)
    Clust.ascender_off(new_idx) = Clust.ascender_off(idx);
end
if ~isempty(Clust.class)
    Clust.class(new_idx) = Clust.class(idx);
end
if Clust.found_true_labels
    Clust.truth_label{new_idx} = Clust.truth_label{idx};
end

%update components by creating new ones
num_new_comps = Clust.num_comps(new_idx);
comp_idcs = Clust.comps{idx};
new_comp_idcs = Comps.max_comp+1:Comps.max_comp+num_new_comps;
Clust.comps{new_idx} = new_comp_idcs';
Comps.max_comp = Comps.max_comp + num_new_comps;
Comps.clust = [Comps.clust; new_idx + zeros(num_new_comps,1,'uint16')];
Comps.pos = [Comps.pos; Comps.pos(comp_idcs,:)];
Comps.pg = [Comps.pg; Comps.pg(comp_idcs) + zeros(num_new_comps,1, 'uint32')];
Comps.nb = [Comps.nb; Comps.nb(comp_idcs,:)];
Comps.nb_dist = [Comps.nb_dist; Comps.nb_dist(comp_idcs,:)];

%adjust component positions, neighbours etc
if keep_right
    Comps.pos(comp_idcs,:) = uint16(double(Comps.pos(comp_idcs,:)) + ...
                             repmat(right_adj, length(comp_idcs),1));
    Comps.pos(new_comp_idcs,:) = uint16(double(Comps.pos(new_comp_idcs,:)) + ...
                                 repmat(left_adj, num_new_comps,1));
    Comps.nb(comp_idcs,1) = new_comp_idcs;
    Comps.nb(new_comp_idcs,3) = comp_idcs;
    Comps.nb_dist(comp_idcs,1) = 1;
    Comps.nb_dist(new_comp_idcs,3) = 1;
    right_adj = abs(right_adj);
    left_adj = abs(left_adj);
    Comps.nb_dist(comp_idcs,:) = Comps.nb_dist(comp_idcs,:) + ...
                                 uint16(repmat(right_adj,length(comp_idcs),1));
    Comps.nb_dist(new_comp_idcs,:) = Comps.nb_dist(new_comp_idcs,:) + ...
                                 uint16(repmat(left_adj, num_new_comps,1));
    %@@neighbours that point to these guys must be updated too
else
    Comps.pos(comp_idcs,:) = uint16(double(Comps.pos(comp_idcs,:)) + ...
                             repmat(left_adj, length(comp_idcs),1));
    Comps.pos(new_comp_idcs,:) = uint16(double(Comps.pos(new_comp_idcs,:)) + ...
                                 repmat(right_adj, num_new_comps,1));
    Comps.nb(comp_idcs,3) = new_comp_idcs;
    Comps.nb(new_comp_idcs,1) = comp_idcs;
    Comps.nb_dist(comp_idcs,3) = 1;
    Comps.nb_dist(new_comp_idcs,1) = 1;
    right_adj = abs(right_adj);
    left_adj = abs(left_adj);
    Comps.nb_dist(comp_idcs,:) = Comps.nb_dist(comp_idcs,:) + ...
                                 uint16(repmat(left_adj,length(comp_idcs),1));
    Comps.nb_dist(new_comp_idcs,:) = Comps.nb_dist(new_comp_idcs,:) + ...
                                 uint16(repmat(right_adj, num_new_comps,1));
    %@@neighbours that point to these guys must be updated too
end

if Comps.found_lines
    Comps.line = [Comps.line; zeros(num_new_comps,1,'uint64')];
    Comps.line(new_comp_idcs) = Comps.line(comp_idcs);
    Comps.ascender_off = [Comps.ascender_off; Comps.ascender_off(comp_idcs)];
    Comps.descender_off = [Comps.descender_off; Comps.descender_off(comp_idcs)];
    Comps.scale_factor = [Comps.scale_factor; Comps.scale_factor(comp_idcs)];
end
