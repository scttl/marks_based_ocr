function [Clust,Comps,newid1] = add_and_reaverage(Clust,Comps,id1,id2,varargin)
%  ADD_AND_REAVERAGE  Merge the contents of one cluster with another.
%
%  [Clust,Comps,newid1] = add_and_reaverage(Clust,Comps,id1,id2, [VA1, VL1]...)
%
%  Clust should be a struct containing cluster information.  See cluster_comps
%
%  Comps should be a struct containing component information.  See get_comps
%
%  id1 should be a scalar indexing the Clust array to which each index in id2
%  will be added (with components updated to have id1 as their associated 
%  cluster).  Each index in id2 will also be removed.
%
%  One can override the default values for any of the LOCAL VARS defined below
%  by passing in name and value pairs.  VA1 should be a string representing the
%  name of the variable to override, and VL1 should be its new value.
%
%  newid1 will be the scalar index corresponding to the Cluster represented by
%  id1 at the start of the call.
%
%  NOTE: since Cluster entries are removed, there are no guarantees that the
%        index into Clusters is preserved upon the completion of this call.
%        Thus newid1 can be used to get the new index to the original Cluster
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: add_and_reaverage.m,v 1.10 2007-02-01 18:10:12 scottl Exp $
%
% REVISION HISTORY
% $Log: add_and_reaverage.m,v $
% Revision 1.10  2007-02-01 18:10:12  scottl
% added new class field that is assigned based on offset information
%
% Revision 1.9  2006-12-19 22:13:42  scottl
% added pos_count field.  Implemented ability to return clusters
% without refining.
%
% Revision 1.8  2006-10-29 17:24:54  scottl
% change to cluster struct, to use descender and ascender offsets, instead
% of a single offset field.
%
% Revision 1.7  2006-10-09 16:35:40  scottl
% changes to argument processing and variable overriding.
%
% Revision 1.6  2006/08/24 21:40:04  scottl
% added ability to use the mode instead of taking the average of cluster
% intensities while refining.
%
% Revision 1.5  2006/08/14 01:33:27  scottl
% remove ability to mark the new cluster as refined.  Since we return the index
% of this cluster, this can be done in the calling function.
%
% Revision 1.4  2006/08/07 21:20:07  scottl
% add ability to return the updated index of the first cluster.
%
% Revision 1.3  2006/07/06 17:50:32  scottl
% added ability to mark the merged cluster as refined.
%
% Revision 1.2  2006/07/05 01:10:24  scottl
% rewritten to process multiple items simultaneously, and handle the new
% Cluster and component structures.
%
% Revision 1.1  2006/06/03 20:55:47  scottl
% Initial check-in.
%

% LOCAL VARS %
%%%%%%%%%%%%%%

%when merging different sized clusters, how should we resize?  see imresize()
resize_method = 'nearest';

%use_avg, when set to true will take the average over all the merged cluster
%average intensity images to create the new average intensity image.  If false,
%then the resulting average intensity image will be set to the cluster that 
%has the most components (no averaging amongst all the clusters being merged 
%will be done).  i.e. the mode will be taken.
use_avg = true;

clust_size = Clust.num;


% CODE START %
%%%%%%%%%%%%%%

if nargin < 4
    error('incorrect number of arguments passed');
elseif nargin > 4
    process_optional_args(varargin{:});
end

if id1 < 1 || id1 > clust_size
    error('first index passed is invalid: %d, clust size: %d', id1, clust_size);
elseif any(id2 < 1) || any(id2 > clust_size)
    error('one of the 2nd indices passed is invalid');
end

%ensure we don't try and add a cluster to itself
id2 = id2(id2 ~= id1);
if length(id2) == 0
    return;
end

%update the components belonging to any of id2's clusters to id1
chg_comps = cell2mat(Clust.comps(id2));
Comps.clust(chg_comps) = id1;
Clust.comps{id1} = [Clust.comps{id1}; chg_comps];
num_c1 = Clust.num_comps(id1);
Clust.num_comps(id1) = num_c1 + sum(Clust.num_comps(id2));

if use_avg
    for cc = id2'
        %update the averages.  Note that the averages may be different sizes, so
        %we take that into account by rescaling the second average
        num_c2 = Clust.num_comps(cc); num_tot = num_c1 + num_c2; 
        avg1 = Clust.avg{id1}; avg2 = Clust.avg{cc};
        if all(size(avg1) == size(avg2))
            Clust.avg{id1} = (num_c1/num_tot.*avg1) + (num_c2/num_tot.*avg2);
        else
            Clust.avg{id1} = (num_c1/num_tot .* avg1) + (num_c2/num_tot .* ...
                             imresize(avg2, size(avg1), resize_method));
        end
    end
    Clust.mode_num(id1) = Clust.mode_num(id1) + sum(Clust.mode_num(id2));
else
    %take the mode
    idcs = [id1; id2];
    [md_val,md_idx] = max(Clust.mode_num(idcs));
    Clust.avg{id1} = Clust.avg{idcs(md_idx)};
    Clust.mode_num(id1) = md_val;
end
%update the squared norm based on this new average.
Clust.norm_sq(id1) = sum(Clust.avg{id1}(:) .^2);

%remove the second list of clusters
keep_list = setdiff(1:Clust.num, id2);
Clust.num = length(keep_list);
Clust.num_comps = Clust.num_comps(keep_list);
Clust.mode_num = Clust.mode_num(keep_list);
Clust.comps = Clust.comps(keep_list);
Clust.avg = Clust.avg(keep_list);
Clust.norm_sq = Clust.norm_sq(keep_list);
Clust.refined = Clust.refined(keep_list);
Clust.changed = Clust.changed(keep_list);
if Clust.found_offsets
    Clust.descender_off = Clust.descender_off(keep_list);
    Clust.ascender_off = Clust.ascender_off(keep_list);
    Clust.class = Clust.class(keep_list);
end
if ~isempty(Clust.bigram)
    Clust.bigram = Clust.bigram(keep_list, keep_list);
end
if ~isempty(Clust.pos_count)
    for ii=1:length(Clust.pos_count)
        Clust.pos_count{ii} = Clust.pos_count{ii}(keep_list,:);
    end
end

%we also must update the cluster id associated with the components 
%since these may have been shifted around by removing clusters
for ii = 1:Clust.num
    Comps.clust(Clust.comps{ii}) = ii;
end

%determine the new index of the original Cluster to which the other clusters
%were merged.
newid1 = id1 - sum(id2 < id1);
