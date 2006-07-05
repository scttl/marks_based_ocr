function [Clust, Comps] = add_and_reaverage(Clust, Comps, id1, id2)
%  ADD_AND_REAVERAGE  Merge the contents of one cluster with another.
%
%  [Clust, Comps] = add_and_reaverage(Clust, Comps, id1, id2)
%
%  Clust should be a struct containing cluster information.  See cluster_comps
%
%  Comps should be a struct containing component information.  See cluster_comps
%
%  id1 should be a scalar indexing the Clust array to which each index in id2
%  will be added (with components updated to have id1 as their associated 
%  cluster).  Each index in id2 will also be removed and the cluster
%  corresponding to id1 will have its refined status set to true.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: add_and_reaverage.m,v 1.2 2006-07-05 01:10:24 scottl Exp $
%
% REVISION HISTORY
% $Log: add_and_reaverage.m,v $
% Revision 1.2  2006-07-05 01:10:24  scottl
% rewritten to process multiple items simultaneously, and handle the new
% Cluster and component structures.
%
% Revision 1.1  2006/06/03 20:55:47  scottl
% Initial check-in.
%

% LOCAL VARS %
%%%%%%%%%%%%%%
resize_method = 'nearest';
clust_size = Clust.num;

% CODE START %
%%%%%%%%%%%%%%

if id1 < 1 || id1 > clust_size
    error('first index passed is invalid: %d, clust size: %d', id1, clust_size);
elseif any(id2 < 1) || any(id2 > clust_size)
    error('one of the 2nd indices passed is invalid');
end

%ensure we don't try and add a cluster to itself
id2 = id2(id2 ~= id1);
if length(id2) == 0
    Clust.refined(id1) = true;
    return;
end

%update the components belonging to any of id2's clusters to id1
chg_comps = cell2mat(Clust.comps(id2));
Comps.clust(chg_comps) = id1;
Clust.comps{id1} = [Clust.comps{id1}; chg_comps];
num_c1 = Clust.num_comps(id1);
Clust.num_comps(id1) = num_c1 + sum(Clust.num_comps(id2));

for cc = id2
    %update the averages.  Note that the averages may be different sizes, so we 
    %take that into account by rescaling the second average
    num_c2 = Clust.num_comps(cc); num_tot = num_c1 + num_c2; 
    avg1 = Clust.avg{id1}; avg2 = Clust.avg{cc};
    if all(size(avg1) == size(avg2))
        Clust.avg{id1} = (num_c1/num_tot .* avg1) + (num_c2/num_tot .* avg2);
    else
        Clust.avg{id1} = (num_c1/num_tot .* avg1) + (num_c2/num_tot .* ...
                         imresize(avg2, size(avg1), resize_method));
    end
end

%update the squared norm based on this new average.
Clust.norm_sq(id1) = sum(sum(Clust.avg{id1} .^2));

%set the status of this cluster to refined.
Clust.refined(id1) = true;

%remove the second list of clusters
keep_list = setdiff(1:Clust.num, id2);
Clust.num = length(keep_list);
Clust.num_comps = Clust.num_comps(keep_list);
Clust.comps = Clust.comps(keep_list);
Clust.avg = Clust.avg(keep_list);
Clust.norm_sq = Clust.norm_sq(keep_list);
Clust.refined = Clust.refined(keep_list);
Clust.offset = Clust.offset(keep_list);
%we also must update the cluster id associated with the components 
%since these may have been shifted around by removing clusters
for ii = 1:Clust.num
    Comps.clust(Clust.comps{ii}) = ii;
end
