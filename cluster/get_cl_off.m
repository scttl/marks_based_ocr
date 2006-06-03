function [clust, off] = get_cl_off(Clust, c)
% GET_CL_OFF    Locate the cluster and offset index of the component passed
%
%  [cluster, offset] = get_cl_off(Clust, component)
%
%  Given a Clust structure and a single component, this routine returns the 
%  cluster index and offset locating this component if it exists, otherwise both
%  values are set to NaN, and a warning is displayed.
%
%  Clust should be an array of structs, each of which is assumed to contain
%  several fields of a particular format.  See cluster_comps for details.

% CVS INFO %
%%%%%%%%%%%%
% $Id: get_cl_off.m,v 1.1 2006-06-03 20:55:48 scottl Exp $
%
% REVISION HISTORY
% $Log: get_cl_off.m,v $
% Revision 1.1  2006-06-03 20:55:48  scottl
% Initial check-in.
%
%

% LOCAL VARS %
%%%%%%%%%%%%%%

% CODE START %
%%%%%%%%%%%%%%
if nargin ~= 2
    error('incorrect number of arguments specified!');
end

for i=1:size(Clust,1)
    idx = find(Clust(i).comp == c,1);
    if ~isempty(idx)
        clust = i;
        off = idx;
        return;
    end
end
clust = NaN;
off = NaN;
warning(sprintf('component %d not found in cluster structure', c));
