function C = sort_clusters(C)
% SORT_CLUSTERS   Sort (by descending number of elements) the Clusters passed
%
%   C = SORT_CLUSTERS(C)
%
%   C should be an array of structs, each of which is assumed to contain an
%   integer valued num field.  This field represent the number of elements
%   belonging to that cluster, and will be our key for sorting.
%
%   The array elements are returned unchanged but in decreasing order
%   according to the num field.
%
%   NOTE: due to the nature of the implementation, all clusters with the same
%   number of elements end up being flipped so that the former first listed
%   element ends up last and vice versa (within the group of clusters with the
%   same number of elements).
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: sort_clusters.m,v 1.1 2006-06-03 20:55:48 scottl Exp $
%
% REVISION HISTORY
% $Log: sort_clusters.m,v $
% Revision 1.1  2006-06-03 20:55:48  scottl
% Initial check-in.
%
%

% LOCAL VARS %
%%%%%%%%%%%%%%


% CODE START %
%%%%%%%%%%%%%%
if nargin ~= 1
    error('incorrect number of arguments passed!');
end

% note we need to flip it to get decreasing order below, but this throws
% elements with the same number out of order.
[Dummy, I] = sort([C.num]);
C = C(fliplr(I));
