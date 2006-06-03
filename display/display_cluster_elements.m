function display_cluster_elements(Comps, Clust, x)
% DISPLAY_CLUSTER_ELEMENTS   Display the first x elements of each Cluster
%
%   DISPLAY_CLUSTER_ELEMENTS(Comps, Clust, x)
%
%   Clust should be an array of structs, each of which is assumed to contain 
%   a pos field which is a four-tuple of co-ordinates (L,T,R,B) specifying a
%   submatrix of Comps to be drawn.
%
%   Comps should be a binary image.
%
%   x should be the number of elements to draw for each cluster (defaults to
%   10 if there are that many elements).  It can also be the string 'all', in
%   which case all elements of only the first element passed are drawn, with 
%   10 elements per row.

% CVS INFO %
%%%%%%%%%%%%
% $Id: display_cluster_elements.m,v 1.1 2006-06-03 20:55:54 scottl Exp $
%
% REVISION HISTORY
% $Log: display_cluster_elements.m,v $
% Revision 1.1  2006-06-03 20:55:54  scottl
% Initial check-in.
%
%

% LOCAL VARS %
%%%%%%%%%%%%%%
num_to_draw = 10;
num_clusts = size(Clust,1);
draw_one_clust = false;


% CODE START %
%%%%%%%%%%%%%%

if nargin >= 3
    if ischar(x)
        draw_one_clust = true;
    else
        num_to_draw = x;
    end
end

%binarize the input matrix so we can draw in black and white (if not already)
Comps = Comps ~= 0;
if draw_one_clust
    num_rows = ceil(Clust(1).num / num_to_draw);
    for j = 1:Clust(1).num
        l = Clust(1).pos(j, 1);
        t = Clust(1).pos(j, 2);
        r = Clust(1).pos(j, 3);
        b = Clust(1).pos(j, 4);
        subplot(num_rows, num_to_draw, j), ...
        imagesc(Comps(t:b,l:r)), axis equal, axis off, colormap gray;
    end
else
    for i=1:num_clusts;
        j = 1;
        while j <= Clust(i).num & j <= num_to_draw
            l = Clust(i).pos(j, 1);
            t = Clust(i).pos(j, 2);
            r = Clust(i).pos(j, 3);
            b = Clust(i).pos(j, 4);
            subplot(num_clusts, num_to_draw, (num_to_draw *(i-1))+j), ...
            imagesc(Comps(t:b,l:r)), axis equal, axis off, colormap gray;
            j = j+1;
        end
    end
end
