function display_cluster_elements(Clust, Comps, idx, num, comp_idx)
% DISPLAY_CLUSTER_ELEMENTS   Display the elements of clusters
%
%   DISPLAY_CLUSTER_ELEMENTS(Clust, Comps, idx, [num, comp_idx])
%
%   Clust and Comps should be as defined in cluster_comps.m
%
%   Comps should be a binary image.
%
%   idx should be the list of Cluster indices whose elements should be drawn.
%
%   num is optional and if specified determines the (maximum) number of elements
%   that should be drawn for each cluster.  If not specified it defaults to 50.
%   The elements are randomly sampled from each cluster (containing > 50 
%   elements, unless comp_idx is passed)
%
%   comp_idx is optional and if specified gives an explicit ordering of the
%   components to display.  If not passed, some random subsample of components
%   is used.  It should be passed as a column vector with elements belonging to
%   the same Cluster, listed in the order they would like to be drawn.  Note
%   that when passed, even though num is required, it will be ignored, and the
%   number of components to draw is calculated based on the list passed.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: display_cluster_elements.m,v 1.4 2006-08-24 21:37:37 scottl Exp $
%
% REVISION HISTORY
% $Log: display_cluster_elements.m,v $
% Revision 1.4  2006-08-24 21:37:37  scottl
% added ability to specify particular component indices to be displayed.
%
% Revision 1.3  2006/08/07 21:19:23  scottl
% remove dependence on imview
%
% Revision 1.2  2006/07/21 20:27:09  scottl
% rewritten based on new Clust and Comps structures.
%
% Revision 1.1  2006/06/03 20:55:54  scottl
% Initial check-in.


% LOCAL VARS %
%%%%%%%%%%%%%%
num_to_draw = 50;
row_pix_border = 10;
col_pix_border = 3;
use_random = true;  %by default use a random sample of elements from the cluster

%set save_elements to true to write the elements image to disk based on the 
%params below it
save_elements = false;
img_prefix = 'results/cluster_elements';
img_format = 'png';


% CODE START %
%%%%%%%%%%%%%%

if nargin < 3 || nargin > 5
    error('incorrect number of arguments specified!');
elseif nargin >= 4
    num_to_draw = num;
    if nargin == 5
        use_random = false;
    end
end

num_clusts = length(idx);
if use_random
    %collect the random Component indices of each cluster
    comp_list = [];
    lrg_idx = find(Clust.num_comps(idx) > num_to_draw);
    lrg_idx = idx(lrg_idx);
    sm_idx = setdiff(idx, lrg_idx);
    for ii=lrg_idx
        this_comps = Clust.comps{ii}(ceil(rand(num_to_draw,1) .* ...
                         Clust.num_comps(ii)));
        comp_list = [comp_list; this_comps];
    end
    comp_list = [comp_list; cell2mat(Clust.comps(sm_idx))];
else
    %use the component indices passed
    comp_list = comp_idx;
    [num_to_draw,num_to_draw] = mode(Comps.clust(comp_idx));
end

%create a cell arrray to hold the component images
M = cell(num_clusts,num_to_draw);
num_added = zeros(num_clusts,1);
size_vals = zeros(num_clusts,num_to_draw,2);
pgs_to_check = unique(Comps.pg(comp_list));

%create a mapping for cluster indices into rows of M
map = zeros(Clust.num,1);
map(idx) = 1:num_clusts;

%extract the required components from each page and add them to M in the
%appropriate position
for pp = pgs_to_check'
    Img = imread(Comps.files{pp});
    this_comps = find(Comps.pg(comp_list) == pp);
    this_comps = comp_list(this_comps);
    for ii = this_comps'
        row = map(Comps.clust(ii));
        num_added(row) = num_added(row) + 1;
        pos = Comps.pos(ii,:);
        M{row, num_added(row)} = ~Img(pos(2):pos(4), pos(1):pos(3));
        size_vals(row,num_added(row),:) = ...
                 reshape(size(M{row, num_added(row)}),1,1,2);
    end
end

%convert M to an appropriately spaced image
col_w = max(max(size_vals(:,:,2)));
MM = zeros(sum(max(size_vals(:,:,1),[],2))+ num_clusts*row_pix_border, ...
           (col_w + col_pix_border) * num_to_draw);
row = 1;
for ii=1:size(M,1)
    col = 1;
    for jj=1:num_added(ii)
        MM(row:row+size_vals(ii,jj,1)-1, col:col+size_vals(ii,jj,2)-1) = ...
           M{ii,jj};
        col = col + col_w + col_pix_border - 1;
    end
    row = row + max(size_vals(ii,:,1)) + row_pix_border - 1;
end
imshow(MM);

%save the image to disk if required.
if save_elements
    fprintf('writing cluster elements image to disk\n');
    imwrite(MM, [img_prefix, '.', img_format], img_format);
end
