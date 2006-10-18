function display_cluster_elements(Clust, Comps, idx, varargin)
% DISPLAY_CLUSTER_ELEMENTS   Display the elements of clusters
%
%   DISPLAY_CLUSTER_ELEMENTS(Clust, Comps, idx, [VAR1, VAL1]...)
%
%   Clust and Comps should be as defined in cluster_comps.m
%
%   Comps should be a binary image.
%
%   idx should be the list of Cluster indices whose elements should be drawn.
%
%   Additional variables specified in LOCAL VARS below, can be overridden by
%   passing the name of the variable to override as a string, followed by its
%   new value (and repeating for each such variable to override).  This can be
%   used to change the default number of elements to draw for each component, 
%   and to specify exactly what components are drawn (instead of a random 
%   subset).
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: display_cluster_elements.m,v 1.7 2006-10-18 15:59:41 scottl Exp $
%
% REVISION HISTORY
% $Log: display_cluster_elements.m,v $
% Revision 1.7  2006-10-18 15:59:41  scottl
% small change to make use of existing component grabbing functionality
%
% Revision 1.6  2006-10-09 16:33:55  scottl
% changed parameter processing, allowed one to display thinned elements,
% and draw the averages as well.
%
% Revision 1.5  2006/09/05 15:50:44  scottl
% made use of MOCR_PATH variable for saving in the results directory.
%
% Revision 1.4  2006/08/24 21:37:37  scottl
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

%by default use a random sample of elements from the cluster.  This can be
%overridden by specifying an explicit listing of components belonging to each of
%the clustters to be drawn.
comp_list = [];

%should we thin the elements before drawing them?
thin_elements = false;

%should we display the cluster average next to its elements?
display_avg = false;

%set save_elements to true to write the elements image to disk based on the 
%params below it
save_elements = false;
global MOCR_PATH;  %make use of the globally defined MOCR_PATH variable
img_prefix = [MOCR_PATH, '/results/cluster_elements'];
img_format = 'png';


% CODE START %
%%%%%%%%%%%%%%

if nargin < 3
    error('incorrect number of arguments specified!');
elseif nargin > 3
    process_optional_args(varargin{:});
end

num_clusts = length(idx);
if isempty(comp_list)
    %collect the random Component indices of each cluster
    lrg_idx = find(Clust.num_comps(idx) > num_to_draw);
    lrg_idx = idx(lrg_idx);
    sm_idx = setdiff(idx, lrg_idx);
    for ii=lrg_idx
        this_comps = Clust.comps{ii}(ceil(rand(num_to_draw,1) .* ...
                         Clust.num_comps(ii)));
        comp_list = [comp_list; this_comps];
    end
    comp_list = [comp_list; cell2mat(Clust.comps(sm_idx))];
    pgs_to_check = unique(Comps.pg(comp_list));
else
    %use the component indices passed
    if size(comp_list,2) > 1
        comp_list = comp_list';
    end
    pgs_to_check = unique(Comps.pg(comp_list));
    num_to_draw = -Inf;
    for pp = pgs_to_check'
        num_to_draw = max(num_to_draw, sum(Comps.pg(comp_list) == pp));
    end
end

if display_avg
    num_to_draw = num_to_draw + 2;  %1 for the avg, and 1 for blankspace
end

%create a cell arrray to hold the component images
M = cell(num_clusts,num_to_draw);
num_added = zeros(num_clusts,1);
size_vals = zeros(num_clusts,num_to_draw,2);

%add the cluster average images if required (along with some blankspace)
if display_avg
    for ii=1:num_clusts
        M{ii, 1} = Clust.avg{idx(ii)};
        M{ii, 2} = zeros(1,col_pix_border);
        num_added(ii) = 2;
        size_vals(ii,1,:) = reshape(size(M{ii, 1}),1,1,2);
        size_vals(ii,2,:) = reshape(size(M{ii, 2}),1,1,2);
    end
end

%create a mapping for cluster indices into rows of M
map = zeros(Clust.num,1);
map(idx) = 1:num_clusts;

%extract the required components from each page and add them to M in the
%appropriate position
imgs = get_comp_imgs(Comps, comp_list);
for ii=1:length(comp_list)
    row = map(Comps.clust(comp_list(ii)));
    num_added(row) = num_added(row) + 1;
    M{row, num_added(row)} = imgs{ii};
    size_vals(row,num_added(row),:) = reshape(size(M{row, num_added(row)}), ...
                                      1,1,2);
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
