function display_mappings(Clust, Syms, map, varargin)
% DISPLAY_MAPPINGS   Display the cluster mappings for each symbol
%
%   DISPLAY_MAPPINGS(Clust, Syms, map, [VAR1, VAL1]...)
%
%   Clust should be a struct that includes an avg field (see cluster_comps())
%
%   Syms should be a struct that contains template images in an img field (see
%   create_alphabet() and generate_templates())
%
%   map should be a cell array with one entry per cluster.  Each entry should
%   list a (possibly empty) vector of symbol indices, representing the
%   symbols that map to that cluster (within the threshold)
%
%   Additional variables specified in LOCAL VARS below, can be overridden by
%   passing the name of the variable to override as a string, followed by its
%   new value (and repeating for each such variable to override).
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: display_mappings.m,v 1.1 2006-12-04 19:21:08 scottl Exp $
%
% REVISION HISTORY
% $Log: display_mappings.m,v $
% Revision 1.1  2006-12-04 19:21:08  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%
row_pix_border = 10;
col_pix_border = 3;

%set this to true to display the symbol images that map to each cluster.
swap_map = false;

%set save_mappings to true to write the elements image to disk based on the 
%params below it
save_mappings = false;
global MOCR_PATH;  %make use of the globally defined MOCR_PATH variable
img_prefix = [MOCR_PATH, '/results/cluster_mappings'];
img_format = 'png';


% CODE START %
%%%%%%%%%%%%%%

if nargin < 3
    error('incorrect number of arguments specified!');
elseif nargin > 3
    process_optional_args(varargin{:});
end

if ~isfield(Syms, 'img') || isempty(Syms.img)
    error('Syms struct must have non-empty img field');
end

if size(map,1) ~= Clust.num
    error('mapping is of incorrect dimensions');
end

num_rows = Syms.num;
row_imgs = Syms.img;
col_imgs = Clust.avg;
%have to invert the mapping so we have a list of cluster id's for each row
%(which represents a symbol)
nmap = cell(Syms.num,1);
for ii=1:Clust.num
    for jj=1:length(map{ii})
        val = map{ii}(jj);
        nmap{val} = [nmap{val}; ii];
    end
end

if swap_map
    num_rows = Clust.num;
    row_imgs = Clust.avg;
    col_imgs = Syms.img;
    nmap = map;
end

num_to_draw = 0;
for ii=1:num_rows
    num_to_draw = max(num_to_draw, length(nmap{ii}));
end
num_to_draw = num_to_draw + 2;  %1 for the first img, and 1 for blankspace

%create a cell arrray to hold the mapping images
M = cell(num_rows,num_to_draw);
num_added = zeros(num_rows,1);
size_vals = zeros(num_rows,num_to_draw,2);

%add the template image (along with some blankspace)
for ii=1:num_rows
    M{ii, 1} = row_imgs{ii};
    M{ii, 2} = zeros(1,col_pix_border);
    num_added(ii) = 2;
    size_vals(ii,1,:) = reshape(size(M{ii, 1}),1,1,2);
    size_vals(ii,2,:) = reshape(size(M{ii, 2}),1,1,2);
end

%add images of each column in the appropriate position
for ii=1:num_rows
    num_col = length(nmap{ii});
    if num_col > 0
        M(ii, num_added(ii)+1:num_added(ii)+num_col) = col_imgs(nmap{ii});
        for jj=num_added(ii)+1:num_added(ii)+num_col
            size_vals(ii,jj,:) = reshape(size(M{ii,jj}), 1,1,2);
        end
        num_added(ii) = num_added(ii) + num_col;
    end
end

%convert M to an appropriately spaced image
col_w = max(max(size_vals(:,:,2)));
MM = zeros(sum(max(size_vals(:,:,1),[],2)) + num_rows*row_pix_border, ...
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
imtool(MM);

%save the image to disk if required.
if save_mappings
    fprintf('writing cluster mapping image to disk\n');
    imwrite(MM, [img_prefix, '.', img_format], img_format);
end
