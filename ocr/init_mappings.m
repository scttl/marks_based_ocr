function map = init_mappings(Clust, Syms, varargin)
% INIT_MAPPINGS  Determine an initial mapping for each cluster to a template img
%
%   map = INIT_MAPPINGS(Clust, Syms, [var1, val1]...)
%   Clust should be a struct like that returned from cluster_comps()
%
%   Syms should be a struct like that returned from create_alphabet.  We
%   require that template images have been created for each symbol too (via
%   generate_templates)
%
%   map returned is a cell array, with one row per cluster, each entry of which
%   contains a vector listing character indices that match this symbol within a
%   predefined threhold (which can be overridden).  We make sure that the
%   closest matching symbol is also included (based on a distance metric)
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: init_mappings.m,v 1.1 2006-10-29 17:28:29 scottl Exp $
%
% REVISION HISTORY
% $Log: init_mappings.m,v $
% Revision 1.1  2006-10-29 17:28:29  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%distance metric to use when matching
match_metric = 'hausdorff';

take_nearest = 1;  %number of closest matchings to take for each cluster

match_thresh = 1.5;  %any matches with a distance less than this are also taken

resize_method = 'nearest';  %method of resize for symbol images



% CODE START %
%%%%%%%%%%%%%%
tic;
if nargin < 2
    error('incorrect number of arguments specified!');
elseif nargin > 2
    process_optional_args(varargin{:});
end

if ~isfield(Syms, 'img') || isempty(Syms.img)
    error('we require template symbol images!');
end

map = cell(Clust.num,1);

cl_heights = zeros(Clust.num,1);
for ii=1:Clust.num
    cl_heights(ii) = size(Clust.avg{ii},1);
end
md_height = mode(cl_heights);
fprintf('%.2fs: determined modal cluster image height = %d\n', toc, md_height);

rsz_syms = cell(Syms.num,1);
for ii=1:Syms.num
    rsz_syms{ii} = imresize(Syms.img{ii}, md_height/size(Syms.img{ii},1), ...
                   resize_method);
end
fprintf('%.2fs: finished resizing symbol images\n', toc, md_height);

%loop through each cluster and determine which symbols map to it.
for ii=1:Clust.num
    if strcmp(match_metric, 'hausdorff')
        D = hausdorff_dist(Clust.avg{ii}, rsz_syms);
        [val,idx] = sort(D);
        num = sum(val <= match_thresh);
        if num < take_nearest
            num = take_nearest;
        end
        map{ii} = idx(1:num);
    end
    fprintf('%.2fs: found %d matches for cluster %d\n',toc,length(map{ii}),ii);
end


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
