function [order,score] = positional_learn_mappings(Clust, Syms, varargin)
% POSITIONAL_LEARN_MAPPINGS  Determine most likely cluster to symbol mappings
%
%   [order, score] = POSITIONAL_LEARN_MAPPINGS(Clust, Syms, [var1,val1]...)
%   This function attempts to learn a good mapping from each cluster to a
%   single character symbol by comparing positional statistics with characters 
%   taken from a large text corpus.
%   
%   Clust should be a struct like that returned from cluster_comps()
%
%   Syms should be a struct like that returned from create_alphabet()
%
%   order returned will be a cell array listing for each cluster the symbol 
%   index order with the closest match in the first column, and the farthest 
%   match in the last column
%
%   score gives the corresponding distance to each symbol
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: positional_learn_mappings.m,v 1.2 2007-01-18 19:14:45 scottl Exp $
%
% REVISION HISTORY
% $Log: positional_learn_mappings.m,v $
% Revision 1.2  2007-01-18 19:14:45  scottl
% changed order to a cell array so that it can differ in length for
% each cluster.
%
% Revision 1.1  2007-01-05 17:16:59  scottl
% initial check-in.
%


% GLOBAL VARS %
%%%%%%%%%%%%%%%


% LOCAL VARS %
%%%%%%%%%%%%%%
dist_metric = 'euc';  %other choice is 'manhattan'

%should we limit the symbols returned to the unique values only?
only_unique = true;

%should we stop including distances after a certain number?  Set to 0 to
%include all distances
keep_best = 0;



% CODE START %
%%%%%%%%%%%%%%
tic;
if nargin < 2
    error('incorrect number of arguments specified!');
elseif nargin > 2
    process_optional_args(varargin{:});
end

%some basic sanity checking
if ~isfield(Syms, 'pos_count') || isempty(Syms.pos_count)
    error('we require symbol position counts from a text corpus!');
end
if ~isfield(Clust, 'pos_count') || isempty(Clust.pos_count)
    error('clusters must have position counts.  See create_cluster_dictionar');
end

%for each cluster, determine the distance to each symbol in symbol
%position feature space.
clust_pts = cell2mat(Clust.pos_count);
sym_pts = cell2mat(Syms.pos_count);
err_dist = zeros(Clust.num, Syms.num);
num_cols = Syms.num;
if keep_best == 0 || keep_best > Syms.num
    keep_best = Syms.num;
end
order = cell(Clust.num, 1);
score = cell(Clust.num, 1);
for ii=1:Clust.num
    if strcmp(dist_metric, 'euc')
        err_dist(ii,:) = sum(...
                 (repmat(clust_pts(ii,:),Syms.num,1) - sym_pts).^2,2)';
    else
        err_dist(ii,:) = sum(...
                 abs(repmat(clust_pts(ii,:),Syms.num,1) - sym_pts),2)';
    end

    %sort the distances to determine the optimal mapping for this cluster
    [val,idx] = sort(err_dist(ii,:));
    %@@@if only_unique
        %trim the indices to ensure that there are no repeat symbol values
        %@@jj = 1;
    %@@end
    order{ii} = idx(1:keep_best);
    score{ii} = val(1:keep_best);
end
fprintf('%.2fs: ordering complete\n', toc);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
