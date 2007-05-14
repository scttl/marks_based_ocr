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
% $Id: positional_learn_mappings.m,v 1.6 2007-05-14 23:17:12 scottl Exp $
%
% REVISION HISTORY
% $Log: positional_learn_mappings.m,v $
% Revision 1.6  2007-05-14 23:17:12  scottl
% spelling typo correction in comments.
%
% Revision 1.5  2007-02-05 22:13:58  scottl
% implemented prior counts.
%
% Revision 1.4  2007-02-01 18:26:41  scottl
% updates to handle weighting.
%
% Revision 1.3  2007-01-29 03:23:39  scottl
% added optional weighting component based on word frequency.
%
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

%how much emphasis should we place on penalizing the distances by weighting
%the positions of each word length based on its expected (symbol) frequency.
%This number should be a probability (between 0 and 1) that will define the
%mixing proportion for the weights (the remainder will be distributed
%uniformly amongst all word lengths)
weight_proportion = 0;

%if assigning a weight matrix, should each symbol be adjusted by its particular
%word-length frequency, or by the overall word-length frequencies found in the
%document?
weight_per_symbol = false;

%how many prior counts should we add to each cluster before estimating the
%closest match?
prior_counts = 0;


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
    error('clusters must have position counts.  See create_cluster_dictionary');
end

%for each cluster, determine the distance to each symbol in symbol
%position feature space.
clust_pts = cell2mat(Clust.pos_count);
sym_pts = cell2mat(Syms.pos_count);
err_dist = zeros(Clust.num, Syms.num);
max_word_len = length(Syms.pos_count);
num_cols = size(sym_pts,2);
e_idx = cumsum(1:max_word_len);
s_idx = [1, e_idx(1:end-1) + 1];

if prior_counts > 0
    clust_pts = clust_pts .* cell2mat(Clust.pos_norms);
    raw_sym_pts = sym_pts .* cell2mat(Syms.pos_norms);
    for ii=unique(Clust.class)'
        sym_idx = find(Syms.class == ii);
        clust_idx = find(Clust.class == ii);
        for jj=1:max_word_len
            pri_pts = sum(raw_sym_pts(sym_idx,s_idx(jj):e_idx(jj)),1);
            norm = sum(pri_pts);
            norm(norm == 0) = 1;  %prevent divide by 0
            clp = clust_pts(clust_idx,s_idx(jj):e_idx(jj)) + ...
                  repmat((prior_counts .* (pri_pts ./ norm)), ...
                         length(clust_idx),1);

            %renormalize the now smoothed counts
            norms = sum(clp,2);
            norms(norms == 0) = 1;  %to prevent dividing by 0
            clust_pts(clust_idx,s_idx(jj):e_idx(jj)) = clp ./ ...
                                repmat(norms, 1, size(clp,2));
        end
    end
end

if keep_best == 0 || keep_best > Syms.num
    keep_best = Syms.num;
end
if weight_proportion > 0
    weights = sym_pts .* cell2mat(Syms.pos_norms);
    if weight_per_symbol
        norms = sum(weights,2);
        %to prevent division by 0 and ensure uniformly distributed weights 
        %where no counts are seen, we adjust these weights appropriately
        zero_idx = norms == 0;
        weights(zero_idx,s_idx) = 1;
        norms(zero_idx) = max_word_len; 
        for ii=1:length(s_idx)
            this_sum = sum(weights(:,s_idx(ii):e_idx(ii)),2);
            weights(:,s_idx(ii):e_idx(ii)) = repmat(this_sum ./ norms, 1, ...
                                             e_idx(ii)-s_idx(ii)+1);
        end
    else
        norm = sum(weights(:));
        freq_vector = sum(weights,1);
        for ii=1:max_word_len
            this_sum = sum(freq_vector(s_idx(ii):e_idx(ii)));
            weights(:,s_idx(ii):e_idx(ii)) = this_sum ./ norm;
        end
    end
else
    weights = 0;
end
%mix the weighted distribution with a uniform according to proportions
W = (weight_proportion .* weights) + ...
    ((1-weight_proportion) .* (max_word_len / num_cols) .* ones(size(sym_pts)));

order = cell(Clust.num, 1);
score = cell(Clust.num, 1);
for ii=1:Clust.num
    if strcmp(dist_metric, 'euc')
        err_dist(ii,:) = sum(W .* ...
                 (repmat(clust_pts(ii,:),Syms.num,1) - sym_pts).^2,2)';
    else
        err_dist(ii,:) = sum(W .* ...
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
