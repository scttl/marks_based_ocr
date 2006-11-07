function [best_map,score] = learn_mappings(Clust, Comps, Syms, Lines, varargin)
% LEARN_MAPPINGS  Determine most likely cluster to symbol mappings
%
%   [map,score] = LEARN_MAPPINGS(Clust, Comps, Syms, Lines, [var1, val1]...)
%   This function attempts to learn a good mapping from each cluster to a
%   single character symbol by finding the maximum-likelihood score mapping
%   when procesing the image lines in Lines.
%   
%   Clust should be a struct like that returned from cluster_comps()
%
%   Comps should be a struct like that returned from get_comps()
%
%   Syms should be a struct like that returned from create_alphabet.  We
%   require that template images have been created for each symbol too (via
%   generate_templates)
%
%   Lines should be a struct like that returned from get_lines()
%
%   map returned is a vector listing which character each cluster index 
%   maps to.
%
%   score gives the overall log probability of the map learned when applied
%   over the sequence of all the lines


% CVS INFO %
%%%%%%%%%%%%
% $Id: learn_mappings.m,v 1.1 2006-11-07 02:55:25 scottl Exp $
%
% REVISION HISTORY
% $Log: learn_mappings.m,v $
% Revision 1.1  2006-11-07 02:55:25  scottl
% initial check-in.
%


% GLOBAL VARS %
%%%%%%%%%%%%%%%
%note that these are used to make function calls more efficient
%they will be cleared upon function termination.
global list num_chars map lg_bg lg_start leaf_node b_map b_score eol_seq ...
       eol_id num_start next_idx prev_idx order;


% LOCAL VARS %
%%%%%%%%%%%%%%
eol_id = 0;   %this is used to create a dummy cluster representing end-of-line


% CODE START %
%%%%%%%%%%%%%%
tic;
if nargin < 3
    error('incorrect number of arguments specified!');
elseif nargin > 3
    process_optional_args(varargin{:});
end

if ~isfield(Syms, 'img') || isempty(Syms.img)
    error('we require template symbol images!');
end
if ~isfield(Syms, 'corpus_files') || isempty(Syms.corpus_files)
    error('we require symbol and word counts from a text corpus!');
end

%check that spaces have been added (and add them if they haven't)
if ~isfield(Clust, 'model_space') || ~Clust.model_space
    warning('MBOCR:ModelSpace', 'Adding space model clusters and comps\n');
    [Clust, Comps] = add_space_model(Clust, Comps);
end

%sort the clusters by frequency
[Clust, Comps] = sort_clusters(Clust, Comps);

%get shortlist of mappings for each cluster
map = init_mappings(Clust, Syms);
fprintf('%.2fs: mappings initialized\n', toc);

%get the sequence of cluster id's for each line we wish to solve
seq = get_cluster_seq(Comps, 1:Lines.num);
fprintf('%.2fs: cluster sequence found\n', toc);


%get the log probabilities of our bigram and start character distributions
num_chars = Syms.num;
lg_start = log(double(Syms.first_count) ./ ...
               repmat(sum(double(Syms.first_count)), Syms.num, 1));
lg_bg = log(Syms.bigram ./ repmat(sum(Syms.bigram, 2), 1, Syms.num));
%augment them with their maximum, for unknown variable calculations
lg_start(num_chars + 1) = max(lg_start);
lg_bg(num_chars+1,:) = max(lg_bg,[],1);
lg_bg(:, num_chars+1) = max(lg_bg,[],2);

%convert the line sequences to a single EOL delimited sequence
eol_seq = [eol_id];
for ii=1:length(seq)
    eol_seq = [eol_seq; seq{ii}; eol_id];
end

%now store for each cluster where in the sequence it is
num_start = zeros(Clust.num);
next_idx = cell(Clust.num,1);
prev_idx = cell(Clust.num,1);
for ii=1:Clust.num
    idx = find(eol_seq == ii);

    pidx = idx - 1;
    prev_idx{ii} = pidx(eol_seq(pidx) ~= eol_id);

    %to prevent double-counting self-transitions for the clust_id, we exclude 
    %them during our next transition counts
    nidx = idx + 1;
    next_idx{ii} = nidx(eol_seq(nidx) ~= eol_id & eol_seq(nidx) ~= ii);

    num_start(ii) = length(pidx) - length(prev_idx{ii});
end
fprintf('%.2fs: cluster sequence neighbours stored\n', toc);

%determine an appropriate ordering
order = 1:Clust.num

%now iterate over mappings in order until the optimal one is found
leaf_node = Clust.num;
b_map = [];
b_score = -Inf;
curr_map = (num_chars+1) .* ones(Clust.num,1);
list = [0, NaN, ...
       sum(score_sequence(clust_to_char_map(seq, curr_map), lg_start, lg_bg))];
while ~isempty(list)
    curr_node = list(end,1);
    curr_char = list(end,2);
    curr_score = list(end,3);
    list = list(1:end-1,:);
    curr_map = iter_bb_solve_map(curr_map, curr_node, curr_char, curr_score);
end
fprintf('%.2fs: mapping complete\n', toc);
best_map = b_map;
score = b_score;

clear map lg_bg lg_start leaf_node b_map b_score eol_seq eol_id num_start ...
      next_idx prev_idx order;




% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this is an iterative DFS that uses branching and bounding to prune non-optimal
%score paths
function curr_map = iter_bb_solve_map(curr_map, curr_idx, curr_char, curr_score)

global list order num_chars leaf_node map b_map b_score;

if curr_idx > 0 && order(curr_idx) == leaf_node
    %compare scores and potentially update the lower bound best score and map
    fprintf('at leaf.  Best score thus far is %.2f\n', b_score);
    if curr_score > b_score
        %new best score
        b_score = curr_score;
        b_map = curr_map;
    end
else
    %iterate through each of this nodes children to determine whether they are
    %worth exploring, and if so, calculate the best lower bound score using
    %them.
    if curr_idx > 0  %don't map the dummy root node
        curr_map(order(curr_idx)) = curr_char;
    end
    for ii=length(map{order(curr_idx+1)}):-1:1
        curr_child = map{order(curr_idx+1)}(ii);
        curr_map(order(curr_idx+1)) = curr_child;
        curr_map(order(curr_idx+2:end)) = num_chars+1;
        delta = calc_delta(curr_map,order(curr_idx+1),curr_child,num_chars+1);
        if curr_score + delta > b_score
            %branch case.  Add this depth, char, and score to the end of list
            list = [list; curr_idx+1, curr_child, curr_score + delta];
        end
    end
end



%this recursive function is used to efficiently find the best mapping using a
%branch and bound approach to enumerating the search space.
%NOTE: it is not currently used since it cannot handle mappings of more than
%100 clusters (Matlab recursion limitation).  Even when resetting the built-in
%variable, the stack still overflows after about 400 or so recursive calls.
function bb_solve_map(curr_node, c_map, curr_score)

global num_chars leaf_node map b_map b_score;

if curr_node == leaf_node
    %compare scores and potentially update the lower bound best score and map
    fprintf('at leaf.  Best score thus far is %.2f\n', b_score);
    orig_char = c_map(leaf_node);
    for ii=1:length(map{leaf_node})
        curr_char = map{leaf_node}(ii);
        c_map(leaf_node) = curr_char;
        delta = calc_delta(c_map, leaf_node, curr_char, orig_char);
        if curr_score + delta > b_score
            %new best score
            b_score = curr_score + delta;
            b_map = c_map;
        end
    end
else
    %iterate through each of this nodes children to determine whether they are
    %worth exploring, and if so, calculate the best lower bound score using
    %them.
    for ii=1:length(map{curr_node+1})
        orig_char = c_map(curr_node+1);
        curr_char = map{curr_node+1}(ii);
        c_map(curr_node+1) = curr_char;
        delta = calc_delta(c_map, curr_node+1, curr_char, orig_char);
        curr_score = curr_score + delta;
        if curr_score > b_score
            %branch case
            bb_solve_map(curr_node+1, c_map, curr_score);
        end
    end
end



%this recursive function is used to enumerate all possible mappings for the
%best one
%NOTE: not currently used due to recursion and efficiency limitations.
%originally used to check correctness of branch and bound implementation
function solve_map(curr_node, c_map, curr_score)

global leaf_node map b_map b_score;

if curr_node == leaf_node
    %compare scores and potentially update the best map
    fprintf('at leaf.  Best score thus far is %.2f\n', b_score);
    orig_char = c_map(leaf_node);
    for ii=1:length(map{leaf_node})
        curr_char = map{leaf_node}(ii);
        c_map(leaf_node) = curr_char;
        delta = calc_delta(c_map, leaf_node, curr_char, orig_char);
        if curr_score + delta > b_score
            %new best score
            b_score = curr_score + delta;
            b_map = c_map;
        end
    end
else
    %iterate through each of this nodes children to find the best score using it
    for ii=1:length(map{curr_node+1})
        orig_char = c_map(curr_node+1);
        curr_char = map{curr_node+1}(ii);
        c_map(curr_node+1) = curr_char;
        delta = calc_delta(c_map, curr_node+1, curr_char, orig_char);
        curr_score = curr_score+delta;
        solve_map(curr_node+1, c_map, curr_score);
    end
end



%this function calculates the delta value for changing the cluster id passed
%from the first character id to the second one.
%NOTE: we assume the cluster sequence is a single vector, with lines separated
%by EOL or dummy cluster id (which should be passed in eol_id).  We also assume
%the cluster sequence both starts with, and ends with the EOL id.
%Finally, if the new or old chararacter id's are outside of the range of valid
%id's (given the size of the start and bigram tables, we treat them as unknowns,
%and calculate their score by choosing the best possible value.
function delta = calc_delta(map, clust_id, newchar_id, oldchar_id)

global num_chars eol_seq eol_id lg_start lg_bg num_start next_idx prev_idx;

if eol_id == clust_id
    error('trying to calculate delta for EOL id!');
end

if ~(1 <= newchar_id && newchar_id <= num_chars) && ...
   ~(1 <= oldchar_id && oldchar_id <= num_chars)
   %both characters unspecified, thus there is no change in the mapping
   delta = 0;
   return;
end

%since we can potentially transfer to clust_id consecutvely, we must ensure we
%use the correct mapping.  We also must fixup the mapping if newchar_id or
%oldchar_id are left unspecified.
map(map < 1 | map > num_chars | isnan(map)) = num_chars + 1;
new_map = map;
new_map(clust_id) = newchar_id;
old_map = map;
old_map(clust_id) = oldchar_id;
if ~(1 <= newchar_id && newchar_id <= num_chars)
    %newchar is unspecified
    newchar_id = num_chars + 1;
    new_map(clust_id) = newchar_id;
end
if ~(1 <= oldchar_id && oldchar_id <= num_chars)
    %oldchar is unspecified
    oldchar_id = num_chars + 1;
    old_map(clust_id) = oldchar_id;
end

%Add the score for the number of times we start a line with this cluster
delta = num_start(clust_id) * (lg_start(newchar_id) - lg_start(oldchar_id));

%now add the scores of transitions *to* this cluster
delta = delta + ...
        sum(lg_bg(new_map(eol_seq(prev_idx{clust_id})), newchar_id)) - ...
        sum(lg_bg(old_map(eol_seq(prev_idx{clust_id})), oldchar_id));

%finally, add the scores of valid transitions *from* this cluster
delta = delta + ....
        sum(lg_bg(newchar_id, new_map(eol_seq(next_idx{clust_id})))) - ...
        sum(lg_bg(oldchar_id, old_map(eol_seq(next_idx{clust_id}))));



%this function takes an input cell array of cluster id's, and a cluster to
%character mapping function, and returns a cell array of correspondung
%character id's.
function char_seq = clust_to_char_map(clust_seq, map)
if ~iscell(clust_seq)
    char_seq = map(clust_seq);
else
    char_seq = cell(size(clust_seq));
    for ii=1:size(clust_seq,1)
        for jj=1:size(clust_seq,2)
            char_seq{ii,jj} = map(clust_seq{ii,jj});
        end
    end
end
