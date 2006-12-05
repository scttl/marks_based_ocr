function [best_map,best_score,best_ppl] = srilm_learn_mappings(Clust, ...
         Comps, Syms, Lines, map, varargin)
% SRILM_LEARN_MAPPINGS  Determine most likely cluster to symbol mappings
%
%   [best_map, best_score, best_ppl] = SRILM_LEARN_MAPPINGS(Clust, Comps, ...
%                                      Syms, Lines, [var1,val1]...)
%   This function attempts to learn a good mapping from each cluster to a
%   single character symbol by finding the maximum-likelihood score mapping
%   when procesing the image lines in Lines.  This function makes use of the 
%   SRILM tookit to take care of most of the heavy lifting
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
%   map should be a cell array the same length as Clust and should list for
%   each cluster, the set of possible indices from Syms to which it best
%   matches (i.e. which character symbols it most closely resembles).
%
%   best_score gives the overall log probability of the best map learned when 
%   applied over the sequence of all the lines
%
%   best_ppl gives the overall sequence perplexity (akin to a branching factor)
%
%   NOTE: the SRILM toolkit must be installed and the ngram tool must be in the 
%   users path


% CVS INFO %
%%%%%%%%%%%%
% $Id: srilm_learn_mappings.m,v 1.3 2006-12-05 16:05:23 scottl Exp $
%
% REVISION HISTORY
% $Log: srilm_learn_mappings.m,v $
% Revision 1.3  2006-12-05 16:05:23  scottl
% major re-write fixing bugs and separating map setup from this function.
%
% Revision 1.2  2006-11-22 17:09:30  scottl
% rewritten based on new scoring methods, got rid of global variables
% since they shouldn't have much an impact on performance.
%
% Revision 1.1  2006-11-13 18:13:05  scottl
% initial revision
%


% GLOBAL VARS %
%%%%%%%%%%%%%%%


% LOCAL VARS %
%%%%%%%%%%%%%%
%by default we do a depth-first search through the mapping candidates.  Set
%this to false to do a breadth-first search.
dfs_ordering = true
;  

unknown_char = '<unk>'; %symbol to use for oov chars

srilm_space_map = '_';  %symbol to use for space characters

%other parameters to pass to SRI-LM's ngram utility
other_params = '';

%this value controls how deep in the tree we will explore.  Set to Inf, all
%levels will be explored
bound_depth = 32;

% CODE START %
%%%%%%%%%%%%%%
tic;
if nargin < 5
    error('incorrect number of arguments specified!');
elseif nargin > 5
    process_optional_args(varargin{:});
end

%some basic sanity checking
if ~isfield(Syms, 'img') || isempty(Syms.img)
    error('we require template symbol images!');
end
if ~isfield(Syms, 'corpus_files') || isempty(Syms.corpus_files)
    error('we require symbol and word counts from a text corpus!');
end
if ~isfield(Clust, 'model_spaces') || ~Clust.model_spaces
    error('space characters required!');
end

%convert the space character in the cluster so we can get proper SRILM counts
idx = find(strcmp(Clust.truth_label, ' '));
if ~isempty(idx)
    Clust.truth_label{idx} = srilm_space_map;
end

%get the sequence of cluster id's for each line we wish to solve
seq = get_cluster_seq(Comps, 1:Lines.num);
fprintf('%.2fs: cluster sequence found\n', toc);

%ensure that multi-char symbols are written vertically with spaces between
%characters and add the newline, space, and unknown characters to the symbols
for ii=1:Syms.num
    if size(Syms.val{ii},2) > 1
        Syms.val{ii} = regexprep(Syms.val{ii}, '(.)', '$1 ')';
        %remove the additional space at the end
        Syms.val{ii} = Syms.val{ii}(1:end-1);
    end
end
eol_sym_id = Syms.num + 1;
eol_clust_id = Clust.num + 1;
Syms.val{eol_sym_id} = char(10); %newline
map{eol_clust_id} = eol_sym_id;
unknown_sym_id = Syms.num + 2;
unknown_clust_id = Clust.num + 2;
Syms.val{unknown_sym_id} = unknown_char';
map{unknown_clust_id} = unknown_sym_id;
space_sym_id = Syms.num + 3;
space_clust_id = Clust.num + 3;
Syms.val{space_sym_id} = ' ';
map{space_clust_id} = space_sym_id;


%convert the line sequences to a single newline delimited sequence with spaces 
%between each character
eol_seq = [];
for ii=1:length(seq)
    nseq = zeros(1,(2*length(seq{ii}))-1);
    nseq(1,1:2:size(nseq,2)) = seq{ii}';
    nseq(1,2:2:size(nseq,2)) = space_clust_id;
    eol_seq = [eol_seq, nseq, eol_clust_id];
end

%now store for each cluster where in the sequence it is
pos_idx = cell(Clust.num,1);
for ii=1:Clust.num
    pos_idx{ii} = find(eol_seq == ii);
end
fprintf('%.2fs: cluster sequence positions stored\n', toc);

%determine an appropriate ordering, trimming nodes with a single mapping
order = 1:Clust.num;
curr_map = unknown_sym_id .* ones(length(map),1);
for ii=1:Clust.num
    if ii > bound_depth
        curr_map(ii) = unknown_sym_id;
        pos = ii - (Clust.num - length(order));
        order = order([1:pos-1,pos+1:end]);
    elseif length(map{ii}) == 1
        curr_map(ii) = map{ii}(1);
        pos = ii - (Clust.num - length(order));
        order = order([1:pos-1,pos+1:end]);
    end
end
curr_map(eol_clust_id) = eol_sym_id;
curr_map(unknown_clust_id) = unknown_sym_id;
curr_map(space_clust_id) = space_sym_id;
char_seq = cell(length(eol_seq),1);
char_seq(:) = Syms.val(curr_map(eol_seq));

[pipe, out_file, out_fid] = srilm_lm_open(Syms.srilm_file);
[curr_score,curr_ppl] = srilm_lm_score(pipe,out_fid,char(char_seq)');
best_map = curr_map;
best_score = -Inf;
best_ppl = -Inf;

if isempty(order)
    %optimal mapping already found since only 1 match for each
    best_score = curr_score;
    best_ppl = curr_ppl;
    best_map = best_map(1:end-3);  %remove the eol, unknown and space maps
    srilm_lm_close(pipe, out_file);
    fprintf('%.2fs: optimal mapping complete\n', toc);
    return;
end

%now iterate over mappings in order until the optimal one is found
leaf_node = order(end);
list = [0, NaN];
if ~dfs_ordering
    %for the breadth-first search, we must first establish the score of the
    %left-most path (for comparison purposes)
    %path taken by exploring the child that gives the lowest score at each
    %level (based on previous mappings)
    b_char_seq = char_seq;
    for ii=order
        %best_map(ii) = map{ii}(1);
        bs = -Inf;
        bi = 0;
        for jj=1:length(map{ii})
            b_char_seq(pos_idx{ii}) = Syms.val(map{ii}(jj));
            ns = srilm_lm_score(pipe, out_fid, char(b_char_seq)');
            if ns > bs
                bs = ns;
                bi = jj;
            end
        end
        best_map(ii) = map{ii}(bi);
        b_char_seq(pos_idx{ii}) = Syms.val(best_map(ii));
    end
    [best_score, best_ppl] = srilm_lm_score(pipe, out_fid, ...
                       char(b_char_seq)');
    fprintf('initial left-most path score: %.2f\n', best_score);
    list = {{0, NaN, []}};
end

while ~isempty(list)
    if dfs_ordering
        curr_idx = list(end,1);
        curr_char = list(end,2);
        list = list(1:end-1,:);
        [list,char_seq,curr_map,best_map,best_score,best_ppl] = ...
         iter_bb_solve_map(list, order, pos_idx, map, char_seq, leaf_node, ...
         Syms, pipe,out_fid, curr_map, curr_idx, curr_char, ...
         best_map, best_score, best_ppl, unknown_char, unknown_sym_id);
    else
        curr_idx = list{1}{1};
        curr_char = list{1}{2};
        prev_map = list{1}{3};
        for ii=1:length(prev_map)
            if curr_map(ii) ~= prev_map(ii)
                curr_map(ii) = prev_map(ii);
                char_seq(pos_idx{ii}) = Syms.val(prev_map(ii));
            end
        end
        list = list(2:end,:);
        [list,char_seq,curr_map,best_map,best_score,best_ppl] = ...
         iter_bb_solve_bfs_map(list, order, pos_idx, map, char_seq, ...
         leaf_node, Syms, pipe,out_fid, curr_map, curr_idx, curr_char, ...
         best_map, best_score, best_ppl, unknown_char, unknown_sym_id);
    end
end
best_map = best_map(1:end-3);  %remove the eol, unknown and space maps
srilm_lm_close(pipe, out_file);
fprintf('%.2fs: mapping complete\n', toc);




% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this is an iterative DFS that uses branching and bounding to prune non-optimal
%score paths
% 
function [list,char_seq,curr_map,b_map,b_score,b_ppl] = iter_bb_solve_map(...
         list, order, pos_idx, map, char_seq, leaf_node, Syms,pipe, out_fid, ...
         curr_map, curr_idx, curr_char, ...
         b_map, b_score, b_ppl, unknown_char, unknown_sym_id)

if curr_idx > 0
    %replace the curr_char in the sequence
    curr_map(order(curr_idx)) = curr_char;
    char_seq(pos_idx{order(curr_idx)}) = Syms.val(curr_char);

    if curr_char == unknown_sym_id
        %this just resets the unknown symbol at this depth, so we can explore
        %earlier levels
        return;
    end

    %re-calucuate the score
    [curr_score,curr_ppl] = srilm_lm_score(pipe, out_fid, char(char_seq)');
    fprintf('depth: %d, score: %.2f, best: %.2f, list length %d\n', ...
                curr_idx, curr_score, b_score, length(list));

    if curr_score < b_score
        %can chop right here
        fprintf('chop at depth %d.  Curr %.2f Best %.2f List length %d\n', ...
                curr_idx, curr_score, b_score, length(list));
        return
    end

    if order(curr_idx) == leaf_node
        %compare scores and potentially update the lower bound best score and 
        %map
        fprintf('at leaf.  New best score is %.2f\n', curr_score);
        b_score = curr_score;
        b_ppl = curr_ppl;
        b_map = curr_map;
        return;
    end
end

%iterate through each of this nodes children to add them to the list.  We
%prepend the list with the unknown symbol so that it gets replaced in the
%sequence after we've explored all the mappings
list = [list; curr_idx+1, unknown_sym_id];
for ii=length(map{order(curr_idx+1)}):-1:1
    curr_child = map{order(curr_idx+1)}(ii);
    list = [list; curr_idx+1, curr_child];
end


%this is an iterative BFS that uses branching and bounding to prune non-optimal
%score paths
% 
function [list,char_seq,curr_map,b_map,b_score,b_ppl] = iter_bb_solve_bfs_map(...
         list, order, pos_idx, map, char_seq, leaf_node, Syms,pipe, out_fid, ...
         curr_map, curr_idx, curr_char, ...
         b_map, b_score, b_ppl, unknown_char, unknown_sym_id)

if curr_idx > 0
    %replace the curr_char (and previous) in the sequence
    curr_map(order(curr_idx)) = curr_char;
    char_seq(pos_idx{order(curr_idx)}) = Syms.val(curr_char);

    %re-calucuate the score
    [curr_score,curr_ppl] = srilm_lm_score(pipe, out_fid, char(char_seq)');
    fprintf('depth: %d, score: %.2f, best: %.2f, list length %d\n', ...
                curr_idx, curr_score, b_score, length(list));

    if curr_score < b_score
        %can chop right here
        fprintf('chop at depth %d.  Curr %.2f Best %.2f List length %d\n', ...
                curr_idx, curr_score, b_score, length(list));
        return
    end

    if order(curr_idx) == leaf_node
        %compare scores and potentially update the lower bound best score and 
        %map
        fprintf('at leaf.  New best score is %.2f\n', curr_score);
        b_score = curr_score;
        b_ppl = curr_ppl;
        b_map = curr_map;
        return;
    end

    %iterate through each of this nodes children to add them to the list.
    for ii=1:length(map{order(curr_idx+1)})
        curr_child = map{order(curr_idx+1)}(ii);
        list = [list; {{curr_idx+1, curr_child, curr_map(1:order(curr_idx))}}];
    end
else
    for ii=1:length(map{order(curr_idx+1)})
        curr_child = map{order(curr_idx+1)}(ii);
        list = [list; {{curr_idx+1, curr_child, []}}];
    end
end

