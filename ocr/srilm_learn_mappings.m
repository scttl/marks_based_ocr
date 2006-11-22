function [best_map,best_score,best_ppl] = srilm_learn_mappings(Clust, ...
         Comps, Syms, Lines, varargin)
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
%   best_map returned is a vector listing which character each cluster index 
%   maps to.
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
% $Id: srilm_learn_mappings.m,v 1.2 2006-11-22 17:09:30 scottl Exp $
%
% REVISION HISTORY
% $Log: srilm_learn_mappings.m,v $
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
dfs_ordering = true;  


unknown_char = '<unk>'; %symbol to use for oov chars

srilm_space_map = '_';  %symbol to use for space characters

%other parameters to pass to SRI-LM's ngram utility
other_params = '';


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
%convert the space character in the cluster so we can get proper SRILM counts
idx = find(strcmp(Clust.truth_label, ' '));
if ~isempty(idx)
    Clust.truth_label{idx} = srilm_space_map;
end

%sort the clusters by frequency
[Clust, Comps] = sort_clusters(Clust, Comps);

%get shortlist of mappings for each cluster
map = init_mappings(Clust, Syms);
fprintf('%.2fs: mappings initialized\n', toc);

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
b_map = [];
b_score = -Inf;
curr_map = unknown_sym_id .* ones(length(map),1);
for ii=1:Clust.num
    if length(map{ii}) == 1
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
    srilm_lm_close(pipe, out_file);
    fprintf('%.2fs: mapping complete\n', toc);
    return;
end

%now iterate over mappings in order until the optimal one is found
leaf_node = order(end);
if ~dfs_ordering
    %for the breadth-first search, we must first establish the score of the
    %left-most path (for comparison purposes)
    b_char_seq = char_seq;
    for ii=1:Clust.num
        best_map(ii) = map{ii}(1);
        b_char_seq(pos_idx{ii}) = Syms.val(b_map(ii));
    end
    [best_score, best_ppl] = srilm_lm_score(pipe, out_fid, ...
                       char(b_char_seq)');
end

list = [0, NaN, curr_score, curr_ppl];
while ~isempty(list)
    if dfs_ordering
        curr_idx = list(end,1);
        curr_char = list(end,2);
        curr_score = list(end,3);
        curr_ppl = list(end,4);
        list = list(1:end-1,:);
    else
        curr_idx = list(1,1);
        curr_char = list(1,2);
        curr_score = list(1,3);
        curr_ppl = list(1,4);
        list = list(2:end,:);
    end
    [list, curr_map, best_map, best_score, best_ppl] = iter_bb_solve_map(...
         list, order, pos_idx, map, char_seq, leaf_node, Syms, pipe,out_fid, ...
         curr_map, curr_idx, curr_char, curr_score, curr_ppl, ...
         best_map, best_score, best_ppl, unknown_char, unknown_sym_id);
end
srilm_lm_close(pipe, out_file);
fprintf('%.2fs: mapping complete\n', toc);




% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this is an iterative DFS that uses branching and bounding to prune non-optimal
%score paths
% 
function [list, curr_map, b_map, b_score, b_ppl] = iter_bb_solve_map(...
         list, order, pos_idx, map, char_seq, leaf_node, Syms,pipe, out_fid, ...
         curr_map, curr_idx, curr_char, curr_score, curr_ppl, ...
         b_map, b_score, b_ppl, unknown_char, unknown_sym_id)

if curr_idx > 0 && order(curr_idx) == leaf_node
    %compare scores and potentially update the lower bound best score and map
    fprintf('at leaf.  Best score thus far is %.2f\n', b_score);
    if curr_score > b_score
        %new best score
        b_score = curr_score;
        b_ppl = curr_ppl;
        b_map = curr_map;
    end
else
    %iterate through each of this nodes children to determine whether they are
    %worth exploring, and if so, calculate the best lower bound score using
    %them.
    if curr_idx > 0  %don't map the dummy root node
        curr_map(order(curr_idx)) = curr_char;
        %char_seq(pos_idx{order(curr_idx)}) = Syms.val(curr_char);
    end
    for ii=length(map{order(curr_idx+1)}):-1:1
    %for ii=1:length(map{order(curr_idx+1)}) %used for BFS ordering
        curr_child = map{order(curr_idx+1)}(ii);
        curr_map(order(curr_idx+1)) = curr_child;
        curr_map(order(curr_idx+2:end)) = unknown_sym_id;
        char_seq(pos_idx{order(curr_idx+1)}) = Syms.val(curr_child);
        for jj=curr_idx+2:length(order)
            char_seq(pos_idx{order(jj)}) = {unknown_char'};
        end
        [curr_score,curr_ppl] = srilm_lm_score(pipe, out_fid, ...
        char(char_seq)');
        if curr_score > b_score
            %branch case.  Add this depth, char, and score to the end of list
            list = [list; curr_idx+1, curr_child, curr_score, curr_ppl];
        else
            fprintf('bound at depth %d\n', curr_idx+1);
        end
    end
end
