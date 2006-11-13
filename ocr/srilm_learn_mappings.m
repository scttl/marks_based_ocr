function [best_map,score] = srilm_learn_mappings(Clust, Comps, Syms, Lines, ...
                            varargin)
% SRILM_LEARN_MAPPINGS  Determine most likely cluster to symbol mappings
%
%   [map,score] = SRILM_LEARN_MAPPINGS(Clust,Comps,Syms,Lines,[var1, val1]...)
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
%   map returned is a vector listing which character each cluster index 
%   maps to.
%
%   score gives the overall log probability of the map learned when applied
%   over the sequence of all the lines
%
%   NOTE: the SRILM toolkit must be installed and the ngram tool must be in the 
%   users path


% CVS INFO %
%%%%%%%%%%%%
% $Id: srilm_learn_mappings.m,v 1.1 2006-11-13 18:13:05 scottl Exp $
%
% REVISION HISTORY
% $Log: srilm_learn_mappings.m,v $
% Revision 1.1  2006-11-13 18:13:05  scottl
% initial revision
%


% GLOBAL VARS %
%%%%%%%%%%%%%%%
%note that these are used to make function calls more efficient
%they will be cleared upon function termination.
global list order map leaf_node unknown_sym_id unknown_char pos_idx ...
       char_seq b_map srilm_seq_file b_score eol_seq eol_id num_start ...
       next_idx prev_idx order;


% LOCAL VARS %
%%%%%%%%%%%%%%
%by default we do a depth-first search through the mapping candidates.  Set
%this to false to do a breadth-first search.
dfs_ordering = true;  

srilm_seq_file = '/tmp/tmp_srilm.txt';

unknown_char = '*'; %character to use for oov chars


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

%add the newline, space, and unknown characters to the symbols
eol_sym_id = Syms.num + 1;
eol_clust_id = Clust.num + 1;
Syms.val(eol_sym_id) = 10; %newline is char 10
map{eol_clust_id} = eol_sym_id;
unknown_sym_id = Syms.num + 2;
unknown_clust_id = Clust.num + 2;
Syms.val(unknown_sym_id) = double(unknown_char);
map{unknown_clust_id} = unknown_sym_id;
space_sym_id = Syms.num + 3;
space_clust_id = Clust.num + 3;
Syms.val(space_sym_id) = 32; %space is char 32
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
leaf_node = order(end);

%now iterate over mappings in order until the optimal one is found
char_seq = char(Syms.val(curr_map(eol_seq)));
fid = fopen(srilm_seq_file, 'w');
fprintf(fid, '%s', char_seq);
fclose(fid);
curr_score = srilm_lm_score(Syms.srilm_file, srilm_seq_file);
if ~dfs_ordering
    %for the breadth-first search, we must first establish the score of the
    %left-most path (for comparison purposes)
    b_map = curr_map;
    b_char_seq = char_seq;
    for ii=1:Clust.num
        b_map(ii) = map{ii}(1);
        b_char_seq(pos_idx{ii}) = char(Syms.val(b_map(ii)));
    end
    fid = fopen(srilm_seq_file, 'w');
    fprintf(fid, '%s', b_char_seq);
    fclose(fid);
    b_score = srilm_lm_score(Syms.srilm_file, srilm_seq_file);
end

list = [0, NaN, curr_score];
while ~isempty(list)
    if dfs_ordering
        curr_node = list(end,1);
        curr_char = list(end,2);
        curr_score = list(end,3);
        list = list(1:end-1,:);
    else
        curr_node = list(1,1);
        curr_char = list(1,2);
        curr_score = list(1,3);
        list = list(2:end,:);
    end
    curr_map = iter_bb_solve_map(Syms, curr_map, curr_node, curr_char, curr_score);
end
fprintf('%.2fs: mapping complete\n', toc);
best_map = b_map;
score = b_score;

clear map leaf_node b_map b_score eol_seq eol_id num_start ...
      next_idx prev_idx order;




% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this is an iterative DFS that uses branching and bounding to prune non-optimal
%score paths
function curr_map = iter_bb_solve_map(Syms, curr_map, curr_idx, curr_char, curr_score)

global list order leaf_node map b_map b_score unknown_sym_id unknown_char ...
       char_seq pos_idx srilm_seq_file;

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
    %for ii=1:length(map{order(curr_idx+1)}) %used for BFS ordering
        curr_child = map{order(curr_idx+1)}(ii);
        curr_map(order(curr_idx+1)) = curr_child;
        curr_map(order(curr_idx+2:end)) = unknown_sym_id;
        char_seq(pos_idx{order(curr_idx+1)}) = char(Syms.val(curr_child));
        for jj=curr_idx+2:length(order)
            char_seq(pos_idx{order(jj)}) = unknown_char;
        end
        fid = fopen(srilm_seq_file, 'w');
        fprintf(fid, '%s', char_seq);
        fclose(fid);
        curr_score = srilm_lm_score(Syms.srilm_file, srilm_seq_file);
        if curr_score > b_score
            %branch case.  Add this depth, char, and score to the end of list
            list = [list; curr_idx+1, curr_child, curr_score];
        else
            fprintf('bound at depth %d\n', curr_idx+1);
        end
    end
end
