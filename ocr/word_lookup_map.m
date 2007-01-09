function map = word_lookup_map(Clust, Comps, Syms, varargin)
% WORD_LOOKUP_MAP  Determine best cluster mappings using word match ratios
%
%   MAP = word_lookup_map(CLUST, COMPS, SYMS, [var1,val1]...)
%   This function exploits word information taken from SYMS to determine an
%   accurate maping from clusters in CLUST to individual character symbols in
%   SYMS.  This function is meant to implement assigning clusters to characters
%   using the ``vp-ratio'' as defined by Ho and Nagy in their 2000 paper:
%   ``OCR with no shape training''
%   
%   CLUST should be a struct like that returned from cluster_comps()
%
%   COMPS should be a truct like that returned from get_comps()
%
%   SYMS should be a struct like that returned from create_alphabet()
%
%   MAP returned will be a vector of length CLUST.num, each entry of which lists
%   the index into SYMS that this cluster represents.
%
%   NOTE that the mappings are determined in order starting from the first 
%   cluster, and checking symbols in order from first to last, unless the local 
%   variable 'order' is overridden.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: word_lookup_map.m,v 1.1 2007-01-09 00:13:50 scottl Exp $
%
% REVISION HISTORY
% $Log: word_lookup_map.m,v $
% Revision 1.1  2007-01-09 00:13:50  scottl
% initial check-in.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%when should a ratio be considered good enough to make the assignment
acceptance_thresh = .75;

%you can define an order to explore the mappings in, or leave blank to
%use default ordering (1:Syms.num)
order = [];


% CODE START %
%%%%%%%%%%%%%%
if nargin < 3
    error('incorrect number of arguments');
elseif nargin > 3
    process_optional_args(varargin{:});
end

if ~isfield(Clust, 'model_spaces') || ~Clust.model_spaces
    error('word lookup requires knowledge of space characters');
end

space_idx = find(strcmp(Clust.truth_label, ' '));
if isempty(space_idx)
    error('unable to locate the cluster representing spaces');
end
space_sym_idx = find(strcmp(Syms.val, ' '));
if isempty(space_sym_idx)
    error('unable to locate the Symbol representing spaces');
end

if isempty(order)
    order = repmat(1:Syms.num, Clust.num, 1);
end

%determine the sequence of cluster words
seq = get_cluster_seq(Comps, unique(Comps.line));
sym_map = cell(Clust.num,1);
map = zeros(Clust.num,1);
[sym_map{:}] = deal('.');
sym_map{space_idx} = ' ';
map(space_idx) = space_sym_idx;
clust_words = cell(Clust.num,1);
word_list = cell(0);
word_num = 1;
for ii=1:length(seq)
    curr_word = [];
    curr_line = seq{ii};
    while ~isempty(curr_line)
        if curr_line(1) == space_idx
            %end of a word.  Update counts if not empty.
            if ~isempty(curr_word)
                word_list{word_num} = curr_word;
                curr_word = [];
                word_num = word_num + 1;
            end
        else
            %continue to grow the current word
            curr_word = [curr_word, curr_line(1)];
            clust_words{curr_line(1)} = [clust_words{curr_line(1)}, word_num];
        end
        curr_line = curr_line(2:end);
    end
end

%determine the length of each word in the lexicon (for faster matching)
lex_num_words = length(Syms.words);
lex_lengths = zeros(lex_num_words,1);
for ii=1:lex_num_words
    lex_length(ii) = length(Syms.words{ii});
end
lex_lists = cell(max(lex_length),1);
for ii=1:length(lex_lists)
    lex_lists{ii} = find(lex_length == ii);
end

%now loop through unmapped clusters, attempting to find those that end up in 
%scores above a pre-defined threshold, using the order to guide when certain 
%symbols are tried.
idx = find(strcmp(sym_map, '.'));
while ~isempty(idx)
    for ii=1:size(order,2)
        sym_map{idx(1)} = Syms.val{order(idx(1),ii)};
        sym_map{idx(1)} = regexprep(sym_map{idx(1)}, ...
                          '([\.\[\]\(\)\|\^\$\*\+\?\{\}])', '\\$1');
        if size(sym_map{idx(1)},2) > 1
            %this is to ensure that characters get glued-together in a single
            %column using char
            sym_map{idx(1)} = sym_map{idx(1)}';
        end
        score = calc_vp_score(sym_map, word_list(clust_words{idx(1)}), ...
                              Syms.words, lex_lists);
        if score >= acceptance_thresh
            %valid mapping!
            map(idx(1)) = order(idx(1),ii);
            idx = idx(2:end);
        end
    end
    %if we reach this point, no valid mapping has been found.
    %@@?
    fprintf('doh!\n');
    sym_map{idx(1)} = '.';
    idx = idx(2:end)'
end



% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this function determines the vp-ratio (proportion of the patterns given that
%match valid dictionary words using the partial mapping passed).
function score = calc_vp_score(map, pattern_list, word_list, word_lengths)
num_pat = length(pattern_list);
match_count = 0;
if num_pat == 0
    score = 0;
    return;
end

for ii=1:num_pat
    this_pat = char(map(pattern_list{ii}))';
    this_len = length(this_pat);
    res = regexp(word_list(word_lengths{this_len}), this_pat, 'once');
    if ~isempty(cell2mat(res))
        match_count = match_count + 1;
    end
end
score = match_count / num_pat;
