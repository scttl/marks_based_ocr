function map = word_lookup_map(Clust, Comps, Syms, varargin)
% WORD_LOOKUP_MAP  Determine best cluster mappings using word match ratios
%
%   MAP = word_lookup_map(CLUST, COMPS, SYMS, [var1,val1]...)
%   This function exploits word information taken from SYMS to determine an
%   accurate mapping from clusters in CLUST to individual character symbols in
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
% $Id: word_lookup_map.m,v 1.4 2007-01-19 19:53:43 scottl Exp $
%
% REVISION HISTORY
% $Log: word_lookup_map.m,v $
% Revision 1.4  2007-01-19 19:53:43  scottl
% ensure that punctuation symbols are compared correctly now that strmatch
% is used.  Use ISRI's expected reject symbol if no matching mapping
% can be found and all scores are 0.
%
% Revision 1.3  2007-01-18 19:15:58  scottl
% changed order to a cell array, updates to handle the case where there
% could be more than one space character.
%
% Revision 1.2  2007-01-13 18:21:36  scottl
% finished implementation.  Added modifications for capital letter
% handling.
%
% Revision 1.1  2007-01-09 00:13:50  scottl
% initial check-in.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%when should a ratio be considered good enough to make the assignment
acceptance_thresh = .75;

%you can define an order to explore the mappings in, or leave blank to
%use default ordering (1:Syms.num for each cluster, going in order from
%1:Clust.num)
order = {};


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

space_idx = find(strcmp(Clust.truth_label, ' '),1);
if isempty(space_idx)
    error('unable to locate the cluster representing spaces');
end
space_sym_idx = find(strcmp(Syms.val, ' '),1);
if isempty(space_sym_idx)
    error('unable to locate the Symbol representing spaces');
end

if isempty(order)
    order = cell(Clust.num,1);
    [order{:}] = deal(1:Syms.num);
end

%determine the sequence of cluster words
seq = get_cluster_seq(Comps, unique(Comps.line));
sym_map = cell(Clust.num,1);
map = zeros(Clust.num,1);
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
    if ~isempty(curr_word)
        word_list{word_num} = curr_word;
        curr_word = [];
        word_num = word_num + 1;
    end
end

%group lexicon words by length, and order by frequency to improve matching speed
lex_num_words = length(Syms.words);
[sort_idx,sort_idx] = sort(Syms.word_count, 'descend');
Syms.word_count = Syms.word_count(sort_idx);
Syms.words = Syms.words(sort_idx);
for ii=1:lex_num_words
    lex_length(ii) = length(Syms.words{ii});
end
lex_lists = cell(max(lex_length),1);
for ii=1:length(lex_lists)
    lex_lists{ii} = cell2mat(Syms.words(find(lex_length == ii)));
end

%now loop through unmapped clusters, attempting to find those that end up in 
%scores above a pre-defined threshold, using the order to guide when certain 
%symbols are tried.
idx = [1:space_idx-1,space_idx+1:Clust.num];
while ~isempty(idx)
    max_score = 0;
    max_idx = 1;
    max_sym = '~';  %reject character
    this_order = order{idx(1)};
    for ii=1:length(this_order)
        sym_map{idx(1)} = Syms.val{this_order(ii)};
        if size(sym_map{idx(1)},2) > 1
            %this is to ensure that characters get glued-together in a single
            %column using char
            sym_map{idx(1)} = sym_map{idx(1)}';
        end
        map(idx(1)) = this_order(ii);
        score = calc_vp_score(sym_map, map, word_list(clust_words{idx(1)}), ...
                              lex_lists);
        fprintf('.');
        if ii == 1 && score >= acceptance_thresh
            %valid mapping!
            map(idx(1)) = this_order(ii);
            fprintf('Score %f, sym: %s\n', score, Syms.val{map(idx(1))});
            idx = idx(2:end);
            break;
        elseif score > max_score
            max_score = score;
            max_idx = ii;
            max_sym = sym_map{idx(1)};
        end
    end
    if ii == length(this_order)
        %if we reach this point, no valid mapping has been found.
        map(idx(1)) = this_order(max_idx);
        fprintf('unable to find valid mapping for: %d\n', idx(1));
        fprintf('Using: Score %f, sym: %s\n', score, Syms.val{map(idx(1))});
        sym_map{idx(1)} = max_sym;
        idx = idx(2:end);
    end
end



% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this function determines the vp-ratio (proportion of the patterns given that
%match valid dictionary words using the partial mapping passed).
function score = calc_vp_score(sym_map, clust_map, pattern_list, word_list)
num_pat = length(pattern_list);
match_count = 0;
if num_pat == 0
    score = 0;
    return;
end

for ii=1:num_pat
    this_len = length(pattern_list{ii});
    idx = clust_map(pattern_list{ii}) ~= 0;
    if all(idx == 0)
        if ~isempty(word_list{this_len})
            match_count = match_count + 1;
        end
        continue;
    end
    str = char(sym_map(pattern_list{ii}(idx)))';
    if idx(1) == 1 && 65 <= str(1) && str(1) <= 90
        %since a lot of words will have capital letters at the start of the 
        %sentence, this will probably be missed by the lexicon, so convert 
        %this character to lower-case before attempting to match
        str(1) = str(1) + 32;
    end
    if this_len <= length(word_list) && ...
       any(strmatch(str, word_list{this_len}(:,idx), 'exact'))
        match_count = match_count + 1;
    elseif idx(end) == 0 && length(idx) > 1 && this_len > 1 && ...
           this_len-1 <= length(word_list)
        %since most of the lexicon will be missing words that end in punctuation
        %also try matching on words with the last cluster symbol missing.
        idx = idx(1:end-1);
        str = char(sym_map(pattern_list{ii}(idx)))';
        if any(strmatch(str, word_list{this_len-1}(:,idx), 'exact'))
            match_count = match_count + 1;
        end
    end
end
score = match_count / num_pat;
