function [map, valid_acc] = word_lookup_map(Clust, Comps, Syms, varargin)
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
%   VALID_ACC is only returned if truth label information is available
%   for the clusters and calc_valid_acc is set below.  It allows us to
%   define an alternate measure of accuracy, dictating what portion of the time
%   clusters are correctly identified when restricting ourselves to those that 
%   either appear in words in the dictionary, or that are made up entirely of 
%   clusters that also appear in the dictionary
%
%   NOTE that the mappings are determined in order starting from the first 
%   cluster, and checking symbols in order from first to last, unless the local 
%   variable 'order' is overridden.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: word_lookup_map.m,v 1.9 2007-04-10 19:38:41 scottl Exp $
%
% REVISION HISTORY
% $Log: word_lookup_map.m,v $
% Revision 1.9  2007-04-10 19:38:41  scottl
% removed calculation of cluster word sequence to separate file
%
% Revision 1.8  2007-04-10 15:49:53  scottl
% implemented ability to use simple shape information as a tie breaking procedure
%
% Revision 1.7  2007-02-02 05:51:28  scottl
% added density based tie breaking.
%
% Revision 1.6  2007-02-01 19:08:57  scottl
% added ability to restrict ordering to those matching a particular offset
% class.  Implemented alternate accuracy calculation based on words
% either being in the dictionary, or composed entirely of clusters that belong
% to dictionary words.
%
% Revision 1.5  2007-01-26 23:25:40  scottl
% rewrote vp ratio calculations to constrain wildcards (if they refer
% to the same cluster).  This currently makes things slower, but perhaps
% profiling can find a few places to tweak things.
%
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

%setting this to true restricts the ordering to just those symbols estimated to
%be in the same class as the cluster currently being mapped
restrict_order_to_class = false;

%setting this to true enables us to return which clusters are valid as well as
%which were predicted correctly (but requires ground truth information)
calc_valid_acc = false;

%should we use density information to break ties?
break_ties_via_density = false;

%should we break ties via shape information from the symbols?
break_ties_via_shape = false;
%if shape information used to break ties, what method should be used to
%normalize symbol image
resize_method = 'nearest';

%should we use the truth label to break ties?  Only works if we have the truth 
%label information, and it happens to be amongst the tied mappings.
break_ties_via_truth = false;


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
valid_acc = NaN;
sym_map = cell(Clust.num,1);
map = zeros(Clust.num,1);
sym_map{space_idx} = ' ';
map(space_idx) = space_sym_idx;
[word_list, clust_words] = get_cluster_word_seq(Comps, space_idx);

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
    if restrict_order_to_class
        this_order = this_order(Syms.class(this_order) == Clust.class(idx(1)));
    end
    for ii=1:length(this_order)
        sym_map{idx(1)} = Syms.val{this_order(ii)};
        if size(sym_map{idx(1)},2) > 1
            %this is to ensure that characters get glued-together in a single
            %column using char
            sym_map{idx(1)} = sym_map{idx(1)}';
        end
        map(idx(1)) = this_order(ii);
        score = calc_vp_score(sym_map, map, word_list(clust_words{idx(1)}), ...
                              lex_lists,false);
        fprintf('.');
        if ii == 1 && score >= acceptance_thresh
            %valid mapping!
            map(idx(1)) = this_order(ii);
            fprintf('Score %f, sym: %s', score, Syms.val{map(idx(1))});
            if ~isempty(Clust.truth_label(idx(1)))
                fprintf(' true val: %s', Clust.truth_label{idx(1)});
            end
            fprintf('\n');
            idx = idx(2:end);
            break;
        elseif score > max_score
            max_score = score;
            max_idx = ii;
            max_sym = sym_map{idx(1)};
        elseif (break_ties_via_density || break_ties_via_shape) && ...
               score == max_score
            %tie!
            max_idx = [max_idx, ii];
            max_sym = [max_sym, sym_map(idx(1))];
        elseif break_ties_via_truth && score == max_score
            true_val = Clust.truth_label(idx(1));
            if ~isempty(true_val) && strcmp(true_val, sym_map{idx(1)})
                max_idx = ii;
                max_sym = sym_map{idx(1)};
            end
        end
    end
    if ii == length(this_order)
        %if we reach this point, no valid mapping has been found.
        if length(max_idx) > 1
            %a tie has occured
            fprintf('tie between: '); fprintf('%s, ',max_sym{:}); fprintf('\n');
            sym_idx = this_order(max_idx);
            if break_ties_via_density
                [best_idx,best_idx] = min(abs(Syms.density(sym_idx) - ...
                                      Clust.density(idx(1))));
            else
                %use shape to break ties, after normalizing candidates aspect
                %ratios
                rsz_syms = cell(length(sym_idx),1);
                clust_height = size(Clust.avg{idx(1)},1);
                for jj=1:length(sym_idx)
                    if sum(sum(Syms.img{sym_idx(jj)})) == 0
                        %trying to resize a blank image
                        rsz_syms{jj} = Syms.img{sym_idx(jj)};
                    else
                        scale = clust_height / size(Syms.img{sym_idx(jj)},1);
                        rsz_syms{jj} = imresize(Syms.img{sym_idx(jj)}, ...
                                                scale, resize_method);
                    end
                end
                [best_idx,best_idx] = min(hausdorff_dist(Clust.avg{idx(1)}, ...
                                          rsz_syms));
            end
            max_idx = max_idx(best_idx);
            max_sym = max_sym{best_idx};
        end
        map(idx(1)) = this_order(max_idx);
        fprintf('unable to find valid mapping for: %d\n', idx(1));
        fprintf('Using: Score %f, sym: %s', max_score, Syms.val{map(idx(1))});
        if ~isempty(Clust.truth_label(idx(1)))
            fprintf(' true val: %s', Clust.truth_label{idx(1)});
        end
        fprintf('\n');
        sym_map{idx(1)} = max_sym;
        idx = idx(2:end);
    end
end

if calc_valid_acc
    valid_clust = logical(zeros(Clust.num,1));
    [dummy,valid_words] = calc_vp_score(Clust.truth_label, 1:Clust.num, ...
                           word_list, lex_lists, false);
    for ii=find(valid_words)'
        valid_clust(word_list{ii}) = true;
    end
    %can also be valid if not in dictionary but made entirely of clusters that
    %are valid (in dictionary)
    for ii=find(~valid_words)'
        if all(valid_clust(word_list{ii}))
            valid_words(ii) = true;
        end
    end

    correct = zeros(Clust.num,1);
    seen = zeros(Clust.num,1);
    for ii=find(valid_words)'
        seq = word_list{ii};
        seen(seq) = seen(seq) + 1;
        actual = char(Clust.truth_label(seq));
        pred = char(sym_map(seq));
        match = actual == pred;
        correct(seq(match)) = correct(seq(match)) + 1;
    end
    valid_acc = [seen, correct, (correct ./ seen)];
end


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this function determines the vp-ratio (proportion of the patterns given that
%match valid dictionary words using the partial mapping passed).
function [score, matches] = calc_vp_score(sym_map, clust_map, pattern_list, ...
                 word_list, list_matches)

num_pat = length(pattern_list);
matches = logical(zeros(num_pat,1));
if num_pat == 0
    score = 0;
    return;
end

for ii=1:num_pat
    this_len = length(pattern_list{ii});
    word_len = length(word_list);
    if this_len > word_len + 1
        continue;
    else
        idx = clust_map(pattern_list{ii}) ~= 0;
        this_pat = pattern_list{ii};
        str = char(sym_map(this_pat(idx)))';
        found_match = false;
        check_lengths = [];
        if this_len <= word_len
            check_lengths = 0;
        end
        %since most of the lexicon will be missing words that end in 
        %punctuation also try matching on words with the last cluster 
        %symbol missing.
        if this_len <= word_len-1 && idx(end) == 0 && length(idx) > 1 && ...
           this_len > 1
            check_lengths = [check_lengths, 1];
        end
        while ~found_match && ~isempty(check_lengths)
            twl = word_list{this_len - check_lengths(1)};
            idx = idx(1:end - check_lengths(1));
            this_pat = this_pat(1:end - check_lengths(1));
            match_rows = 1:size(twl,1);
            wild = find(idx == 0);
            while length(wild) >= 2
                %look for constraints amongst these columns
                match_cols = this_pat(wild(1)) == this_pat(wild);
                mc = wild(match_cols);
                while length(mc) > 1
                    mr = twl(match_rows,mc(1)) == twl(match_rows,mc(end));
                    match_rows = match_rows(mr);
                    mc = mc(1:end-1);
                end
                wild = wild(~match_cols);
            end
            
            if all(idx == 0)
                if ~isempty(twl(match_rows,:))
                    found_match = true;
                    continue;
                end
            end
            if any(strmatch(str, twl(match_rows,idx), 'exact'))
                found_match = true;
            elseif idx(1) == 1 && 65 <= str(1) && str(1) <= 90
                %since a lot of words will have capital letters at the start of 
                %the sentence, this will probably be missed by the lexicon, so 
                %convert this character to lower-case and attempt to match
                tmp_str = str;
                tmp_str(1) = tmp_str(1) + 32;
                if any(strmatch(tmp_str, twl(match_rows,idx), 'exact'))
                    found_match = true;
                end
            end
            check_lengths = check_lengths(2:end);
        end
        if found_match
            matches(ii) = true;
            if list_matches
                fprintf('%s\n', str);
            end
        end
    end
end
score = sum(matches) / num_pat;
