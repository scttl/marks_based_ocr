function [map, valid_acc] = simul_word_lookup_map(Clust, Comps, Syms, varargin)
% SIMUL_WORD_LOOKUP_MAP Determine best cluster mappings using constrained lookup
%
%   [MAP, VALID_ACC] = simul_word_lookup_map(CLUST, COMPS, SYMS, [var1,val1]...)
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


% CVS INFO %
%%%%%%%%%%%%
% $Id: simul_word_lookup_map.m,v 1.1 2007-05-08 01:00:46 scottl Exp $
%
% REVISION HISTORY
% $Log: simul_word_lookup_map.m,v $
% Revision 1.1  2007-05-08 01:00:46  scottl
% initial check-in
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%when should a ratio be considered good enough to make the assignment
%@@to implement!
acceptance_thresh = .75;

%should we add counts to the first position of mappings by replacing the first
%lowercase letter with its uppercase equivalent?
add_first_pos_up_let = false;

%should we try extending words whose last symbol isn't punctuation by adding
%counts using each valid punctuation symbol?
valid_punct_syms = '.,:;!?';
add_last_pos_punct_syms = false;


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

valid_acc = NaN;
sym_map = cell(Clust.num,1);
map = zeros(Clust.num,1);
sym_map{space_idx} = ' ';
map(space_idx) = space_sym_idx;

%group lexicon words by length, and order by frequency to improve matching speed
lex_num_words = length(Syms.words);
[sort_idx,sort_idx] = sort(Syms.word_count, 'descend');
Syms.word_count = Syms.word_count(sort_idx);
Syms.words = Syms.words(sort_idx);
for ii=1:lex_num_words
    lex_length(ii) = length(Syms.words{ii});
end
lex_max_wordlen = max(lex_length);
lex_lists = cell(lex_max_wordlen,1);
lex_num_lists = cell(lex_max_wordlen,1);
lex_wc = cell(lex_max_wordlen,1);
for ii=1:length(lex_lists)
    idx = find(lex_length == ii);
    lex_lists{ii} = cell2mat(Syms.words(idx));
    lex_wc{ii} = Syms.word_count(idx);
    if add_first_pos_up_let
        %add uppercase symbols in the appropriate places
        vals = double(lex_lists{ii}(:,1));
        low_idx = find(vals >= 97 & vals <= 122);
        if ~isempty(low_idx)
            new_items = lex_lists{ii}(low_idx,:);
            new_items(:,1) = char(vals(low_idx) - 32);
            new_wc = lex_wc{ii}(low_idx);
            lex_lists{ii} = [lex_lists{ii}; new_items];
            lex_wc{ii} = [lex_wc{ii}; new_wc];
        end
    end
    if ii > 1 && add_last_pos_punct_syms && ~isempty(valid_punct_syms)
        %add punctation marks in the appropriate positions, using shorter words
        vals = double(lex_lists{ii-1}(:,end));
        np_idx = find(~ismember(vals, double(valid_punct_syms)));
        if ~isempty(np_idx)
            new_items = lex_lists{ii-1}(np_idx,:);
            new_wc = lex_wc{ii-1}(np_idx);
            for jj=1:length(valid_punct_syms);
                this_new = new_items;
                this_new(:,end) = valid_punct_syms(jj);
                lex_lists{ii} = [lex_lists{ii}; [new_items, ...
                           repmat(valid_punct_syms(jj),size(new_items,1),1)]];
                lex_wc{ii} = [lex_wc{ii}; new_wc];
            end
        end
    end
    lex_num_lists{ii} = numerize(lex_lists{ii});
end

%determine the sequence of unique cluster words, numerize them, group them by
%length, then find the associated list of matching lexicon words
[word_list, clust_words] = get_cluster_word_seq(Comps, space_idx);
clust_num_words = length(word_list);
for ii=1:clust_num_words
    cw_length(ii) = length(word_list{ii});
end
cw_max_wordlen = max(cw_length);
cw_uniq_lists = cell(cw_max_wordlen,1);
cw_uniq_num_lists = cell(cw_max_wordlen,1);
cw_uniq_idx = zeros(clust_num_words,1);
lex_matches = cell(cw_max_wordlen,1);
lex_match_length = zeros(clust_num_words,1);
cw_uniq_length = zeros(clust_num_words,1);
cw_uniq_valid_idx = true(clust_num_words,1);
for ii=1:length(cw_uniq_lists)
    idx = find(cw_length == ii);
    cw_uniq_lists{ii} = cell2mat(word_list(idx));
    [cw_uniq_lists{ii},uniq_idx,uniq_idx] = unique(cw_uniq_lists{ii}, 'rows');
    cw_uniq_idx(idx) = uniq_idx;
    cw_uniq_num_lists{ii} = numerize(cw_uniq_lists{ii});
    this_length = size(cw_uniq_lists{ii},1);
    cw_uniq_length(ii) = this_length;
    lex_matches{ii} = cell(this_length,1);
    if ii <= length(lex_lists)
        for jj=1:this_length
            %loosen the lexicon match requirements so that the same value can
            %be matched for distinct corpus string values.  e.g. if corpus
            %string is 1 2 3 2, then this should also match lexicon words of the
            %form 1 1 1 1, 1 1 2 1, 1 2 1 2, and 1 2 2 2 (i.e any string in
            %which the second and fourth symbols are the same
            lex_idx = 1:size(lex_num_lists{ii},1);
            for val=unique(cw_uniq_num_lists{ii}(jj,:));
                pos=find(cw_uniq_num_lists{ii}(jj,:) == val);
                if length(pos) > 1
                    comp_val = lex_num_lists{ii}(lex_idx,pos(1));
                    lex_i2 = find(all(lex_num_lists{ii}(lex_idx,pos(2:end)) ...
                             ==  repmat(comp_val, 1, length(pos)-1),2));
                    lex_idx = lex_idx(lex_i2);
                end
            end
            lex_matches{ii}{jj} = lex_lists{ii}(lex_idx,:);
            if isempty(lex_idx)
                cw_uniq_valid_idx(idx(uniq_idx == jj)) = false;
            else
                lex_match_length(idx(uniq_idx == jj)) = ...
                                 size(lex_matches{ii}{jj},1);
            end
        end
    else
        cw_uniq_valid_idx(idx) = false;
    end
end

%now loop through unmapped clusters, attempting to constrain matches until a
%single valid mapping remains
unmapped_idx = [1:space_idx-1,space_idx+1:Clust.num];
while ~isempty(unmapped_idx)
    %find the cluster word containing the first unmapped identifier that is
    %valid and has the fewest matching entries
    val_map = false(Clust.num, Syms.num);
    for ii=1:length(map)
        if map(ii) == 0
            val_map(ii,[1:space_sym_idx-1,space_sym_idx+1:end]) = true;
        else
            val_map(ii,map(ii)) = true;
        end
    end
    seq_idx = clust_words{unmapped_idx(1)};
    seq_idx = seq_idx(cw_uniq_valid_idx(seq_idx));
    lens = lex_match_length(seq_idx);

    while(true)
        %examine the matching lexicon words to ensure constraints are still
        %satisfied.
        if isempty(lens)
            %no more valid words containing this cluster, arbitrarily choose
            %the first symbol that is still valid
            sym_idx = find(val_map(unmapped_idx(1),:),1);
            map(unmapped_idx(1)) = sym_idx;
            sym_map{unmapped_idx(1)} = Syms.val{sym_idx};
            fprintf('no valid left, mapping %d to %s\n', unmapped_idx(1), ...
                    Syms.val{sym_idx});
            unmapped_idx = unmapped_idx(2:end);
            break;
        end
        [min_len,min_len_idx] = min(lens);
        min_len_idx = seq_idx(min_len_idx);
        if min_len == 0
            %each of these words is invalid
            cw_uniq_valid_idx(min_len_idx) = false;
            [seq_idx,pos] = setdiff(seq_idx, min_len_idx);
            lens = lens(pos);
            continue;
        end
        min_len_idx = min_len_idx(1);
        wl = cw_length(min_len_idx);
        off = cw_uniq_idx(min_len_idx);
        [cl,pos] = unique(cw_uniq_lists{wl}(off,:));
        this_lex = lex_matches{wl}{off};
        lex_idx = 1:size(this_lex,1);
        %reduce this list of words to those that contain matches to already
        %mapped cluster values
        map_idx = find(map(cl) ~= 0);
        while ~isempty(lex_idx) && ~isempty(map_idx)
            map_pos = pos(map_idx(1));
            map_val = sym_map{cl(map_idx(1))};
            new_idx = find(this_lex(lex_idx,map_pos) == map_val);
            lex_idx = lex_idx(new_idx);
            map_idx = map_idx(2:end);
        end
        if isempty(lex_idx)
            %mark this word as invalid
            cw_uniq_valid_idx(min_len_idx) = false;
            [seq_idx,pos] = setdiff(seq_idx, min_len_idx);
            lens = lens(pos);
            continue;
        else
            this_lex = this_lex(lex_idx,:);
        end
        %check if there are uniquely satisfied lexicon matches for the
        %remaining unmapped clusters
        u_idx = find(map(cl) == 0);
        u_cl = cl(u_idx);
        u_pos = pos(u_idx);
        u_map = false(length(u_cl),Syms.num);
        for ii=1:length(u_cl)
            syms = unique(this_lex(:,u_pos(ii)));
            for jj=1:length(syms)
                o = find(strcmp(Syms.val, syms(jj)));
                if val_map(u_cl(ii),o)
                    u_map(ii,o) = true;
                end
            end
        end
        if any(sum(u_map,2) == 1)
            %at least one identifier has a unique mapping
            for ii=1:size(u_map,1)
                if sum(u_map(ii,:)) == 1
                    map(u_cl(ii)) = find(u_map(ii,:));
                    sym_map{u_cl(ii)} = Syms.val{map(u_cl(ii))};
                    val_map(u_cl(ii),:) = u_map(ii,:);
                    unmapped_idx = setdiff(unmapped_idx, u_cl(ii));
                    fprintf('mapping %d to symbol %s\n', u_cl(ii), ...
                            sym_map{u_cl(ii)});
                elseif sum(u_map(ii,:)) > 1
                    val_map(u_cl(ii),:) = val_map(u_cl(ii),:) & u_map(ii,:);
                %elseif no valid lexicon matchings for a cluster, don't refine
                %the set of valid symbols
                end
            end
            break;
        elseif all(sum(u_map,2) ~= 0)
            %no invalid and no unique mappings found.  Restrict the valid
            %mapping choices for the unmapped clusters and set this word to
            %invalid (since examined)
            for ii=1:size(u_map,1)
                val_map(u_cl(ii),:) = val_map(u_cl(ii),:) & u_map(ii,:);
            end
        end
        %mark this now examined word as invalid
        cw_uniq_valid_idx(min_len_idx) = false;
        [seq_idx,pos] = setdiff(seq_idx,min_len_idx);
        lens = lens(pos);
    end
end
