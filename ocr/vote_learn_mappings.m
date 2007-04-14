function [order,score] = vote_learn_mappings(Clust, Comps, Syms, varargin)
% VOTE_LEARN_MAPPINGS  Determine most likely cluster to symbol mappings
%
%   [order, score] = VOTE_LEARN_MAPPINGS(Clust, Comps, Syms, [var1,val1]...)
%   This function attempts to learn a good mapping from each cluster to a
%   single character symbol by a weighted vote scheme based on corpus words
%   that match the format of each cluster word.
%   
%   Clust should be a struct like that returned from cluster_comps()
%
%   Comps should be a struct like that returned from get_comps()
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
% $Id: vote_learn_mappings.m,v 1.1 2007-04-14 20:16:29 scottl Exp $
%
% REVISION HISTORY
% $Log: vote_learn_mappings.m,v $
% Revision 1.1  2007-04-14 20:16:29  scottl
% initial check-in
%


% GLOBAL VARS %
%%%%%%%%%%%%%%%


% LOCAL VARS %
%%%%%%%%%%%%%%

%should we limit matching words, to those that take into account the current
%partial mapping?
limit_to_map = true;

%what percentage of the counts should be counted as additional votes?  Set this
%to 0, and each matching word gets a single vote (regardless of the number of
%times it occurs in the document).  Set this to 1 and each occurence of each
%word is counted as a vote.
word_count_weight_pct = 0;

%should we add counts to the first position of mappings by replacing the first
%lowercase letter with its uppercase equivalent?
add_first_pos_up_let = false;

%should we try extending words whose last symbol isn't punctuation by adding
%counts using each valid punctuation symbol?
valid_punct_syms = '.,:;!?';
add_last_pos_punct_syms = false;


% CODE START %
%%%%%%%%%%%%%%
tic;
if nargin < 2
    error('incorrect number of arguments specified!');
elseif nargin > 2
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

%numerize the lexicon words, creating a separate matrix for each length
lex_length = zeros(length(Syms.words),1);
for ii=1:length(Syms.words)
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
    lex_num_lists{ii} = numerize(lex_lists{ii});
end

%get the list of cluster words, and numerize them
[word_list, word_idx] = get_cluster_word_seq(Comps, space_idx);
cl_num_list = numerize(word_list);

%create an inverse map from character values to their corresponding symbol
%entry
inv_sym_map = [];
for ii=1:Syms.num
    if length(Syms.val{ii}) == 1
        inv_sym_map(double(Syms.val{ii}),1) = ii;
    end
end

%go through the clusters, counting votes for each valid mapping word
order = cell(Clust.num, 1);
score = zeros(Clust.num, Syms.num);
score(space_idx, space_sym_idx) = Inf;
order{space_idx} = [space_sym_idx,1:space_sym_idx-1,space_sym_idx+1:Syms.num];
for ii=[1:space_idx-1,space_idx+1:Clust.num]
    for jj=1:length(word_idx{ii})
        idx = word_idx{ii}(jj);
        w = word_list{idx};
        len = length(w);
        pos = find(w == ii);
        if len <= length(lex_num_lists)
            lex_idx = ismember(lex_num_lists{len}, cl_num_list{idx}, 'rows');
            if limit_to_map
                lex_idx = find(lex_idx);
                map_pos = find(w < ii);
                while ~isempty(map_pos)
                    lex_idx =lex_idx(find(lex_lists{len}(lex_idx,map_pos(1)) ...
                             == Syms.val{order{w(map_pos(1))}(1)}));
                    map_pos = map_pos(2:end);
                end
            end
            %@@initially counts were multiplied
            %weights = word_count_weight_pct*lex_wc{len}(lex_idx);

            num_match = length(lex_idx);
            counts = lex_wc{len}(lex_idx);
            sym_idx = inv_sym_map(double(lex_lists{len}(lex_idx, pos(1))));
            if add_first_pos_up_let && pos(1) == 1
                low_idx = find(lex_lists{len}(lex_idx,pos(1)) <= 122 & ...
                               lex_lists{len}(lex_idx,pos(1)) >= 97);
                num_match = num_match + length(low_idx);
                counts = [counts; counts(low_idx)];
                sym_idx = [sym_idx; ...
                           inv_sym_map(double(...
                           lex_lists{len}(lex_idx(low_idx),pos(1))) - 32)];
            end
            if add_last_pos_punct_syms && pos(end) == len && len > 1
                new_idx = ismember(lex_num_lists{len-1}, ...
                          cl_num_list{idx}(1:end-1), 'rows');
                if limit_to_map
                    new_idx = find(new_idx);
                    map_pos = find(w(1:end-1) < ii);
                    while ~isempty(map_pos)
                        new_idx = new_idx(find(lex_lists{len-1}(new_idx, ...
                                  map_pos(1)) == ...
                                  Syms.val{order{w(map_pos(1))}(1)}));
                        map_pos = map_pos(2:end);
                    end
                end
                punct_idx = find(lex_lists{len-1}(new_idx,end) >= 65);
                num_match = num_match + ...
                            length(punct_idx)*length(valid_punct_syms);
                counts = [counts; repmat(lex_wc{len-1}(new_idx(punct_idx)), ...
                          length(valid_punct_syms),1)];
                sym_idx = [sym_idx; reshape(repmat(inv_sym_map(double(...
                           valid_punct_syms)), ...
                           length(punct_idx), 1), ...
                           length(punct_idx)*length(valid_punct_syms), 1)];
            end
            weights = word_count_weight_pct .* (counts ./ sum(counts));

            for kk=1:length(sym_idx)
                score(ii,sym_idx(kk)) = score(ii,sym_idx(kk)) + weights(kk);
            end
        end
    end
    [score(ii,:), order{ii}] = sort(score(ii,:),'descend');
    fprintf('%.2fs: cluster %d --> %s, score of %d\n', toc, ii, ...
             Syms.val{order{ii}(1)}, score(ii,1));
end

%since we want score to represent a distance, invert the sorted scores before
%converting to a cell array
score(score == 0) = eps;
score = mat2cell(1 ./ score, ones(Clust.num,1));

fprintf('%.2fs: ordering complete\n', toc);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%return a numerized representation of each string or numeric sequence passed.
%Each occurence of the same symbol is given the same numerical value, and the
%values are numbered starting at 1 (no two different symbols may receive the
%same value).  ex. 'sassy' gets 12113, and [4 1 10 1 9] gets 12324
%Inputs should list one entry per row (or cell item).  The output will be of
%the same type and size as the input (though strings are converted to double
%arrays), but with values changed to vectors of numerization values.  If the
%input is not a cell array, each column must be the same length
function out = numerize(in)

if iscell(in)
    out = cell(size(in));
    for ii=1:length(in)
        str = in{ii};
        nm = zeros(size(str));
        curr = 1;
        for jj=1:length(nm)
            pos = find(str(1:jj-1) == str(jj),1);
            if isempty(pos)
                nm(jj) = curr;
                curr = curr+1;
            else
                nm(jj) = nm(pos);
            end
        end
        out{ii} = nm;
    end
else
    out = zeros(size(in));
    for ii=1:size(in,1)
        curr = 1;
        for jj=1:size(in,2)
            pos = find(in(ii,1:jj-1) == in(ii,jj),1);
            if isempty(pos)
                out(ii,jj) = curr;
                curr = curr+1;
            else
                out(ii,jj) = out(ii,pos);
            end
        end
    end
end
