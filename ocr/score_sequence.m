function score = score_sequence(seq, lg_start, lg_bg, varargin)
% SCORE_SEQUENCE  Calculate the overall log-score of the sequence passed
%
%   SCORE = SCORE_SEQUENCE(SEQ, LG_START_PROBS, LG_BIGRAM, [VAR1, VAL1]...)
%   SEQ should either be a vector containing a single sequence of integers
%   (representing indices of characters), or it should be a cell array
%   containing multiple such sequence vectors (1 vector per row).
%
%   LG_START_PROBS should be a vector listing the log probability of seeing the
%   character to which it indexes as the first character of a line.  It should
%   have length n, where n is the number of characters.
%
%   LG_BIGRAM should be a matrix listing the log probability of transitioning
%   from a character i at one point on a sequence, to a character j at the next
%   point.  The rows should represent the characters i, and the columns j.
%   There should be n rows and n columns, where n is the number of characters.
%
%   SCORE returned is a corresponding vector (or scalar if only 1 sequence
%   passed) giving the total score (log prob.) of seeing the entire sequence
%   based on bigram and start probabilities.  If any of the sequence is
%   missing, an upper-bound on the partial sequence is computed.
%
%   NOTE: one can get the perplexity (instead of log prob), by overridding the
%   calc_perplexity variable below
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: score_sequence.m,v 1.2 2006-11-07 02:53:37 scottl Exp $
%
% REVISION HISTORY
% $Log: score_sequence.m,v $
% Revision 1.2  2006-11-07 02:53:37  scottl
% calculate upper bounds for partial mappings.
%
% Revision 1.1  2006-10-29 17:27:37  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%set this to true to return perplexities instead of log probs.
calc_perplexity = false;


% CODE START %
%%%%%%%%%%%%%%
if nargin < 3
    error('incorrect number of arguments specified!');
elseif nargin > 3
    process_optional_args(varargin{:});
end

%some simple error checking
if isempty(seq)
    error('empty sequence passed');
elseif any(size(lg_bg) ~= length(lg_start))
    error('incorrect dimensions for bigram or start distributions');
end

if ~iscell(seq)
    %single row sequence passed.  Convert to a cell with one entry
    seq = {seq};
end

num_seq = length(seq);
score = zeros(num_seq,1);
max_char = length(lg_start);
for ii=1:length(seq)

    if isempty(seq{ii})
        warning('MBOCR:emptySeq', 'empty single sequence passed!\n');
        score(ii) = 0;
        continue;
    end
    unknown_idx = find(~(1 <= seq{ii} & seq{ii} <= max_char));

    %for the unknowns we calculate the best possible mapping based on its known
    %neighbours.
    rem_idx = [];
    for jj=1:length(unknown_idx)
        prev = unknown_idx(jj)-1;
        next = unknown_idx(jj)+1;
        if prev >= 1 && all(unknown_idx ~= prev) && ...
           next <= length(seq{ii}) && all(unknown_idx ~= next)
           %unknown sandwiched between 2 known values
           max_score = -Inf;
           max_idx = 0;
           for kk=1:max_char
               curr_score = lg_bg(seq{ii}(prev),kk) + lg_bg(kk,seq{ii}(next));
               if curr_score > max_score
                   max_score = curr_score;
                   max_idx = kk;
               end
           end
           seq{ii}(unknown_idx(jj)) = max_idx;
        elseif prev >= 1 && all(unknown_idx ~= prev)
            %transition *to* an unknown from a known
            [max_idx, max_idx] = max(lg_bg(seq{ii}(prev),:));
            seq{ii}(unknown_idx(jj)) = max_idx;
        elseif next <= length(seq{ii}) && all(unknown_idx ~= next)
            %transition *from* an unknown to a known
            [max_idx, max_idx] = max(lg_bg(:,seq{ii}(next)));
            seq{ii}(unknown_idx(jj)) = max_idx;
        else
            %consecutive blank transition, remove this index
            rem_idx = [rem_idx; unknown_idx(jj)];
        end
    end

    %since we assume log domain, probability of sequence is the *sum* of the
    %individual log probabilities
    keep_idx = setdiff(1:length(seq{ii}), rem_idx);
    if length(keep_idx) > 0 && keep_idx(1) == 1
        score(ii) = lg_start(seq{ii}(1));
        if calc_perplexity
            score(ii) = 1 /exp(score(ii));
        end
    end

    if length(keep_idx) > 1
        from_idcs = keep_idx;
        if from_idcs(end) == length(seq{ii})
            from_idcs = from_idcs(1:end-1);
        end
        to_idcs = keep_idx;
        if to_idcs(1) == 1
            to_idcs = to_idcs(2:end);
        end
        if ~isempty(rem_idx)
            %prune transition over gaps caused by the remaining unknowns
            from_idcs = setdiff(from_idcs, rem_idx-1);
            to_idcs = setdiff(to_idcs, rem_idx+1);
        end
        %also prune transitions between consecutive kept original unknowns
        %consider 1 -> a -> b -> 4  (a and b are originally unknown)
        %we estimate a as the best possible score transitioning from 1, and
        %b is estimated as the best possible score transition to 4.  But the
        %overall best transition may not include the optimal values a and b
        %when we factor in the cost of the transition between them.
        curr_pos = 1;
        while curr_pos <= length(from_idcs)
            if any(unknown_idx == from_idcs(curr_pos)) && ...
               any(unknown_idx == to_idcs(curr_pos))
               from_idcs = from_idcs([1:curr_pos-1, curr_pos+1:end]);
               to_idcs = to_idcs([1:curr_pos-1, curr_pos+1:end]);
           else
               curr_pos = curr_pos + 1;
           end
        end
        from_idcs = seq{ii}(from_idcs);
        to_idcs = seq{ii}(to_idcs);

        idcs = sub2ind(size(lg_bg), from_idcs, to_idcs);
        if calc_perplexity
            score(ii) = (score(ii) * prod(1 ./ exp(lg_bg(idcs))))^...
                        (1/length(keep_idx));
        else
            score(ii) = score(ii) + sum(lg_bg(idcs));
        end
    end
end


if calc_perplexity
    tot_score = mean(score);
else
    tot_score = sum(score);
end


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
