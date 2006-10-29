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
%   based on bigram and start probabilities.
%
%   NOTE: once can get the perplexity (instead of log prob), by overridding the
%   calc_perplexity variable below
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: score_sequence.m,v 1.1 2006-10-29 17:27:37 scottl Exp $
%
% REVISION HISTORY
% $Log: score_sequence.m,v $
% Revision 1.1  2006-10-29 17:27:37  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%set this to true to return perplexities instead of log probs.
calc_perplexity = false;


% CODE START %
%%%%%%%%%%%%%%
tic;
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
score = -Inf(num_seq,1);
for ii=1:length(seq)

    if isempty(seq{ii})
        warning('MBOCR:emptySeq', 'empty single sequence passed!\n');
        score(ii) = 0;
        continue;
    end

    %since we assume log domain, probability of sequence is the *sum* of the
    %individual log probabilities
    score(ii) = lg_start(seq{ii}(1));
    if calc_perplexity
        score(ii) = 1 /exp(score(ii));
    end

    if length(seq{ii}) > 1
        from_idcs = seq{ii}(1:end-1);
        to_idcs = seq{ii}(2:end);
        idcs = sub2ind(size(lg_bg), from_idcs, to_idcs);
        if calc_perplexity
            score(ii) = (score(ii) * prod(1 ./ exp(lg_bg(idcs))))^...
                        (1/length(seq{ii}));
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

fprintf('%.2fs: score for all sequences is %.2f\n',toc,tot_score);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
