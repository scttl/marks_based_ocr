function [log_prob, ppl] = srilm_lm_score(pipe, fid, seq, varargin)
% SRILM_LM_SCORE Score seqeunce passed using the pipe handle
%
%   [LOG_PROB, PERPLEXITY] = SRILM_LM_SCORE(PIPE, OUT_FID, SEQ, [VAR1,VAL1]...)
%
%   PIPE should be a numeric handle denoting the descriptor of the pipe that we
%   will write the sequence passed to.  This pipe is setup (and returned) in a
%   call to srilm_lm_open
%
%   OUT_FID should be the file descriptor of the output file that we will read
%   sequence scores from.  This file descriptor is created via a call to
%   srilm_lm_open.
%
%   SEQ should be a long character array (newlines delimiting each line) that
%   will be scored via SRI-LM's language model
%
%   LOG_PROB returned will be a scalar listing the overall log 
%   probability of the sequence passed in SEQ
%
%   PERPLEXITY returned will be a vector (or scalar) listing the total 
%   perpplexity of each sequence passed in FILES
%
%   NOTE: We require that SRILM is installed, and the tool ngram is in the
%   users path
%
%   NOTE: any default values defined below in the LOCAL VARS can be overridden 
%   by passing the name of parameter, and its new value on the command line.
% 

% CVS INFO %
%%%%%%%%%%%%
% $Id: srilm_lm_score.m,v 1.2 2006-11-22 17:11:13 scottl Exp $
%
% REVISION HISTORY
% $Log: srilm_lm_score.m,v $
% Revision 1.2  2006-11-22 17:11:13  scottl
% rewritten to use popen for efficient scoring.
%
% Revision 1.1  2006-11-13 17:58:23  scottl
% initial check-in.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%sequence of characters used to denote the end of a sequence of lines (and that
%log probability and perplexity score should be calculated
%NOTE: this value *must* match the escape_seq in srilm_lm_open
escape_seq = 'EOS';


% CODE START %
%%%%%%%%%%%%%%

if nargin < 3
    error('incorrect number of arguments');
elseif nargin > 3
    process_optional_args(varargin{:});
end

%some sanity checking on the arguments passed
if ~ischar(seq)
    error('sequence passed must be a character array');
end

log_prob = NaN;
ppl = NaN;

if ~strcmp(seq(end), char(10))
    seq = [seq, char(10)];
end
num_bytes = popenw(pipe, [seq, escape_seq, char(10)], 'char');
line = fgetl(fid);
ii = 0;
while ii < 10 && ((ischar(line) && isempty(regexp(line, 'logprob'))) || ...
      (isnumeric(line) && line == -1))
    line = fgetl(fid);
    ii = ii+1;
end
if ii < 10
    %if run correctly, line should be something like the following:
    %0 zeroprobs, logprob= -1865.95 ppl= 23.497 ppl1= 24.987
    tmp_res = regexp(line, 'logprob= (\S+) ppl= (\S+)', 'tokens');
    log_prob = str2num(tmp_res{1}{1});
    ppl = str2num(tmp_res{1}{2});
else
    fprintf('unable to read valid line\n');
end
