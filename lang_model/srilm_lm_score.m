function [log_prob, ppl] = srilm_lm_score(lm, Files, varargin)
% SRILM_LM_SCORE Score seqeunces passed using the language model passed.
%
%   [LOG_PROB, PERPLEXITY] = SRILM_LM_SCORE(MODEL, FILES, [VAR1, VAL1]...)
%
%   MODEL should be a string denoting where to find a SRILM language model file
%   (path and name).  The file should be in the ARPA format SRILM expects.  This
%   is returned in srilm_lm_init().
%
%   FILES should either be a string or a cell array listing the path and name of
%   1 or more plain-text sequence files.  Since SRILM operates on word tokens, 
%   we assume that if characters are to be used as tokens, that the file has 
%   already been converted to have each character space-separated (and spaces 
%   mapped to another symbol)
%
%   LOG_PROB returned is will be a vector (or scalar) listing the overall log 
%   probability of each sequence passed in FILES
%
%   PERPLEXITY returned will be a vector (or scalar) listing the total 
%   perpplexity of each sequence passed in FILES
%
%   NOTE: any default values defined below in the LOCAL VARS can be overridden 
%   by passing the name of parameter, and its new value on the command line.
% 
%   NOTE: We require that SRILM is installed, and the tool ngram is in the
%   users path
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: srilm_lm_score.m,v 1.1 2006-11-13 17:58:23 scottl Exp $
%
% REVISION HISTORY
% $Log: srilm_lm_score.m,v $
% Revision 1.1  2006-11-13 17:58:23  scottl
% rewritten to use popen for efficient scoring calculations.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%what order n-grams should we use?
order = '3';

%use a specific input vocabulary?  Note that if doing so, the space character 
%should be mapped to some other character: like '_'
input_vocab = '';

%should we skip unknowns in our score calculation?  By default, we backoff to 
%lower models and try to estimate them
skip_oovs = false;

%any other parameters to pass to ngram?  see man ngram for details
other_params = '';



% CODE START %
%%%%%%%%%%%%%%

if nargin < 1
    error('incorrect number of arguments');
elseif nargin > 1
    process_optional_args(varargin{:});
end

%some sanity checking on the argument passed
if ~ (ischar(Files) || iscell(Files))
    error('Files passed are in the wrong format');
end

if ischar(Files)
    Files = {Files};
end

log_prob = zeros(length(Files),1);
ppl = zeros(length(Files),1);
for ii=1:length(Files)
    cmd = ['ngram -order ', order, ' -lm ', lm, ' -ppl ', Files{ii}, ...
           ' ', other_params, ' -debug 0'];
    if ~isempty(input_vocab)
        cmd = [cmd, '-vocab ', input_vocab];
    end
    if skip_oovs
        cmd = [cmd, '-skipoovs'];
    end
    [s,w] = unix(cmd);
    if s ~= 0
        error('error running ngram on file %s:\n%s', Files{ii}, cmd);
    end

    %if run correctly, output should be something like the following:
    %file file.srilm: 26 sentences, 1335 words, 0 OOVs
    %0 zeroprobs, logprob= -1865.95 ppl= 23.497 ppl1= 24.987
    tmp_res = regexp(w, 'logprob= (\S+) ', 'tokens');
    log_prob(ii) = str2num(tmp_res{1}{1});
    tmp_res = regexp(w, 'ppl= (\S+) ', 'tokens');
    ppl(ii) = str2num(tmp_res{1}{1});
end
