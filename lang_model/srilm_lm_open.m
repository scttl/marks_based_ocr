function [pipe, srilm_seq_file, srilm_fid] = srilm_lm_open(lm, varargin)
% SRILM_LM_OPEN Setup a pipe to SRI-LM to allow sequences to be scored
%
%   [PIPE, OUT_FILE, OUT_FID] = SRILM_LM_OPEN(MODEL, [VAR1, VAL1]...)
%
%   MODEL should be a string denoting where to find a SRILM language model file
%   (path and name).  The file should be in the ARPA format SRILM expects.  This
%   model is returned in srilm_lm_init().
%
%   PIPE returned is handle to the input pipe setup to efficiently score
%   sequences of text lines.  This pipe should then be passed to srilm_lm_score
%
%   OUT_FILE returned will be a char array specifying the path and name of the
%   output file that the ngram stats will be written to by SRI-LM
%
%   OUT_FID returned is a descriptor to OUT_FILE (so we can just pass this to
%   the scoring routine to read the sequence probability)
%
%   NOTE: any default values defined below in the LOCAL VARS can be overridden 
%   by passing the name of parameter, and its new value on the command line.
% 
%   NOTE: We require that SRILM is installed, and the tool ngram is in the
%   users path
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: srilm_lm_open.m,v 1.2 2006-11-25 20:09:37 scottl Exp $
%
% REVISION HISTORY
% $Log: srilm_lm_open.m,v $
% Revision 1.2  2006-11-25 20:09:37  scottl
% bugfix to ensure the correct line gets read (see srilm_lm_score too)
%
% Revision 1.1  2006-11-22 17:11:13  scottl
% initial revision.
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

%default path and name of the output file
srilm_seq_file = '/tmp/tmp_srilm.txt';

%sequence of characters used to denote the end of a sequence of lines (and that
%log probability and perplexity score should be calculated
%NOTE: this value *must* match the escape_seq in srilm_lm_score
escape_seq = 'EOS';


% CODE START %
%%%%%%%%%%%%%%

if nargin < 1
    error('incorrect number of arguments');
elseif nargin > 1
    process_optional_args(varargin{:});
end

cmd = ['ngram -order ', order, ' -lm ', lm, ' -ppl - -escape ', escape_seq, ...
       ' ', other_params, ' -debug 0 '];
if ~isempty(input_vocab)
    cmd = [cmd, '-vocab ', input_vocab];
end
if skip_oovs
    cmd = [cmd, '-skipoovs'];
end
cmd = [cmd, ' > ', srilm_seq_file];
pipe = popenw(cmd);

if pipe < 0
    error('problems creating pipe handle via popenw');
end  

popenw(pipe, [escape_seq, char(10)], 'char');
pause(1);  %ensure the file will exist before we try and open it
srilm_fid = fopen(srilm_seq_file, 'r');

if srilm_fid == -1
    popenw(pipe, []);  %close the program
    delete(srilm_seq_file);
    error('problems opening output file: %s', srilm_seq_file);
end
