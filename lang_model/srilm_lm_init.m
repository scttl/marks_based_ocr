function [lm_file, vocab] = srilm_lm_init(Files, varargin)
% SRILM_LM_INIT Create a SRILM format language model
%
%   [LM_FILE, VOCAB_LIST] = SRILM_LM_INIT(FILES, [VAR1, VAL1]...)
%
%   FILES should either be a string or a cell array listing the path and name of
%   1 or more plain-text corpus files.  Since SRILM operates on word tokens, we 
%   assume that if characters are to be used as tokens, that the file has 
%   already been converted to have each character space-separated (and spaces 
%   mapped to another symbol)
%
%   LM_FILE returned is a string listing the path and name of the ARPA formatted
%   text file that can be used in subsequent language model tasks.
%
%   VOCAB_LIST returned is a cell array listing the symbols (or words) found 
%   in the copora passed.
%
%   NOTE: any default values defined below in the LOCAL VARS can be overridden 
%   by passing the name of parameter, and its new value on the command line.
% 
%   NOTE: We require that SRILM is installed, and the tools ngram-count and
%   ngram-merge are in the users path
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: srilm_lm_init.m,v 1.2 2006-11-22 17:11:13 scottl Exp $
%
% REVISION HISTORY
% $Log: srilm_lm_init.m,v $
% Revision 1.2  2006-11-22 17:11:13  scottl
% strip unknown, pause, and start and end of sentence tokens from vocab.
%
% Revision 1.1  2006-11-13 17:58:23  scottl
% initial check-in.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%create files relative to the MOCR_PATH
global MOCR_PATH;

%default language model file path and name
lm_file = [MOCR_PATH, '/results/lm.arpa'];

%what order n-grams should we create?
order = '3';

%what sort of discounting should we do?
discount = '';  %this defaults to good-turing backoff

%use a specific input vocabulary file?  Note that if doing so, the space 
%character should be mapped to some other character: like '_', and there should %be 1 character per line
input_vocab = '';

%any other parameters to pass to ngram-count?  see man ngram-count for details
other_params = '';

%tmp directory to store intermediate files
tmp_dir = '/tmp';

%temporary vocab file name
tmp_vocab = [tmp_dir, '/vocab.temp'];

%which vocabulary tokens should be stripped?  (present the list as a cell array)
vocab_strip = {'-pau-'; '<s>'; '</s>'; '<unk>'};


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

if length(Files) > 1
    %must merge corpora
    file_list = '';
    for ii=1:length(Files)
        tmp_file = [tmp_dir, sprintf('/tmp_nc.%d', ii)];
        cmd = ['ngram-count -sort -text ', Files{ii}, ' -write ', tmp_file];
        if ~isempty(input_vocab)
            cmd = [cmd, ' -limit-vocab -vocab ', input_vocab];
        end
        [s,w] = unix(cmd);
        if s ~= 0
            error('error running ngram count: %s', cmd);
        end
        file_list = [file_list, ' ', tmp_file];
    end
    cmd = ['ngram-merge -write ', lm_file, ' ', file_list];
    [s,w] = unix(cmd);
    if s ~= 0
        error('error running ngram merge: %s', cmd);
    end

    %cleanup the temp files
    cmd = ['rm -f ', file_list];
    [s,w] = unix(cmd);
    if s ~= 0
        error('problems deleting the following files: %s', file_list);
    end
    input = ['-read ', lm_file];
else
    input = ['-text ', Files{1}];
end
cmd = ['ngram-count -order ', order, ' ', discount, ' ', other_params, ...
       ' -write-vocab ', tmp_vocab, ' ', input, ' -lm ' lm_file];
if ~isempty(input_vocab)
   cmd = [cmd, ' -limit-vocab -vocab ', input_vocab];
end
[s,w] = unix(cmd);
if s ~= 0
    error('error running ngram-count: %s', cmd);
end

%read and store the vocabulary file contents
fid = fopen(tmp_vocab);
if fid == -1
    error('problems reading vocab file: %s', tmp_vocab);
end
vocab = textscan(tmp_vocab, '%s', 'delimiter', '\n', 'whitespace', '');
vocab = vocab{1};
cmd = ['rm -f ', tmp_vocab];
[s,w] = unix(cmd);
if s ~= 0
    error('error removing temp vocab file: %s', cmd);
end

%strip any unwanted vocabulary characters (like <unk> and others added by SRILM)
if ~isempty(vocab_strip)
    rem_idx = [];
    for ii=1:length(vocab_strip)
        rem_idx = [rem_idx; find(strcmp(vocab, vocab_strip{ii}))];
    end
    vocab = vocab(setdiff(1:length(vocab), rem_idx));
end
