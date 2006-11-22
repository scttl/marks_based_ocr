function [total_logp, total_ppl, num_s, num_sym, num_oov] = ...
         srilm_lm_close(pipe, file, varargin)
% SRILM_LM_CLOSE Close the pipe to SRI-LM's scoring utility
%
%   [TOTAL_LOGP, TOTAL_PPL, NUM_SENT, NUM_SYMS, NUM_OOV] = 
%   SRILM_LM_CLOSE(PIPE, OUTPUT_FILE, [VAR1, VAL1]...)
%
%   PIPE should be a descriptor to a running version of SRI-LM's ngram utility.
%
%   OUTPUT_FILE should be the path and name of the file that SRI-LM is writing
%   its output to.
%
%   TOTAL_LOGP is the total log probability of all the sequences examined thus
%   far.
%
%   TOTAL_PPL is a scalar giving the total perplexity taken overa ll the
%   sequences examined thus far.
%
%   NUM_SENT will be a scalar specifying the total number of sentences seen
%
%   NUM_SYMS will be a scalar specifying the total number of symbols seen
%
%   NUM_OOVS will be a scalar specifying the total number of out-of-vocabulary
%   i.e. unknown symbols seen
%
%   NOTE: This method is intended to be called after calling srilm_lm_open
%   (to get the pipe handle and output file name)
% 
%   NOTE: We require that SRILM is installed, and the tool ngram is in the
%   users path
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: srilm_lm_close.m,v 1.1 2006-11-22 17:11:13 scottl Exp $
%
% REVISION HISTORY
% $Log: srilm_lm_close.m,v $
% Revision 1.1  2006-11-22 17:11:13  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%by default, delete the output file
rm_output_file = true;

total_logp = NaN;
total_ppl = NaN;
num_s = NaN;
num_sym = NaN;
num_oov = NaN;


% CODE START %
%%%%%%%%%%%%%%

if nargin < 2
    error('incorrect number of arguments');
elseif nargin > 2
    process_optional_args(varargin{:});
end

%close the pipe to the command (this will flush final stats to the file)
popenw(pipe, []);

%open and read the last two lines of the file (these contain the stats we need)
fid = fopen(file);
if fid == -1
    warning('MBOCR:CantOpenFile', ...
            'Problems opening: %s, no stats collected\n', file);
else
    %first get sentence, symbol, and oov counts
    D = textscan(fid, 'file%*s%d%*s%d%*s%d%*[^\n]');
    num_s = D{1};
    num_sym = D{2};
    num_oov = D{3};

    %now get logp and ppl counts
    frewind(fid);
    D = textscan(fid, '%*dzeroprobs, logprob=%fppl=%f%*[^\n]');
    if ~isempty(D) && ~isempty(D{1}) && ~isempty(D{2})
        total_logp = D{1}(end);
        total_ppl = D{2}(end);
    end
end
fclose(fid);

%during the srilm_lm_open call we write a dummy word to the sequence to ensure
%buffering works correctly.  Subtract these values from the counts
num_s = num_s - 1;
num_sym = num_sym - 1;
num_oov = num_oov - 1;

%delete the output file if required
if rm_output_file
    delete(file);
end



