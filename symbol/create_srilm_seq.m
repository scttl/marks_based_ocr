function create_srilm_seq(seq, file, varargin)
% CREATE_SRILM_SEQ Write the sequence out in SRILM format to the file passed
%
%   CREATE_SRILM_SEQ(SEQ, FILE, ['var1', new_val1]...)
%
%   SEQ should be a newline delimited string, or a cell array of strings
%   (each of which will be treated as a single line)
%
%   FILE should be a string listing the path and name of the file containing
%   the alphabet of output symbols to choose amongst during character 
%   recognition.  We assume the file is encoded in the default Locale (though
%   this should be able to be overridden)


% CVS INFO %
%%%%%%%%%%%%
% $Id: create_srilm_seq.m,v 1.1 2006-11-13 18:10:34 scottl Exp $
%
% REVISION HISTORY
% $Log: create_srilm_seq.m,v $
% Revision 1.1  2006-11-13 18:10:34  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%character to map spaces to
space_map = '_';


% CODE START %
%%%%%%%%%%%%%%

%process arguments
if nargin < 2
    error('must specify a file to process!');
elseif nargin > 2
    process_optional_args(varargin{:});
end

if iscell(seq)
    %convert to a large newline delimited string
    newseq = '';
    for ii=1:length(seq)
        if size(seq{ii},1) > 1
            seq{ii} = seq{ii}';
        end
        newseq = [newseq, seq{ii}, char(10)]; %newline
    end
    seq = newseq;
end

%replace any spaces with the appropariate character
seq(seq == ' ') = space_map;

%add spaces between each character so SRILM interprets them as words
if size(seq,1) > 1
    seq = seq';
end
seq = regexprep(seq, '([^\n])', '$1 ');

fid = fopen(file, 'w');
fprintf(fid,'%s', seq);
fclose(fid);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
