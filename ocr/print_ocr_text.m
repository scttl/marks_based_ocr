function text = print_ocr_text(line_nums, Comps, Syms, map, varargin)
%  PRINT_OCR_TEXT  Display the mapped OCR symbol text of the line number passed
%
%   print_ocr_text(line_nums, Comps, Syms, map)
%
%   line_nums should be a vector listing valid line numbers to be printed (in
%   the order they are to be printed).  Any invalid lines are skipped when
%   printed to the screen (an empty string is saved to text for that entry)
%
%   Comps is the comp struct like that returned in get_comps()
%
%   Syms is the symbol structure.  See create_alphabet()
%
%   map is a learned cluster to symbol index map vector.  For each cluster, the
%   index of the mapped symbol is returned.  This can be found by running
%   learn_mappings() or srilm_learn_mappings()
%
%   text returned will be a cell array, with one entry per valid line, each of
%   which will be a character array representing the characters found.  By 
%   default, the charcter sequences are also printed to the screen (though this 
%   can be overridden
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: print_ocr_text.m,v 1.1 2006-11-25 20:11:09 scottl Exp $
%
% REVISION HISTORY
% $Log: print_ocr_text.m,v $
% Revision 1.1  2006-11-25 20:11:09  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%set save_results to true to write the results to disk based on the params
%below it
save_results = false;
global MOCR_PATH;  %make use of the globally defined MOCR_PATH variable
save_file = [MOCR_PATH, '/results/ocr_text_res.txt'];

display_text = true;  %print the character sequences to the terminal?

%setting this to true prints a space between each character (useful for
%testing with SRI-LM's ngram utility)
add_spaces = false;

%leaving this empty does not convert characters to spaces, however you can
%specify particular characters to turn into spaces (useful for converting the
%'_' into space if modified for use in SRI-LM)
%the format for this parameter should be a regexp pattern
map_to_space_pattern = '[_]';


% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 4
    error('incorrect number of arguments specified!');
elseif nargin > 4
    process_optional_args(varargin{:});
end

text = cell(length(line_nums),1);

%first ensure each of the multi-char symbols are row vectors
for ii=1:Syms.num
    if size(Syms.val{ii},2) > 1
        Syms.val{ii} = Syms.val{ii}';
    end
end

%get the cluster sequences corresponding to the lines passed
seq = get_cluster_seq(Comps, line_nums);

for ii=1:length(seq)
    text{ii} = cell2mat(Syms.val(map(seq{ii})))';
    if ~isempty(text{ii}) && ~isempty(map_to_space_pattern)
        idx = regexp(text{ii}, map_to_space_pattern);
        text{ii}(idx) = ' ';
    end
    if ~isempty(text{ii}) && add_spaces 
        text{ii} = regexprep(text{ii}, '(.)', '$1 ')';
        %remove the additional space at the end
        text{ii} = text{ii}(1:end-1);
    end
end

if display_text
    fprintf('%s\n', text{:});
end

%save the results to disk if required.
if save_results
    fprintf('\n%.2f: writing ocr text results to disk\n', toc);
    [fid, err_msg] = fopen(save_file, 'w');
    if fid == -1
        error('problems opening save file %s.  Error:\n%s', save_file, err_msg);
    end
    fprintf(fid, '%s\n', text{:});
    fclose(fid);
end

fprintf('Elapsed time: %f\n', toc);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
