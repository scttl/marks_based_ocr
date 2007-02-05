function Syms = create_alphabet(file, varargin)
% CREATE_ALPHABET   Initialize a struct containing output alphabet symbols
%
%   SYMS = CREATE_ALPHABET(FILE, ['var1', new_val1]...)
%   FILE should be a string listing the path and name of the file containing
%   the alphabet of output symbols to choose amongst during character 
%   recognition.  We assume the file is encoded in the default Locale (though
%   this should be able to be overridden)
%
%   SYMS is a struct meant to hold information on each output symbol.
%   It contains the following fields:
%     num - integer specifying the largest symbol index number.
%     input_file - this specifies the path and name of the file used to find
%                  the symbols
%     encoding - character encoding used for the symbols
%     val - this cell array will hold string representations of each encoded
%           symbol.
%     img - this is optionally created (depending on parameter settings below),
%           and will store a template image representation of the symbol.  The
%           size and font used are defined below.  If present this is used to
%           initially guide initial mappings between cluster images and
%           symbols.
%     use_srilm - this field determines how we build or n-gram model (if we
%                 build one).  Defaults to false
%     corpus_files - this is optionally created as a cell array, and used to
%                    store the list of text files read to determine bigram
%                    counts and words etc.
%     NOTE: we assume the symbol alphabet is a superset of the symbols found in
%           the corpora (otherwise an error is returned)
%     count - a vector of values listing how many times each symbol value
%             appears in the corpus.
%     pos_count - A cell array, each entry of which lists how many times each
%                 symbol occurs in each position of words of length the same as
%                 the index of the entry.  Thus the 10th entry lists positional
%                 counts for each symbol in words of length 10
%     pos_total - A count listing the total number of occurences of each 
%                 character in each position (up to the dimensions of pos_count)
%     trigram - A distribution over which symbol is likely to follow the
%               previous two symbols based on counts in the corpora (stored as
%               a 3 dimensional matrix)
%     bigram - A distribution over which symbols are likely to follow other
%              symbols based on counts in the corpora (stored as a matrix).
%              The entries are (possibly smoothed) counts listing how many
%              times one character was seen following another.
%     first_bg_count - a matrix lsting how many times each symbol follows the
%                      start symbol at the beginning of a line (possibly 
%                      smoothed)
%     first_count - a vector values listing how many times each symbol value
%                   appears at the start of a line (possibly smoothed)
%     words - a cell array of all the words grouped together from symbols in
%             the alphabet, and found in the corpora
%     word_count - a vector listing how many times each word appears in the
%                  corpora
%     srilm_file - if use_srilm is set to true, we estimate our n-gram counts
%                  and save them on disk in srilm_file
%     class - categorizes each symbol according to its type


% CVS INFO %
%%%%%%%%%%%%
% $Id: create_alphabet.m,v 1.10 2007-02-05 22:16:57 scottl Exp $
%
% REVISION HISTORY
% $Log: create_alphabet.m,v $
% Revision 1.10  2007-02-05 22:16:57  scottl
% added density field.
%
% Revision 1.9  2007-02-01 17:59:27  scottl
% replaced specific class type fields with a single general purpose
% class field.
%
% Revision 1.8  2007-01-30 01:37:49  scottl
% added grouping based on ascender and descender offsets.
%
% Revision 1.7  2007-01-25 18:49:09  scottl
% changed normalization value.
%
% Revision 1.6  2007-01-08 22:11:35  scottl
% added positional count information to each symbol.
%
% Revision 1.5  2006-12-01 23:02:39  scottl
% implemented ability to use TeX and dvipng to generate character images.
%
% Revision 1.4  2006-11-22 16:59:45  scottl
% updates to handle ligatures (multi-char symbols)
%
% Revision 1.3  2006-11-14 22:50:41  scottl
% changed using character values to strings (to accomodate ligatures etc)
%
% Revision 1.2  2006-11-13 18:10:09  scottl
% added ability to use SRILM based files, and internal trigram models.
%
% Revision 1.1  2006-10-29 17:06:21  scottl
% initial revision
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%set default encoding to this machines locale
encoding = feature('DefaultCharacterSet');  

%by default we don't generate templates.  Change this to a string giving the 
%fontname (or full path to a pk font) to generate templats.  Or one can use the
%string specified by tex_pattern below to use TeX to render images in its
%default font
template_font = []; 

%in what point size should the templates be generated.  This has no effect if
%templates are not being generated, or if we are using a pk file in a
%particular DPI
font_ptsize = 12;

%how can we tell if a template font passed is the path to a pk font?
pk_pattern = 'pk$';  %match *pk

%how can we tell if a template font passed requires us to use TeX?
tex_pattern = 'USE_TEX';

%by default we don't process any corpora looking for word and character counts
corpora_files = {};

%by default estimate the model internally
use_srilm = false;

%which order n-gram model should we estimate if using SRILM
order = 3;

%which character should we map space to if srilm being used?
srilm_space_map = '_';
srilm_vocab_file = '/tmp/srilm_input.vocab';

%the following parameters are only set if we aren't using SRILM
%if processing corpora, how many pseudocounts should we add to each trigram
%entry?
tg_pseudo = 1;

%if processing corpora, how many pseudocounts should we add to each bigram
%entry?
bg_pseudo = 1;

%if processing corpora, how many pseudocounts should we add to each first
%two characters in line appearance?
fc_bg_pseudo = 1;

%if processing corpora, how many pseudocounts should we add to each first
%character in line appearance?
fc_pseudo = 1;

%what font and size should we use for creating densities?
density_font = 'helvetica';
density_font_size = '36';

%when creating positional counts, up to what length word should be included?
max_word_len = 10;


% CODE START %
%%%%%%%%%%%%%%
tic;

%process arguments
if nargin < 1
    error('must specify a file to process!');
elseif nargin > 1
    process_optional_args(varargin{:});
end

fid = fopen(file,'r','n');
if fid == -1
    error('cannot open file %s', file);
end
Syms = init_symbols(file, encoding);

%read all input symbols (including whitespace, but ignoring trailing newline)
sym_list = textscan(fid, '%s','delimiter', '\n', 'whitespace', '');
%this creates a single cell containing the cell array
sym_list = sym_list{1};
fclose(fid);

Syms.val = sym_list;
Syms.num = length(Syms.val);

if ~isempty(template_font)
    fprintf('%.2fs: generating template images\n', toc);

    if ~isempty(regexp(template_font, pk_pattern))
        %attempt to generate templates using pk2bm
        Syms.img = generate_pk_templates(template_font, 'charnames', sym_list);
    elseif ~isempty(regexp(template_font, tex_pattern))
        %use TeX to render image templates (note that these won't be scaled,
        %and we must mutiply the pointsize by 10 to get the appropriate dpi for
        %dvipng
        Syms.img = cell(Syms.num,1);
        for ii=1:Syms.num
            if strcmp(Syms.val{ii}, ' ')
                %just build a dummy blank image
                Syms.img{ii} = zeros(font_ptsize);
            elseif ~isempty(regexp(Syms.val{ii}, '[\{\}\%&\$#_\^~\\]'))
                %have to treat these characters, since they are special in TeX
                %use \char and the number (need to get numbers)
            else
                Syms.img{ii} = build_tex_image(Syms.val{ii}, 'png_dpi', ...
                               font_ptsize*10);
            end
        end
    else
        %use imagemagick to render image templates in the font passed
        Syms.img = generate_templates(template_font, sym_list, ...
                   'ptsize', font_ptsize);
    end
end

if ~isempty(corpora_files)
    fprintf('%.2fs: reading text corpus\n', toc);
    if use_srilm
        Syms.use_srilm = true;
        %first create the temporary input vocab file and convert any spaces to 
        %the appropriate map character
        idx = find(strcmp(Syms.val, ' '));
        if ~isempty(idx)
            Syms.val{idx} = srilm_space_map;
        end
        fid = fopen(srilm_vocab_file, 'w');
        if fid == -1
            error('problem creating srilm vocab file');
        end
        %ensure that the srilm vocabulary only contains single characters (not
        %ligatures since these will not be found)
        for ii=1:Syms.num
            if length(Syms.val{ii}) == 1
                fprintf(fid, '%s\n', Syms.val{ii});
            end
        end
        fclose(fid);

        %now create the language model file
        Syms.srilm_file = srilm_lm_init(corpora_files, 'order', order, ...
                          'input_vocab', srilm_vocab_file);

        %cleanup the temp vocab file
        delete(srilm_vocab_file);
    else
        D = create_word_dictionary(corpora_files, 'max_word_len', max_word_len);
    
        if length(D.char) > Syms.num
            error('more symbols in corpus than in specified symbol files');
        end
        Syms.count = zeros(Syms.num, 1, 'uint32');
        Syms.pos_count = cell(size(D.pos_count));
        Syms.pos_total = 0;
        Syms.bigram = zeros(Syms.num, Syms.num, 'uint32');
        Syms.trigram = zeros(Syms.num, Syms.num, Syms.num, 'uint32');
        Syms.first_count = zeros(Syms.num, 1, 'uint32');
        Syms.first_bg_count = zeros(Syms.num, Syms.num, 'uint32');
    
        %must find the mapping between characters found in the corpus, and our
        %symbol characters
        idx = zeros(length(D.char),1);
        for ii=1:length(D.char)
            match_idx = find(strcmp(Syms.val, D.char(ii)));
            if length(match_idx) > 1
                warning('MBOCR:MultMatch', 'Found multiple matches for %c', ...
                D.char(ii));
                idx(ii) = match_idx(1);
            elseif length(match_idx) == 1
                idx(ii) = match_idx;
            else
                error('symbol: "%c" in corpus, not in symbol list', D.char(ii));
            end
        end
        Syms.count(idx) = D.char_count;
    
        %count the total number of occurences of each char in each position
        for ii=1:length(D.pos_count)
            Syms.pos_total = Syms.pos_total + sum(D.pos_count{ii}(:));
        end
        %Syms.pos_total = sum(Syms.count);

        %now add the normalized positional counts
        for ii=1:length(D.pos_count)
            Syms.pos_count{ii}(idx,:) = D.pos_count{ii} ./ Syms.pos_total;
        end

        %smooth the bigram by adding pseudo-counts
        Syms.bigram(idx,idx) = D.char_bigram;
        Syms.bigram = Syms.bigram + bg_pseudo;
    
        %smooth the trigram by adding pseudo-counts
        Syms.trigram(idx,idx,idx) = D.char_trigram;
        Syms.trigram = Syms.trigram + tg_pseudo;

        %smooth the start counts by adding pseudo-counts
        Syms.first_count(idx) = D.first_count;
        Syms.first_count = Syms.first_count + fc_pseudo;
        Syms.first_bg_count(idx,idx) = D.first_bg_count;
        Syms.first_bg_count = Syms.first_bg_count + fc_bg_pseudo;
    
        Syms.words = D.word;
        Syms.word_count = D.word_count;
    end
    Syms.corpus_files = corpora_files;
end

%now assign each value to the appropriate class
fprintf('%.2fs: assigning symbols to their respective classes\n', toc);
Syms.class = assign_class(Syms.val);

%now calculate densitys for each symbol
Syms.density = assign_density(char(Syms.val), 'img_font', density_font, ...
               'img_font_sz', density_font_size);

fprintf('%.2fs: Symbols initialized\n', toc);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this function creates an empty symbol stucture
function Syms = init_symbols(file, enc)

Syms.num = uint16(0);
Syms.input_file = file;
Syms.encoding = enc;
Syms.val = cell(0);
Syms.img = cell(0);
Syms.corpus_files = cell(0);
Syms.use_srilm = false;
Syms.count = uint32([]);
Syms.pos_count = cell(0);
Syms.pos_total = 0;
Syms.trigram = uint32([]);
Syms.bigram = uint32([]);
Syms.first_bg_count = uint32([]);
Syms.first_count = uint32([]);
Syms.words = cell(0);
Syms.word_count = uint32([]);
Syms.srilm_file = '';
Syms.class = uint16([]);
Syms.density = [];
