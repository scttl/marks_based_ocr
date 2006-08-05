function [vals, imgs] = create_synth_training_data(Files, chars, num)
% CREATE_SYNTH_TRAINING_DATA  Parse & create images of ASCII textfile sentences.
%
%   [vals,imgs] = CREATE_SYNTH_TRAINING_DATA(Files, chars, [num])
%   Files should be a cell array (or single string), and give the full 
%   path and name of the file(s) to be used for creating the training data.
%
%   chars is a char arry giving the valid ASCII characters
%
%   num is optional and if specified gives the maximum number of new training
%   cases to create (there may be fewer than num cases created if the Files
%   contain fewer sentences).  Set to Inf to guarantee processing the entire
%   file.
%
%   vals is a cell array with one sentence (a char array) per row
%
%   imgs is the corresponding cell array containing image versions of the 
%   sentences in vals (stored as logical arrays, with a value of 1 representing
%   'on' pixels, and 0 'off').
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: create_synth_training_data.m,v 1.2 2006-08-05 17:39:33 scottl Exp $
%
% REVISION HISTORY
% $Log: create_synth_training_data.m,v $
% Revision 1.2  2006-08-05 17:39:33  scottl
% small grammatical fix.
%
% Revision 1.1  2006/06/10 21:01:47  scottl
% Initial revision.
%

% LOCAL VARS %
%%%%%%%%%%%%%%
max_cases = inf;
num_files = size(Files,1);

%this variable defines the regular expression used to pick out the end of
%each sentence.  We want to find the appropriate punctuation character followed
%by at least one whitespace character, followed by a capital letter
sentence_delim_pat = '[\.\?!]\s+[A-Z]';

%this variable defines the type of characters to keep.  It will be passed to
%the unix tr utility so it should match that format
%tr_keep_list = '[:lower:][:upper:] ';  %keep letters and the space char
tr_keep_list = chars;

%don't include any sentences that are too short or long (in characters)
min_sentence_len = 10;
max_sentence_len = 95;  %this ensure it fits on a single line when in TeX

%these parameters control how noisy the images become.  See add_img_noise
%also, setting rem_pct to a number less than 0, disables adding noise to the 
%image
rem_pct = .01;
add_pct = .01;
thresh  = .5;
fspec_args = {'gaussian', [3 3], .5};


% CODE START %
%%%%%%%%%%%%%%
tic;
if nargin < 2 || nargin > 3
    error('incorrect number of arguments specified');
elseif nargin == 3
    max_cases = num;
end

if num_files == 1
    fid = fopen(Files);
    if fid == -1
        error('unable to open file: %s', Files);
    end
    S = fread(fid);
    fclose(fid);
else
    for i = 1:num_files
        fid = fopen(Files{i});
        if fid == -1
            error('unable to open file: %s', Files{i});
        end
        S = [S; fread(fid)];
        fclose(fid);
    end
end
S = char(S)';
fprintf('%.2fs: finished reading all files\n', toc);

[endidx, startidx] = regexp(S, sentence_delim_pat, 'start', 'end');
startidx = [1, startidx];
endidx = [endidx, length(S)];

fprintf('%.2fs: generating sentences\n', toc);
vals = cell(0);
count = 0;
while length(endidx) > 0 && count < max_cases
    len = endidx(1) - startidx(1);
    if len <= max_sentence_len && len >= min_sentence_len
        count = count + 1;
        fprintf('  item: %d\r', count);
        vals{end+1,1} = S(startidx(1):endidx(1));
        %now strip out anything not on the keep list, converting 
        %newlines and tabs to spaces.  If the sentence contains a double
        %quote character we have to escape it.
        quote_idx = strfind(vals{end,1}, '"');
        if quote_idx
            num_ins = 0;
            while quote_idx
                pos = quote_idx(1) + num_ins;
                vals{end,1} = [vals{end,1}(1:pos-1), '\', vals{end,1}(pos:end)];
                num_ins = num_ins + 1;
                quote_idx = quote_idx(2:end);
            end
        end

        [s,vals{end,1}] = unix(['echo "', vals{end,1}, ...
           '" | tr -s [:space:] '' '' | ', 'tr -d -c ''', tr_keep_list, '''']);
        if s ~= 0
            fprintf(S(startidx(1):endidx(1)));
            error('problem converting sentence: %s', vals{end,1});
        end
    end
    startidx = startidx(2:end);
    endidx = endidx(2:end);
end
clear S;
fprintf('\n%.2fs: finished processing sentences\n', toc);

%now create images of each sentence
for i = 1:count
    fprintf('  item: %d\r', i);
    %imgs{i,1} = build_bitmap_image(vals{i}, names, bitmaps);
    imgs{i,1} = build_tex_image(vals{i});
end
%we must ensure that all images are of the same height, so add zeros to the
%bottom of shorter rows, and chop off any leading rowspace
first_row = inf;
for i = 1:count
    row = 1;
    while ~ any(imgs{i}(row,:))
        row = row + 1;
    end
    first_row = min(first_row, row);
end
for i = 1:count
    imgs{i} = imgs{i}(first_row:end,:);
end
max_height = 0;
for i = 1:count
    max_height = max(size(imgs{i},1), max_height);
end
for i = 1:count
    if size(imgs{i},1) ~= max_height
        imgs{i}(end+1:max_height,:) = 0;
    end
end
fprintf('\n%.2fs: finished creating images\n', toc);

%add noise to the images (using the GUI)
[rem_pct, add_pct, thresh, fspec_args] = noise_gui(imgs{1});
if rem_pct >= 0
    imgs = add_img_noise(imgs, rem_pct, add_pct, thresh, fspec_args{:});
    fprintf('%.2fs: finished adding noise to images\n', toc);
end



% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%build_tex_image: given an input string, create a TeX formatted file around it,
%compile it, and pass it to the dvipng utility to create an image representation
%then load that image back into matlab, and clean up created files
%NOTE: currently cannot change font etc.
function img = build_tex_image(str)

tex_file='/tmp/tmp_ocr.tex';
log_file='/tmp/tmp_ocr.log';
dvi_file='/tmp/tmp_ocr.dvi';
png_file='/tmp/tmp_ocr.png';
png_dpi='300';

img = [];

fid = fopen(tex_file, 'w');
fprintf(fid, '\\nopagenumbers\n%s\n\\end', str);
fclose(fid);

[s,w] = unix(['tex -output-directory=/tmp ', tex_file]);
if s ~= 0
    unix(['rm -f ' tex_file, ' ', dvi_file, ' ', log_file]);
    error('problem running TeX: %s', w);
end

[s,w] = unix(['dvipng -o ', png_file, ' ', ' -D ', png_dpi, ' ' dvi_file]);
if s ~= 0
    unix(['rm -f ' tex_file, ' ', dvi_file, ' ', log_file, ' ', png_file]);
    error('problem running dvipng: %s', w);
end

img = (imread(png_file) > 0);
[r,c] = find(img == 1, 1, 'first');
img = img(:,c:end);
[r,c] = find(img == 1, 1, 'last');
img = img(:,1:c);
unix(['rm -f ' tex_file, ' ', dvi_file, ' ', log_file, ' ', png_file]);


%build_bitmap_image: given an input string as well as character and associated 
%bitmaps, convert the string into an image representation
function img = build_bitmap_image(str, names, bitmaps)
img = [];
for i = 1:length(str)
    idx = strfind(names, str(i));
    img = [img, bitmaps(idx)];
end
img = cell2mat(img);
