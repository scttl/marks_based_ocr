function [vals, segs] = do_ocr(data, char_bitmaps, char_offsets, lang_model)
% DO_OCR   Run language-model based OCR on the line images of data passed
%
%   [vals, segs] = DO_OCR(data, char_bitmaps, char_offsets, lang_model)
%   data should either be a logical array, or a cell array of logical arrays
%   (one per row) each of which is an image representation of a sentence/line
%   upon which we will perform OCR.
%
%   char_bitmaps should be a cell array of logical arrays representing the 
%   individual characters in our alphabet (1 per cell entry).
%
%   char_offsets should be an array of values listing how far up/down the bottom
%   of the associated character bitmap should be placed from the base line.
%   Typically, most characters will have an offset of 0 (i.e. they should have 
%   their bottom aligned with the baseline, however characters like 'y' will
%   hang below the baseline).
%
%   lang_model should be a matrix representing character transition 
%   probabilities based on a language (like English).  It should have the
%   same number of rows as columns, and both dimensions should be the same
%   size as the char_bitmaps cell array.
%
%   vals will either be an ASCII character array, or a cell array of char arrays
%   containing the sequence of recognized characters found in each image.
%
%   segs will either be a vector, or a cell array of vectors each entry 
%   specifying a cut point used to match characters
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: do_ocr.m,v 1.2 2006-06-19 21:48:29 scottl Exp $
%
% REVISION HISTORY
% $Log: do_ocr.m,v $
% Revision 1.2  2006-06-19 21:48:29  scottl
% implemented baselines in character models and training case lines.
%
% Revision 1.1  2006/06/10 21:01:44  scottl
% Initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%
ins_prob = .1;  %1e-10;
del_prob = .05; %1e-10;
min_window_width = 1;
max_window_width = 4;



% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin ~= 4
    error('incorrect number of arguments specified!');
end

num_chars = length(char_bitmaps);
if any(size(lang_model) ~= [num_chars num_chars])
    error('char_bitmaps length does not match lang_model dimensions!');
end
if num_chars ~= length(char_offsets)
    error('char_bitmaps length does not match char_offsets length!');
end

%first get the baseline offsets of each training case (to realign and resize
%the char bitmaps)
fprintf('%.2fs: acquiring training data image baselines\n', toc);
data_offsets = get_baselines(data);

if ~ iscell(data)
    fprintf('%.2fs: augmenting character bitmaps for line 1\n', toc);
    aug_char_bitmaps = augment_bitmaps(char_bitmaps, char_offsets, ...
                       size(data,1), data_offsets);
    fprintf('%.2fs: attempting to solve line 1\n', toc);
    [vals, segs] = solveline(data, aug_char_bitmaps, lang_model, ...
    ins_prob, del_prob, min_window_width, max_window_width);
else
    num_cases = size(data,1);
    vals = cell(0);
    segs = cell(0);
    for i=1:num_cases
        fprintf('%.2fs: augmenting character bitmaps for line %d\n', toc, i);
        aug_char_bitmaps = augment_bitmaps(char_bitmaps, char_offsets, ...
                           size(data{i},1), data_offsets(i));
        fprintf('%.2fs: attempting to solve line %d\n', toc, i);
        [vals{i,1},segs{i,1}] = solveline(data{i}, aug_char_bitmaps, ...
           lang_model, ins_prob, del_prob, min_window_width, max_window_width);
    end
end

fprintf('\n%.2fs: all lines completed\n');

% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%augment_bitmaps - this function resizes the character bitmaps passed, placing
%them at their appropriate offset, and padding the necessary lines above and 
%below with blanks
function aug_maps = augment_bitmaps(bitmaps, bitmap_offs, line_height, baseline)

num_chars = length(bitmaps);

for i=1:num_chars
    [char_h, char_w] = size(bitmaps{i});
    l=1;
    r=char_w;
    b=baseline+bitmap_offs(i);
    t=b-char_h+1;

    if t < 1  || b > line_height %bitmap at its offset is taller than the line
        %@@ how to handle this?  For now just resize the bitmap
        bottom = max(b, line_height);
        top = 1;
        if t < 1
            height = bottom + 1 + (-t);
            t = 1;
            b = char_h;
        else
            height = bottom;
        end
        aug_maps{i} = zeros(height, char_w);
        aug_maps{i}(t:b,l:r) = bitmaps{i};
        aug_maps{i} = imresize(aug_maps{i}, [line_height, ...
                      ceil((line_height/height)*char_w)], 'bilinear');
    else
        aug_maps{i} = zeros(line_height, char_w);
        aug_maps{i}(t:b,l:r) = bitmaps{i};
    end
end


