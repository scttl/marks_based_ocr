function [vals, segs] = do_ocr(data, char_bitmaps, lang_model)
% DO_OCR   Run language-model based OCR on the line images of data passed
%
%   [vals, segs] = DO_OCR(data, char_bitmaps, lang_model)
%   data should either be a logical array, or a cell array of logical arrays
%   (one per row) each of which is an image representation of a sentence/line
%   upon which we will perform OCR.
%
%   char_bitmaps should be a cell array of logical arrays representing the 
%   individual characters in our alphabet (1 per cell entry).
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
% $Id: do_ocr.m,v 1.1 2006-06-10 21:01:44 scottl Exp $
%
% REVISION HISTORY
% $Log: do_ocr.m,v $
% Revision 1.1  2006-06-10 21:01:44  scottl
% Initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%
ins_prob = 1e-10;
del_prob = 1e-10;
min_window_width = 1;
max_window_width = 4;



% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin ~= 3
    error('incorrect number of arguments specified!');
end

num_chars = length(char_bitmaps);
if any(size(lang_model) ~= [num_chars num_chars])
    error('char_bitmaps length does not match lang_model dimensions!');
end

if iscell(data)
    num_cases = size(data,1);
else
    num_cases = 1;
end

if num_cases == 1
    fprintf('attempting to solve line 1\n');
    [vals, segs] = solveline(data, char_bitmaps, lang_model, ...
    ins_prob, del_prob, min_window_width, max_window_width);
else
    vals = cell(0);
    segs = cell(0);
    for i=1:num_cases
        fprintf('%.2fs: attempting to solve line %d\n', toc, i);
        [vals{i,1},segs{i,1}] = solveline(data{i}, char_bitmaps, lang_model, ...
                   ins_prob, del_prob, min_window_width, max_window_width);
    end
end

fprintf('\n%.2fs: all lines completed\n');

% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
