function [vals, segs] = do_ocr(data, char_bitmaps, char_offsets, lang_model,pmt)
% DO_OCR   Run language-model based OCR on the line images of data passed
%
%   [vals, segs] = DO_OCR(data, char_bitmaps, char_offsets, lang_model,
%                  [prematch_thresh])
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
%   prematch_thresh is optional and if passed it will be used as a threhsold 
%   to pre-scan the line images, and attempt to do a straight Euclidean distance
%   match with char_bitmaps that fit within the threshold.  This constrains the
%   line solver, and should hopefully speed things up.
%
%   vals will either be an ASCII character array, or a cell array of char arrays
%   containing the sequence of recognized characters found in each image.
%
%   segs will either be a vector, or a cell array of vectors each entry 
%   specifying a cut point used to match characters
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: do_ocr.m,v 1.8 2006-08-30 17:37:32 scottl Exp $
%
% REVISION HISTORY
% $Log: do_ocr.m,v $
% Revision 1.8  2006-08-30 17:37:32  scottl
% implemented ability to pre-scan lines and place matches below a passed threshold
%
% Revision 1.7  2006/08/24 21:36:09  scottl
% fixed bug in descender matching (was off by 1 pixel).
%
% Revision 1.6  2006/08/14 16:39:42  scottl
% small fix to ensure matching datatypes are passed to imresize.
%
% Revision 1.5  2006/08/14 01:21:40  scottl
% updated default parameter settings based on some test performed.
%
% Revision 1.4  2006/08/05 17:31:54  scottl
% changed default minimum width.
%
% Revision 1.3  2006/07/05 00:57:33  scottl
% small fix to use length instead of size (remove dependence on data being a
% row vector)
%
% Revision 1.2  2006/06/19 21:48:29  scottl
% implemented baselines in character models and training case lines.
%
% Revision 1.1  2006/06/10 21:01:44  scottl
% Initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%
ins_prob = 1e-6;
del_prob = 1e-6;
min_window_width = 1;
max_window_width = 10;

%scan the lines for bitmap matches to constrain the line solver?
use_prematch = false;  
prematch_thresh = 0;  %only take perfect matches

return_cell = true; %return a cell structure unless there is only one input line

% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 4 || nargin > 5
    error('incorrect number of arguments specified!');
elseif nargin == 5
    use_prematch = true;
    prematch_thresh = pmt;
end

%imresize has problems when char_offsets are not doubles.
char_offsets = double(char_offsets);

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
    return_cell = false;
    data = {data};
end

num_cases = length(data);
vals = cell(0);
segs = cell(0);
prematch_count = 0;
for ii=1:num_cases
    fprintf('%.2fs: augmenting character bitmaps for line %d\n', toc, ii);
    aug_char_bitmaps = augment_bitmaps(char_bitmaps, char_offsets, ...
                       size(data{ii},1), data_offsets(ii));

    if use_prematch
        fprintf('%.2fs: attempting to pre-match bitmaps against line %d\n', ...
                toc, ii);
        norms = zeros(length(aug_char_bitmaps),1);
        for jj=1:length(aug_char_bitmaps)
            norms(jj) = sum(aug_char_bitmaps{jj}(:));
        end
        [initv,inits] = prematch_line(data{ii}, aug_char_bitmaps, norms, ...
                        prematch_thresh);
        prematch_count = prematch_count + length(initv);
        fprintf('%.2fs: found %d matches for this line\n', toc, length(initv));
        fprintf('%.2fs: attempting to solve line %d\n', toc, ii);
        %use the matched constraints above to repeatedly solve line chunks
        %between matches
        last_col = size(data{ii},2);
        while sum(data{ii}(:,last_col)) == 0
            last_col = last_col - 1;
        end
        if length(initv) > 0 && inits(2,end) < last_col
            %last match occurs before the end of the line.  Solve up to that
            %point
            [vals{ii,1},this_segs] = solveline(...
                                      data{ii}(:,inits(2,end)+1:end), ...
                                      aug_char_bitmaps, lang_model, ...
                                      ins_prob, del_prob, min_window_width, ...
                                      max_window_width);
            segs{ii,1} = this_segs + inits(2,end);
        else
            vals{ii,1} = [];
            segs{ii,1} = [];
        end
        while(length(initv) > 1)
            vals{ii,1} = [initv(end), vals{ii,1}];
            segs{ii,1} = [inits(1,end), segs{ii,1}];
            end_col = inits(1,end) - 1;
            start_col = inits(2,end-1) + 1;
            [vv,ss] = solveline(data{ii}(:, start_col:end_col), ...
                      aug_char_bitmaps, lang_model, ins_prob, del_prob, ...
                      min_window_width, max_window_width, initv(end));
            vals{ii,1} = [vv, vals{ii,1}];
            segs{ii,1} = [ss+start_col-1, segs{ii,1}];
            initv = initv(1:end-1);
            inits = inits(:,1:end-1);
        end
        if length(initv) == 1
            vals{ii,1} = [initv(1), vals{ii,1}];
            segs{ii,1} = [inits(1,1), segs{ii,1}]; 
            if inits(1,1) ~= 1
                %first match occurs before the first column
                [vv,ss] = solveline(data{ii}(:, 1:inits(1,1)-1), ...
                          aug_char_bitmaps, lang_model, ins_prob, del_prob, ...
                          min_window_width, max_window_width, initv(1));
                vals{ii,1} = [vv, vals{ii,1}];
                segs{ii,1} = [ss, segs{ii,1}];
            end
        else
            %we can reach this case if there weren't any prematches within the
            %threshold.  In such a case, just solve the entire line as normal
            [vals{ii,1},segs{ii,1}] = solveline(data{ii}, aug_char_bitmaps, ...
                                      lang_model, ins_prob, del_prob, ...
                                      min_window_width, max_window_width);
        end
    else
        fprintf('%.2fs: attempting to solve line %d\n', toc, ii);
        [vals{ii,1},segs{ii,1}] = solveline(data{ii}, aug_char_bitmaps, ...
        lang_model, ins_prob, del_prob, min_window_width, max_window_width);
    end
end

if ~ return_cell
    vals = vals{1};
    segs = segs{1};
end

if use_prematch
    fprintf('%.2fs: found %d prematches\n', toc, prematch_count);
end
fprintf('\n%.2fs: all lines completed\n',toc);

% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%augment_bitmaps - this function resizes the character bitmaps passed, placing
%them at their appropriate offset, and padding the necessary lines above and 
%below with blanks
function aug_maps = augment_bitmaps(bitmaps, bitmap_offs, line_height, baseline)

num_chars = length(bitmaps);

for ii=1:num_chars
    [char_h, char_w] = size(bitmaps{ii});
    l=1;
    r=char_w;
    if bitmap_offs(ii) == 0
        b=baseline+bitmap_offs(ii);
    else
        b=baseline+bitmap_offs(ii)+1;
    end
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
        aug_maps{ii} = zeros(height, char_w);
        aug_maps{ii}(t:b,l:r) = bitmaps{ii};
        aug_maps{ii} = imresize(aug_maps{ii}, [line_height, ...
                      ceil((line_height/height)*char_w)], 'nearest');
    else
        aug_maps{ii} = zeros(line_height, char_w);
        aug_maps{ii}(t:b,l:r) = bitmaps{ii};
    end
end


%prematch_line - this function takes an image line and augmented character
%bitmaps and returns vectors containing the model id's, and their left most
%position, for bitmaps that match the data within the threshold passed.
function [initv, inits] = prematch_line(line, bitmaps, norms, threshold)

initv = []; inits = [];

%since the bitmaps are as tall as the line, we should just need to pull out
%components from the data (based on a vertical projection sum), and calculate
%the Euclidean distance between these found components and each bitmap.
lsegs = find(sum(line) == 0);
lsegs = [0, lsegs, size(line,2)+1];

while(length(lsegs) > 1)
    start = lsegs(1) + 1;
    finish = lsegs(2) - 1;
    img_seg = line(:,start:finish);
    D = euc_dist(img_seg, bitmaps, sum(img_seg(:)), norms);
    %if there are multiple matches within the threshold, the first found lowest
    %match is returned.
    idx = find(D <= threshold,1);
    if ~isempty(idx)
        initv = [initv, idx];
        inits = [inits, [start; finish]];
    end
    lsegs = lsegs(2:end);
    while(length(lsegs) > 1 && lsegs(2) == lsegs(1) + 1)
        lsegs = lsegs(2:end);
    end
end


