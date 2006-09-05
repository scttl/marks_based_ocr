function [vals, segs] = do_ocr(data, char_bitmaps, char_offsets, lang_model, ...
                        use_sl, pmt)
% DO_OCR   Run language-model based OCR on the line images of data passed
%
%   [vals, segs] = DO_OCR(data, char_bitmaps, char_offsets, lang_model,
%                  [use_shortlist, prematch_thresh])
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
%   use_shortlist is an optional boolean that when set to true, will limit the
%   number of character models considered for each column.  It defaults to
%   false.
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
% $Id: do_ocr.m,v 1.9 2006-09-05 15:53:56 scottl Exp $
%
% REVISION HISTORY
% $Log: do_ocr.m,v $
% Revision 1.9  2006-09-05 15:53:56  scottl
% implemented short-lists to constrain line solver.
%
% Revision 1.8  2006/08/30 17:37:32  scottl
% implemented ability to pre-scan lines and place matches below a passed 
% threshold
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
max_window_width = 15;

%scan the lines for bitmap matches to constrain the line solver?
use_prematch = false;  
prematch_thresh = -1;  %by default don't constrain the solver
max_sl_size = 5;  %max. number of short-list candidates to return.

return_cell = true; %return a cell structure unless there is only one input line

% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 4 || nargin > 6
    error('incorrect number of arguments specified!');
elseif nargin >= 5
    if use_sl
        use_prematch = true;
    else
        max_sl_size = length(char_bitmaps);
    end
    if nargin == 6
        use_prematch = true;
        prematch_thresh = pmt;
    end
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
        [candidates,nm] = prematch_line(data{ii}, aug_char_bitmaps, ...
                                        max_sl_size, prematch_thresh);
        prematch_count = prematch_count + nm;
        fprintf('%.2fs: found %d matches for this line\n', toc, nm);
        fprintf('%.2fs: attempting to solve line %d\n', toc, ii);
        [vals{ii,1},segs{ii,1}] = solveline(data{ii}, aug_char_bitmaps, ...
                                  lang_model, ins_prob, del_prob, ...
                                  min_window_width, max_window_width, ...
                                  candidates);
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
%bitmaps and returns a cell array containing vectors for each column of the
%image as well as a scalar giving the number of pre-matches.  Each vector 
%gives a list of possible candidate indices and may be a single candidate if it 
%is below or equal to the threshold passed.
function [candidates,nm] = prematch_line(line, bitmaps, max_models, threshold)

pct = .75;  %percentage of perfect pixel match for the model to the data

num_cols = size(line,2);
candidate_scores = Inf(max_models, num_cols);
candidates = zeros(max_models, num_cols);

for ii=1:length(bitmaps)
    bm_sum = sum(bitmaps{ii}(:));
    bm_width = size(bitmaps{ii},2);
    %perform a convolutional sweep across the data to find columns of large
    %overlap with this bitmap.  Since this returns the overlap about the center
    %of the bitmap we must shift the results to align over the left-most column.
    %Since the sweep doesn't give values for the last few columns, set their
    %score to a maximum to ensure Euclidean distance calculations are carried
    %out for those columns.
    scores = max(filter2(bitmaps{ii}, line));
    scores = [scores(ceil(bm_width/2):end), Inf(1,ceil(bm_width/2)-1)];
    %perform Euclidean distance calculations at the columns that have at least
    %pct of their total pixels matching
    if bm_sum == 0
        %this happens if we pass in a completely blank model
        chk_cols = 1:num_cols;
    else
        chk_cols = find(scores/bm_sum >= pct);
    end
    num_chk = length(chk_cols);
    sq_norms = zeros(num_chk,1);
    data_segs = cell(num_chk,1);
    for jj=1:num_chk
        start = chk_cols(jj);
        finish = start + bm_width -1;
        if finish > num_cols
            finish = num_cols;
        end
        data_segs{jj} = line(:,start:finish);
        sq_norms(jj) = sum(data_segs{jj}(:));  %no need to square since logical
    end
    scores(chk_cols) = euc_dist(bitmaps{ii}, data_segs, ...
                       sum(bitmaps{ii}(:).^2), sq_norms);

    row = 1;
    while length(chk_cols) > 0 && row <= max_models
        rplc_cols = find(scores(chk_cols) < candidate_scores(row, chk_cols));
        rplc_cols = chk_cols(rplc_cols);
        if row < max_models
            %push previous bests down 1 row to make room for the new scores
            candidate_scores(row+1:max_models, rplc_cols) = candidate_scores(...
                                               row:max_models-1, rplc_cols);
            candidates(row+1:max_models, rplc_cols) = candidates(...
                                               row:max_models-1, rplc_cols);
        end
        candidate_scores(row, rplc_cols) = scores(rplc_cols);
        candidates(row, rplc_cols) = ii;
        chk_cols = setdiff(chk_cols, rplc_cols);
        row = row+1;
    end
end

candidates = mat2cell(candidates, max_models, ones(1,num_cols));

%prune the candidate list for perfect matches (those within the threshold)
match_idcs = find(candidate_scores(1,:) <= threshold);
for ii=match_idcs
    candidates{ii} = candidates{ii}(1);
end
nm = length(match_idcs);

%also prune any non-matching candidates
row = 1;
chk_cols = 1:num_cols;
chk_cols = setdiff(chk_cols, match_idcs);
while length(chk_cols) > 0 && row <= max_models
    match_idcs = find(candidate_scores(row,chk_cols) == Inf);
    match_idcs = chk_cols(match_idcs);
    for jj=match_idcs
        candidates{jj} = candidates{jj}(1:row-1);
    end
    chk_cols = setdiff(chk_cols, match_idcs);
    row = row+1;
end
