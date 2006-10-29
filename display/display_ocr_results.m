function display_ocr_results(indices, start_pos, bitmaps, base_offs, num, ...
         scores, imgs)
%  DISPLAY_OCR_RESULTS  Construct an image of the OCR'd indices passed
%
%   display_ocr_results(indices,start_pos,bitmaps,base_offs,[num,scores,imgs])
%
%   indices should either be a vector of index values into the bitmaps cell 
%   array, or it should be a cell array of index vectors (one per row).  Each
%   entry is assumed to lie within the length of the bitmaps array.  Multiple
%   vectors are drawn in order (from top-to-bottom).
%
%   start_pos should either be a vector of column positions at which to start
%   placing each bitmap in the corresponding indices entry, or it should be a
%   cell array of such vectors (one per row).
%
%   bitmaps should be a cell array of logical arrays, each of which represents
%   a single character in the alphabet, which are used to construct the output
%   image.
%
%   base_offs should be a vector listing the offset relative to the baseline of
%   the bottom of the corresponding bitmap.  This ensure characters are lined
%   up correctly.
%
%   num is optional and if specified, determines the number of OCR'd
%   lines to display.  If not specified, all indices are drawn.
%
%   scores is optional and if specified (and the same length as indices) will
%   be used to create a "heat map" which will be used to colour the characters
%   based on their score.
%
%   imgs is optional and if specified, will draw the original image line
%   directly above the OCR'd result line for comparison purposes.  If not
%   present, this will be left out.  There should be at least as many images
%   passed as num.  And both should be at least as large as the length of the
%   indices cell array.

% CVS INFO %
%%%%%%%%%%%%
% $Id: display_ocr_results.m,v 1.11 2006-10-29 17:20:54 scottl Exp $
%
% REVISION HISTORY
% $Log: display_ocr_results.m,v $
% Revision 1.11  2006-10-29 17:20:54  scottl
% small change in parameter name.
%
% Revision 1.10  2006-09-18 21:01:28  scottl
% no change.
%
% Revision 1.9  2006/09/16 22:30:32  scottl
% implemented heat maps based on passed scores.
%
% Revision 1.8  2006/09/05 15:50:45  scottl
% made use of MOCR_PATH variable for saving in the results directory.
%
% Revision 1.7  2006/08/24 21:13:17  scottl
% fix bug in display routine.  Use of logical or meant display of any non-zero
% pixel in the cluster average.
%
% Revision 1.6  2006/08/14 17:37:45  scottl
% added ability to interleve original image amongst OCR results.
%
% Revision 1.5  2006/08/14 01:32:09  scottl
% used column segments to space characters as per the OCR results.
%
% Revision 1.4  2006/08/05 17:35:29  scottl
% added ability to display segment points, removed dependence on imview.
%
% Revision 1.3  2006/07/05 01:07:46  scottl
% don't save results by default.
%
% Revision 1.2  2006/06/19 21:02:44  scottl
% aligned result bitmaps based on baseline offsets.
%
% Revision 1.1  2006/06/10 21:01:32  scottl
% Initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%set save_averages to true to write the results to disk based on the params
%below it
save_results = false;
global MOCR_PATH;  %make use of the globally defined MOCR_PATH variable
img_prefix = [MOCR_PATH, '/results/ocr_res'];
img_format = 'png';

display_segments = true;  %draw the segement lines in a different colour?
segment_col = reshape([255,0,0],1,1,3);  %this is red

row_margin = 10;  %number of pixels between consecutive rows in the image

%interleave original image line between OCR lines?  requires images to be
%passed as a parameter
display_original = false;

%use heat map based on passed scores to colour characters? requires scores to
%be passed as a paramter
map_scores = false; 
map_colormap = cool(64);  %type: help colormapeditor for other choices
map_stddev = 0.5;  %number of standard devs of data to map.  Set to 0 to map all
map_length = size(map_colormap,1);

on_thresh = 0.5;  %at what point are Cluster average intensities considered 'on'


% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 4 || nargin > 7
    error('incorrect number of arguments specified!');
elseif nargin >= 5
    num_rows = num;
    if num_rows > size(indices,1)
        num_rows = size(indices,1);
    end
    if nargin >= 6
        map_scores = true;
        if nargin == 7
            display_original = true;
            if ~ iscell(imgs)
                imgs = {imgs};
            end
        end
    end
else
    num_rows = size(indices,1);
end

if length(bitmaps) ~= length(base_offs)
    error('number of bitmaps must be the same as the number of base_offs!');
elseif size(bitmaps,1) > 1
    %convert to a column array
    bitmaps = bitmaps';
end

if display_original
    num_rows = 2*num_rows;
end
M = cell(num_rows,1);

if ~ iscell(indices)
    indices = {indices};
end
if ~ iscell(start_pos)
    start_pos = {start_pos};
end
if map_scores 
    if ~ iscell(scores)
        scores = {scores};
    end
    if any(size(indices) ~= size(scores))
        error('scores passed must be the same dimensions as the indices');
    else
        for ii=1:length(indices)
            if length(indices{ii}) ~= length(scores{ii})
                error('num of scores for row %d mismatches num of indices', ii);
            end
        end
    end
    all_scores = cell2mat(scores');
    max_score = max(all_scores);
    min_score = min(all_scores);
    if map_stddev ~= 0
        avg_score = mean(all_scores);
        std_score = std(all_scores);
        min_score = max(min_score, (avg_score - map_stddev*std_score));
        max_score = min(max_score, (avg_score + map_stddev*std_score));
    end
end

max_width = 0;
for ii=1:num_rows
    val = ceil(ii/2);
    if display_original && rem(ii,2) ~= 0
        %add the original image
        M{ii} = imgs{val};
    else
        if ~ display_original
            val = ii;
        end
        mm = bitmaps(indices{val});
        moffs = base_offs(indices{val});
        sp = start_pos{val};
        if map_scores
            ss = scores{val};
        end
        max_height = 0;
        for jj=1:length(mm)
            max_height = max(max_height, size(mm{jj},1)-moffs(jj));
        end
        max_base = max(moffs);
        for jj=1:length(mm)
            [h w] = size(mm{jj});
            mm{jj} = [zeros(max_height - (h-moffs(jj)), w); ...
                     mm{jj} > on_thresh; ...
                     zeros(max_base - moffs(jj), w)];
            if map_scores
                if ss(jj) >= max_score
                    mm{jj} = map_length * mm{jj};
                elseif ss(jj) <= min_score
                    mm{jj} = 1 * mm{jj};
                else
                    mm{jj} = ceil((ss(jj)-min_score)/(max_score-min_score) ...
                                  * map_length) * mm{jj};
                end
            end
        end
        M{ii} = [];
        for jj=1:length(mm)
            w = size(M{ii},2);
            extra_width = sp(jj) + size(mm{jj},2) - 1 - w;
            if extra_width > 0
                M{ii} = [M{ii}, zeros(max_height+max_base,extra_width)];
            end
            if map_scores
                %since larger scores get larger label values, we can take the
                %max when displaying overlapping mismatches i.e. display the
                %'hotter' or more mismatched for overlapping pixels
                M{ii}(:,sp(jj):sp(jj)+size(mm{jj},2)-1) = max(mm{jj}, ...
                      M{ii}(:,sp(jj):sp(jj)+size(mm{jj},2)-1));
            else
                %just take the logical or, since all we care about is having
                %the pixel on or off (regardless of it overlapping with
                %multiple bitmaps)
                M{ii}(:,sp(jj):sp(jj)+size(mm{jj},2)-1) = mm{jj} | ...
                      M{ii}(:,sp(jj):sp(jj)+size(mm{jj},2)-1);
            end
        end
    end
    max_width = max(max_width, size(M{ii},2));
end

%convert each image into one big array
%first pad the images to have the same number of columns, and optionally draw
%the segment line.
for ii=1:num_rows
    [h w] = size(M{ii});
    %zero pad missing columns from this row
    M{ii} = double([M{ii}, zeros(h, max_width - w)]);

    %convert the row to RGB (if neccessary)
    if display_segments || map_scores
        if map_scores && ((display_original && rem(ii,2) == 0) || ...
                          ~display_original)
            M{ii} = label2rgb(M{ii}, map_colormap, 'k');
        else
            M{ii} = label2rgb(M{ii}, 'white', 'k');
        end
    end
    if display_segments
        %add the segment line to the image
        if display_original && rem(ii,2) == 0
            sp = start_pos{ii/2};
        elseif ~ display_original
            sp = start_pos{ii};
        else
            sp = [];
        end
        for jj=1:length(sp)
            M{ii}(:,sp(jj),:) = repmat(segment_col, [h,1,1]);
        end
    end
end
M = cell2mat(M);

%now view the image
imshow(M);

%save the image to disk if required.
if save_results
    fprintf('writing ocr image to disk\n');
    imwrite(M, [img_prefix, '.', img_format], img_format);
end

fprintf('Elapsed time: %f\n', toc);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
