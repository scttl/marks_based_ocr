function display_ocr_results(indices, start_pos, bitmaps, base_offs, num)
%  DISPLAY_OCR_RESULTS  Construct an image of the OCR'd indices passed
%
%   display_ocr_results(indices, bitmaps, base_offs, [num])
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

% CVS INFO %
%%%%%%%%%%%%
% $Id: display_ocr_results.m,v 1.4 2006-08-05 17:35:29 scottl Exp $
%
% REVISION HISTORY
% $Log: display_ocr_results.m,v $
% Revision 1.4  2006-08-05 17:35:29  scottl
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

%set save_averages to true to write the averages to disk based on the params
%below it
save_averages = false;
img_prefix = 'results/nips5_line1to50_ocr';
img_format = 'png';

display_segments = true;  %draw the segement lines in a different colour?
segment_col = reshape([255,0,0],1,1,3);  %this is red

row_margin = 10;  %number of pixels between consecutive rows in the image


% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 4 || nargin > 5
    error('incorrect number of arguments specified!');
elseif nargin == 5
    num_rows = num;
    if num_rows > size(indices,1)
        num_rows = size(indices,1);
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

M = cell(num_rows,1);

if ~ iscell(indices)
    indices = {indices};
end



max_width = 0;
for ii=1:num_rows
    mm = bitmaps(indices{ii});
    moffs = base_offs(indices{ii});
    max_height = 0;
    for jj=1:length(mm)
        max_height = max(max_height, size(mm{jj},1)-moffs(jj));
    end
    max_base = max(moffs);
    for jj=1:length(mm)
        [h w] = size(mm{jj});
        mm{jj} = [zeros(max_height - (h-moffs(jj)), w); ...
                 mm{jj}; ...
                 zeros(max_base - moffs(jj), w)];
    end
    M{ii} = [];
    for jj=1:length(mm)
        w = size(M{ii},2);
        extra_width = start_pos(jj) + size(mm{jj},2) - 1 - w;
        if extra_width > 0
            M{ii} = [M{ii}, zeros(max_height+max_base,extra_width)];
        end
        M{ii}(:,start_pos(jj):start_pos(jj)+size(mm{jj},2)-1) = mm{jj} | ...
              M{ii}(:,start_pos(jj):start_pos(jj)+size(mm{jj},2)-1);
    end
    max_width = max(max_width, size(M{ii},2));
end

%convert each image into one big array
%first pad the images to have the same number of columns, and optionally draw
%the segment line.
for ii=1:num_rows
    [h w] = size(M{ii});
    M{ii} = [M{ii}, zeros(h, max_width - w)];
    if display_segments
        %add the segment line to the image (after converting to RGB)
        M{ii} = label2rgb(M{ii}, 'white', 'k');
        for jj=1:length(start_pos)
            M{ii}(:,start_pos(jj),:) = repmat(segment_col, ...
                                       [max_height+max_base,1,1]);
        end
    end
end
M = cell2mat(M);

%now view the image
imshow(M);

%save the image to disk if required.
if save_averages
    fprintf('writing ocr image to disk\n');
    imwrite(M, [img_prefix, '.', img_format], img_format);
end

fprintf('Elapsed time: %f\n', toc);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
