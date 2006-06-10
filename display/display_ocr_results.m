function display_ocr_results(indices, bitmaps, num)
%  DISPLAY_OCR_RESULTS  Construct an image of the OCR'd indices passed
%
%   display_ocr_results(indices, bitmaps, [num])
%
%   indices should either be a vector of index values into the bitmaps cell 
%   array, or it should be a cell array of index vectors (one per row).  Each
%   entry is assumed to lie within the length of the bitmaps array.  Multiple
%   vectors are drawn in order (from top-to-bottom).
%
%   bitmaps should be a cell array of logical arrays, each of which represents
%   a single character in the alphabet, which are used to construct the output
%   image.
%
%   num is optional and if specified, determines the number of OCR'd
%   lines to display.  If not specified, all indices are drawn.

% CVS INFO %
%%%%%%%%%%%%
% $Id: display_ocr_results.m,v 1.1 2006-06-10 21:01:32 scottl Exp $
%
% REVISION HISTORY
% $Log: display_ocr_results.m,v $
% Revision 1.1  2006-06-10 21:01:32  scottl
% Initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%set save_averages to true to write the averages to disk based on the params
%below it
save_averages = false;
img_prefix = 'results/aa01_ocr';
img_format = 'png';

row_margin = 10;  %number of pixels between consecutive rows in the image

% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 2 || nargin > 3
    error('incorrect number of arguments specified!');
elseif nargin == 3
    num_rows = num;
    if num_rows > size(indices,1)
        num_rows = size(indices,1);
    end
else
    num_rows = size(indices,1);
end

M = cell(num_rows,1);

if ~ iscell(indices)
    tmp = cell(1);
    tmp{1} = indices;
    indices = tmp;
end

max_width = 0;
for i=1:num_rows
    M{i} = cell2mat(bitmaps(indices{i}));
    max_width = max(max_width, size(M{i},2));
end

%convert each image into one big array
%first pad the images to have the same number of columns
for i=1:num_rows
    [h w] = size(M{i});
    M{i} = [M{i}, zeros(h, max_width - w)];
end
M = cell2mat(M);

%now view the image
imview(M);

%save the image to disk if required.
if save_averages
    fprintf('writing ocr image to disk\n');
    imwrite(M, [img_prefix, '.', img_format], img_format);
end

fprintf('Elapsed time: %f\n', toc);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
