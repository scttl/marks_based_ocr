function display_offset_bitmaps(bitmaps, base_offs)
%  DISPLAY_OFFSET_BITMAPS  Construct an image of the aligned character bitmaps
%
%   display_offset_bitmaps(bitmaps, base_offs)
%
%   bitmaps should be a cell array of logical arrays, each of which represents
%   a single character in the alphabet, which are used to construct the output
%   image.
%
%   base_offs should be a vector listing the offset relative to the baseline of
%   the bottom of the corresponding bitmap.  This ensure characters are lined
%   up correctly.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: display_offset_bitmaps.m,v 1.2 2006-08-24 21:44:06 scottl Exp $
%
% REVISION HISTORY
% $Log: display_offset_bitmaps.m,v $
% Revision 1.2  2006-08-24 21:44:06  scottl
% remove dependence on imview
%
% Revision 1.1  2006/06/19 20:59:38  scottl
% Initial check-in.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%set save_bitmaps to true to write the bitmap image to disk based on the params
%below it
save_bitmaps = false;
img_prefix = 'results/nips5_bitmaps';
img_format = 'png';

img_spacing = 5;  %number of pixels between bitmaps

display_baseline = true;  %if this is draw we draw the baseline overtop the 
                          %bitmaps
baseline_col = reshape([255,0,0],1,1,3);  % this is red


% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin ~= 2
    error('incorrect number of arguments specified!');
end

num_bitmaps = length(bitmaps);
if num_bitmaps ~= length(base_offs)
    error('number of bitmaps must be the same as the number of base_offs!');
end


max_height = 0;
max_base = max(base_offs);
M = [];

for i=1:num_bitmaps
    max_height = max(max_height, size(bitmaps{i},1)-base_offs(i));
end

for i=1:num_bitmaps
    [h, w] = size(bitmaps{i});
    X = [zeros(max_height - (h-base_offs(i)), w); ...
         bitmaps{i}; ...
         zeros(max_base - base_offs(i), w)];
    M = [M, X, zeros(max_height+max_base, img_spacing)];
end

if display_baseline
    %add the baseline to M (after converting to RGB)
    M = label2rgb(M, 'white', 'k');
    num_cols = size(M,2);
    M(max_height+1,:,:) = repmat(baseline_col, [1,num_cols,1]);
end

%now view the image
imshow(M);

%save the image to disk if required.
if save_bitmaps
    fprintf('writing bitmaps image to disk\n');
    imwrite(M, [img_prefix, '.', img_format], img_format);
end

fprintf('Elapsed time: %f\n', toc);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
