function display_lines(Lines, idx, varargin)
%  DISPLAY_LINES  Display the lines passed by index
%
%   display_lines(LINES, IDX, [VAR1, VAL1]...)
%
%   LINES should be a struct like that returned in get_Lines()
%
%   idx should be a vector listing the indices of the LINES struct to be drawn.
%
%   optional LOCAL VARS values below can be overriden specifying the name and
%   new value for the variable to be overwritten as additinoal parameters.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: display_lines.m,v 1.2 2006-11-07 02:52:41 scottl Exp $
%
% REVISION HISTORY
% $Log: display_lines.m,v $
% Revision 1.2  2006-11-07 02:52:41  scottl
% small bugfix for displaying x-height line.
%
% Revision 1.1  2006-09-22 17:56:06  scottl
% initial check-in
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%number of pixels of spacing between consecutive lines
row_spacing = 15;

%should we draw the baseline?
draw_baseline = true;
base_col = reshape([0,0,255],1,1,3);  % this is blue

%should we draw the x-height?
draw_xheight = true;
xheight_col = reshape([0,255,0],1,1,3);  %this is green

%set save_lines to true to write the line images to disk based on the params
%below it
save_lines = false;
global MOCR_PATH;  %make use of the globally defined MOCR_PATH variable
img_prefix = [MOCR_PATH, '/results/line_image'];
img_format = 'png';


% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 2
    error('incorrect number of arguments specified!');
elseif nargin > 2
    process_optional_args(varargin{:});
end

%first ensure the idx is a row vector (for if statement processing)
if size(idx,1) > 1
    idx = idx';
end

num_lines = length(idx);
M = cell(num_lines,1);
max_width = 0;

for pp = unique(Lines.pg(idx))';
    img = ~imread(Lines.files{pp});
    for ii = find(Lines.pg(idx) == pp)';
        pos = Lines.pos(idx(ii),:);
        M{ii} = img(pos(2):pos(4), pos(1):pos(3));
        max_width = max(max_width, size(M{ii},2));
    end
end

%pad each lines right-side with 0's to make them equal length, and add the
%border pixels
for ii=1:num_lines
    [h,w] = size(M{ii});
    M{ii} = [M{ii}, zeros(h, max_width - w); zeros(row_spacing, max_width)];
end

if draw_baseline || draw_xheight
    for ii=1:num_lines
        M{ii} = label2rgb(M{ii}, 'white', 'k');
        if draw_baseline
            M{ii}(1+Lines.baseline(idx(ii)),:,:)=repmat(base_col,1,max_width);
        end
        if draw_xheight
            M{ii}(1+Lines.xheight(idx(ii)),:,:)=repmat(xheight_col,1,max_width);
        end
    end
end

M = cell2mat(M);
imshow(M);

%save the image to disk if required.
if save_lines
    fprintf('writing line image to disk\n');
    imwrite(M, [img_prefix, '.', img_format], img_format);
end

fprintf('Elapsed time: %f\n', toc);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = write_nums(M, X, Y, Txt)
% helper function to convert strings of digits in Txt to images and place them
% centered about the corresponding X, and Y offsets in M which is an image
% matrix.
% This is meant to be a crude replacement for Matlab's 'text' function that 
% will be resized with the image, saved to disk, and work with imview
%
% Each digit image is sized to be about 9 pixels, so row_pix_border should be
% set larger than this to prevent overlap
%
% we also assume that all X, Y co-ordinates fall within the range of the number
% of rows and columns of M (no sanity checking done)
if length(X) ~= length(Y) || length(X) ~= length(Txt)
    error('X, Y, and Txt must all contain the same number of items');
end

%these represent each digit.  We assume 0 for background, and 1 for foreground
one = [0 0 0 1 1 0 0;
       0 0 1 0 1 0 0;
       0 1 0 0 1 0 0;
       0 0 0 0 1 0 0;
       0 0 0 0 1 0 0;
       0 0 0 0 1 0 0;
       0 0 0 0 1 0 0;
       0 0 0 0 1 0 0;
       0 1 1 1 1 1 1];

two = [0 1 1 1 1 1 0;
       1 0 0 0 0 0 1;
       0 0 0 0 0 0 1;
       0 0 0 0 0 1 0;
       0 0 0 1 1 0 0;
       0 0 1 1 0 0 0;
       0 0 1 0 0 0 0;
       0 1 0 0 0 0 0;
       0 1 1 1 1 1 1];

three = [0 0 1 1 1 1 0;
         0 1 0 0 0 0 1;
         0 0 0 0 0 0 1;
         0 0 0 0 0 0 1;
         0 0 0 1 1 1 0;
         0 0 0 0 0 0 1;
         0 0 0 0 0 0 1;
         0 1 0 0 0 0 1;
         0 0 1 1 1 1 0];

four = [0 0 0 0 0 1 0;
        0 0 0 1 1 1 0;
        0 0 1 0 0 1 0;
        0 1 0 0 0 1 0;
        1 0 0 0 0 1 0;
        1 1 1 1 1 1 1;
        0 0 0 0 0 1 0;
        0 0 0 0 0 1 0;
        0 0 0 0 0 1 0];

five = [0 1 1 1 1 1 1;
        0 1 0 0 0 0 0;
        0 1 0 0 0 0 0;
        0 1 1 1 1 1 0;
        0 0 0 0 0 0 1;
        0 0 0 0 0 0 1;
        0 0 0 0 0 0 1;
        0 1 0 0 0 0 1;
        0 0 1 1 1 1 0];

six =[0 0 0 1 1 0 0;
      0 0 1 0 0 0 0;
      0 1 0 0 0 0 0;
      0 1 1 1 1 0 0;
      0 1 0 0 0 1 0;
      0 1 0 0 0 1 1;
      0 1 0 0 0 1 1;
      0 1 0 0 0 1 0;
      0 0 1 1 1 0 0];

seven =[0 1 1 1 1 1 1;
        0 0 0 0 0 0 1;
        0 0 0 0 0 1 0;
        0 0 0 0 0 1 0;
        0 0 0 1 1 0 0;
        0 0 0 1 1 0 0;
        0 0 1 0 0 0 0;
        0 0 1 0 0 0 0;
        0 0 1 0 0 0 0];

eight =[0 0 1 1 1 1 0;
        0 1 0 0 0 0 1;
        0 1 0 0 0 0 1;
        0 0 1 1 1 1 0;
        0 1 0 0 0 0 1;
        0 1 0 0 0 0 1;
        0 1 0 0 0 0 1;
        0 1 0 0 0 0 1;
        0 0 1 1 1 1 0];

nine = [0 0 1 1 1 1 0;
        0 1 0 0 0 0 1;
        0 1 0 0 0 0 1;
        0 1 0 0 0 0 1;
        0 0 1 1 1 1 1;
        0 0 0 0 0 0 1;
        0 0 0 0 1 1 0;
        0 0 0 0 1 1 0;
        0 0 1 1 0 0 0];

zero = [0 0 1 1 1 1 0;
        0 1 0 0 0 0 1;
        0 1 0 0 1 1 1;
        0 1 0 1 0 0 1;
        0 1 0 1 0 0 1;
        0 1 1 0 0 0 1;
        0 1 0 0 0 0 1;
        0 1 0 0 0 0 1;
        0 0 1 1 1 1 0];

sz = size(M);
for ii=1:length(X)
    %convert the string into a large image matrix
    img = [];
    for cc=1:length(Txt{ii})
        switch Txt{ii}(cc)
        case '1'
            img = [img, one];
        case '2'
            img = [img, two];
        case '3'
            img = [img, three];
        case '4'
            img = [img, four];
        case '5'
            img = [img, five];
        case '6'
            img = [img, six];
        case '7'
            img = [img, seven];
        case '8'
            img = [img, eight];
        case '9'
            img = [img, nine];
        case '0'
            img = [img, zero];
        otherwise
            warning('MBOCR:charIgnored', ...
                    'this method only displays digits 0-9, character ignored');
        end
    end

    %overlay img centered about the X,Y position.  We may have to truncate/
    %overwrite other parts of M
    [h,w] = size(img);
    x = X(ii) - ceil(w/2) + 1;
    if x < 1
        diff = -x + 1;
        img = img(:,1+diff:end);
        w = w - diff;
        x = 1;
    end
    if x+w-1 > sz(2)
        diff = x+w-1 - sz(2);
        img = img(:,1:end-diff);
        w = w - diff;
    end
    y = Y(ii) - ceil(h/2) + 1;
    if y < 1
        diff = -y + 1;
        img = img(1+diff:end,:);
        h = h - diff;
        y = 1;
    end
    if y+h-1 > sz(1)
        diff = y+h-1 - sz(1);
        img = img(1:end-diff,:);
        h = h - diff;
    end
    M(y:y+h-1, x:x+w-1) = img;
end
