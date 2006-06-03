function display_cluster_averages(Clust, num)
%  DISPLAY_CLUSTER_AVERAGES  Display clusters and their # of items as an image
%
%   display_cluster_averages(Clust, num)
%
%   Clust should be an array of structs, each of which is assumed to contain 
%   a avg field which is a matrix giving the average pixel intensity
%   corresponding to the elements of that cluster.
%
%   num is optional and if specified, determines the number of cluster
%   averages to display.  If not specified, all clusters averages are shown.

% CVS INFO %
%%%%%%%%%%%%
% $Id: display_cluster_averages.m,v 1.1 2006-06-03 20:55:54 scottl Exp $
%
% REVISION HISTORY
% $Log: display_cluster_averages.m,v $
% Revision 1.1  2006-06-03 20:55:54  scottl
% Initial check-in.
%
%

% LOCAL VARS %
%%%%%%%%%%%%%%
row_pix_border = 15;
col_pix_border = 3;
cl_h = 0;
cl_w = 0;

%set save_averages to true to write the averages to disk based on the params
%below it
save_averages = false;
img_prefix = 'images/aa01_clust';
img_format = 'png';

% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 1 || nargin > 2
    error('incorrect number of arguments specified!');
elseif nargin == 2
    num_clust = num;
else
    num_clust = size(Clust,1);
end

num_cols = ceil(sqrt(num_clust));
num_rows = ceil(num_clust / num_cols);

%first get the dimensions of the largest cluster
for i=1:num_clust
    sz = size(Clust(i).avg);
    if sz(1) > cl_h
        cl_h = sz(1);
    end
    if sz(2) > cl_w
        cl_w = sz(2);
    end
end

M = zeros(num_rows * (cl_h + row_pix_border), ...
          num_cols * (col_pix_border + cl_w + col_pix_border));

row = 1;
col = 1;
X =[];
Y =[];
Txt = cell(1, num_clust);
for i=1:num_clust
    %center the average of i at the current row and col position in M
    sz = size(Clust(i).avg);
    tp = (row-1) * (cl_h + row_pix_border) + ceil(cl_h/2) - ceil(sz(1)/2) +1;
    bt = tp + sz(1) -1;
    lf = (col-1) * (cl_w + 2*col_pix_border) + ceil(cl_w/2) - ceil(sz(2)/2) +1;
    rg = lf + sz(2) -1;
    M(tp:bt,lf:rg) = Clust(i).avg;

    %build the text area which will be overlaid on the image, showing the size 
    %of each cluster (number of elements)
    X = [X, ((col-1) * (cl_w + 2 * col_pix_border) + ceil(cl_w/2))];
    Y = [Y, ((row) * (cl_h + row_pix_border) - ceil(row_pix_border/2))];
    Txt{i} = num2str(Clust(i).num);

    col = col+1;
    if col > num_cols
        row = row+1;
        col = 1;
    end
end
M = write_nums(M, X, Y, Txt);
imview(M);

%save the image to disk if required.
if save_averages
    fprintf('writing cluster averages to disk\n');
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

for i=1:length(X)
    %convert the string into a large image matrix
    img = [];
    for c=1:length(Txt{i})
        switch Txt{i}(c)
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
            warning('this method only displays digits 0-9, character ignored');
        end
    end

    %overlay img centered about the X,Y position
    [h,w] = size(img);
    x = X(i) - ceil(w/2) + 1;
    y = Y(i) - ceil(h/2) + 1;
    M(y:y+h-1, x:x+w-1) = img;
end

