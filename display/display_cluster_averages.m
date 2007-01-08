function h = display_cluster_averages(Clust, varargin)
%  DISPLAY_CLUSTER_AVERAGES  Display clusters and their # of items as an image
%
%   h = display_cluster_averages(Clust, [VAR1, VAL1]...)
%
%   Clust should be a struct containing a cell array field labelled average, 
%   each of which is assumed to contain a matrix giving the average pixel 
%   intensity corresponding to the elements of that cluster.
%
%   h is the associated figure handle created when displaying the image
%
%   LOCAL VARS parameters specified below can have their defaults overriden by
%   specifiyng name and value pairs.  VAR1 should be a string specifying the
%   name of the variable to override.  VAL1 should be the new value for that
%   variable.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: display_cluster_averages.m,v 1.9 2007-01-08 22:05:49 scottl Exp $
%
% REVISION HISTORY
% $Log: display_cluster_averages.m,v $
% Revision 1.9  2007-01-08 22:05:49  scottl
% add ability to display reverse video (black text on whitebackground)
% image.  Also return the figure handle to the caller.
%
% Revision 1.8  2006-10-09 16:32:37  scottl
% change to allow optional parameter processing and variable overriding.
%
% Revision 1.7  2006-09-22 18:02:52  scottl
% added MSGID to warning message
%
% Revision 1.6  2006-09-18 21:01:47  scottl
% no change
%
% Revision 1.5  2006/09/05 15:50:44  scottl
% made use of MOCR_PATH variable for saving in the results directory.
%
% Revision 1.4  2006/08/07 21:19:11  scottl
% remove dependence on imview
%
% Revision 1.3  2006/07/21 20:17:26  scottl
% made textual character display more robust to small row/column widths.
%
% Revision 1.2  2006/07/05 01:05:34  scottl
% updated based on new Cluster structure.
%
% Revision 1.1  2006/06/03 20:55:54  scottl
% Initial check-in.
%
%

% LOCAL VARS %
%%%%%%%%%%%%%%
row_pix_border = 15;
col_pix_border = 5;
cl_h = 0;
cl_w = 0;

%how many clusters to display?  Inf = display all
num_clust = Inf;

%set save_averages to true to write the averages to disk based on the params
%below it
save_averages = false;
global MOCR_PATH;  %make use of the globally defined MOCR_PATH variable
img_prefix = [MOCR_PATH, '/results/cluster_averages'];
img_format = 'png';

%set this to true to reverse display (black text on white foreground)
reverse_display = false;


% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 1
    error('incorrect number of arguments specified!');
elseif nargin > 1
    process_optional_args(varargin{:});
end

if num_clust == Inf
    num_clust = Clust.num;
end

num_cols = ceil(sqrt(num_clust));
num_rows = ceil(num_clust / num_cols);

%first get the dimensions of the largest cluster, as well as the number of
%elements
cl_h = 0;
cl_w = 0;
num_comps = zeros(1,num_clust);
for ii=1:num_clust
    sz = size(Clust.avg{ii});
    cl_h = max(cl_h, sz(1));
    cl_w = max(cl_w, sz(2));
end

M = zeros(num_rows * (cl_h + row_pix_border), ...
          num_cols * (col_pix_border + cl_w + col_pix_border));

%get a sorted index of clusters in descending order based on their number of 
%components
[Dummy, sorted_clust_idx] = sort(Clust.num_comps,1,'descend');

row = 1;
col = 1;
X =[];
Y =[];
Txt = cell(1, num_clust);
txt_idx = 1;
half_h = ceil(cl_h/2);
half_w = ceil(cl_w/2);
for ii=sorted_clust_idx(1:num_clust)'
    %center the average of i at the current row and col position in M
    sz = size(Clust.avg{ii});
    tp = (row-1) * (cl_h + row_pix_border) + half_h - ceil(sz(1)/2) +1;
    bt = tp + sz(1) -1;
    lf = ((col-1) * (cl_w + col_pix_border)) + col * col_pix_border + half_w ...
         - ceil(sz(2)/2) +1;
    rg = lf + sz(2) -1;
    M(tp:bt,lf:rg) = Clust.avg{ii};

    %build the text area which will be overlaid on the image, showing the size 
    %of each cluster (number of elements)
    X = [X, ((col-1) * (cl_w + col_pix_border)) + col * col_pix_border + ...
             ceil(cl_w/2)];
    Y = [Y, ((row) * (cl_h + row_pix_border) - ceil(row_pix_border/2))];
    Txt{txt_idx} = num2str(Clust.num_comps(ii));
    txt_idx = txt_idx + 1;

    col = col+1;
    if col > num_cols
        row = row+1;
        col = 1;
    end
end
M = write_nums(M, X, Y, Txt);

if reverse_display
    M = 1-M;
end
h = imshow(M);

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
