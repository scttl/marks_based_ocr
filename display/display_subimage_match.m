function display_subimage_match(sub_img, img, t_row, l_col, varargin)
% DISPLAY_SUBIMAGE_MATCH   Display where the subimage matches the image
%
%   DISPLAY_SUBIMAGE_MATCH(SUB_IMG, IMG, ROWS, COLS, [VAR1, VAL1]...)
%
%   SUB_IMG and IMG should be image matrices that use the value 0 for backgound
%   and must have foreground pixels with a value no larger than 1.
%
%   ROWS, COLS are the top and left corner offsets in IMG that the SUB_IMG
%   matches at.  These are typically calculated in subimage_match()
%
%   VAR1 and VAL1 are optional, and can be used to override the default values
%   for the LOCAL VARS defined below.  Each VAR1 should be a string giving the
%   exact name of the variable to override.  Each VAL1 should be the updated
%   value to use for that variable.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: display_subimage_match.m,v 1.2 2006-12-17 19:50:26 scottl Exp $
%
% REVISION HISTORY
% $Log: display_subimage_match.m,v $
% Revision 1.2  2006-12-17 19:50:26  scottl
% bugfix in processing optional arguments, use imtool instead of imshow
%
% Revision 1.1  2006-12-04 23:13:35  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%list the rgb colours to display the images in.  The first should be the
%foreground pixel colour, and the second the matching pixel colour.
map = [1 1 1; 1 0 0];  %use white for the non-matches and red for the matches

%the background colour (see label2rgb)
bg_col = 'k';

%set save_matches to true to write the line images to disk based on the params
%below it
save_matches = false;
global MOCR_PATH;  %make use of the globally defined MOCR_PATH variable
img_prefix = [MOCR_PATH, '/results/match_image'];
img_format = 'png';


% CODE START %
%%%%%%%%%%%%%%

if nargin < 4
    error('incorrect number of arguments specified');
elseif nargin > 4
    process_optional_args(varargin{:});
end

sz = size(sub_img);
img_sz = size(img);

img = double(img);
img(img ~= 0) = 1;  %need to binarize the image (for labelling to work)
for ii=1:length(t_row)
    sr = t_row(ii);
    er = t_row(ii)+sz(1)-1;
    sc = l_col(ii);
    ec = l_col(ii)+sz(2)-1;
    if sr < 1
        er = er + sr - 1;
        sr = 1;
    end
    if er > img_sz(1)
        er = img_sz(1);
    end
    if sc < 1
        ec = ec + sc - 1;
        sc = 1;
    end
    if ec > img_sz(2)
        ec = img_sz(2);
    end
    xx = img(sr:er, sc:ec);
    xx(xx ~= 0) = 2;
    img(sr:er, sc:ec) = xx;
end
img = label2rgb(img, map, bg_col);
imtool(img);

if save_matches
    fprintf('writing image to disk\n');
    imwrite(img, [img_prefix, '.', img_format], img_format);
end
