function [img,col1,row1,col2,row2] = crop_image(img, varargin)
% CROP_IMAGE   Remove extra pixels to create tight bounding box around image
%
%   [IMG, L_COL, T_ROW, R_COL, B_ROW] = CROP_IMAGE(IMG, [VAR1, VAL1]...)
%   IMG should be a matrix of pixel values (between 0 and 1), with 0 denoting
%       background
%
%   The positions L_COL, T_ROW, R_COL, and B_ROW returned refer to the left,
%   top, right, and bottom positions where the original image was cropped.
%
%   VAR1 and VAL1 are optional, and can be used to override the default values
%   for the LOCAL VARS defined below.  Each VAR1 should be a string giving the
%   exact name of the variable to override.  Each VAL1 should be the updated
%   value to use for that variable.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: crop_image.m,v 1.2 2006-12-17 19:59:56 scottl Exp $
%
% REVISION HISTORY
% $Log: crop_image.m,v $
% Revision 1.2  2006-12-17 19:59:56  scottl
% updated to include returning the updated row and column positions.
%
% Revision 1.1  2006-12-04 19:19:39  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%
%after what threshold is a pixel considered foreground
thresh = 0.5;

%by default create a tight bounding box.  Changing these leaves either the
%horizontal or vertical direction uncropped.
crop_vert = true;
crop_horiz = true;


% CODE START %
%%%%%%%%%%%%%%

%process input arguments
if nargin < 1
    error('input image must be passed as first argument');
elseif nargin > 1
    process_optional_args(varargin{:});
end

if crop_horiz
    [col1,col1] = find(img >= thresh, 1, 'first');
    [col2,col2] = find(img >= thresh, 1, 'last');
    img = img(:,col1:col2);
end

if crop_vert
    [row1,row1] = find(img' >= thresh, 1, 'first');
    [row2,row2] = find(img' >= thresh, 1, 'last');
    img = img(row1:row2,:);
end
