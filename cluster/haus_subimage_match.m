function [t_row, l_col, score] = haus_subimage_match(sub_img, dist_sub, img, ...
                                 dist_img, varargin)
% HAUS_SUBIMAGE_MATCH   Look for copies of the subimage inside the image passed.
%
%   [TOP_IDX, LEFT_IDX, SCORE] = HAUS_SUBIMAGE_MATCH(SUB_IMG, SUB_DIST, IMG, 
%                                IMG_DIST, [VAR1, VAL1]...)
%
%   SUB_IMG and IMG should be binary image matrices, with SUB_IMG smaller in 
%   size than IMG.
%
%   SUB_DIST and IMG_DIST should be corrsponding Euclidean distance transforms
%   of SUB_IMG and IMG above.  See bwdist() for further details
%
%   TOP_IDX, LEFT_IDX, are vectors containing the top and left pixel
%   co-ordinate offsets in IMG that SUB_IMG matches at.
%
%   SCORE is the score of the match.  It will be the maximal Euclidean distance
%   separation between the points of the sub-image to the image, and the image
%   to the sub-image (in the bouding box region surrounding the sub-image).
%
%   VAR1 and VAL1 are optional, and can be used to override the default values
%   for the LOCAL VARS defined below.  Each VAR1 should be a string giving the
%   exact name of the variable to override.  Each VAL1 should be the updated
%   value to use for that variable.
%   
%   If the 'thresh' parameter is not overridden, then it takes on a value of 1
%   by default
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: haus_subimage_match.m,v 1.1 2007-02-01 18:08:32 scottl Exp $
%
% REVISION HISTORY
% $Log: haus_subimage_match.m,v $
% Revision 1.1  2007-02-01 18:08:32  scottl
% initial revision
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%what is the maximal Euclidean distance between the pixels of the sub-image
%and the image before we cannot consider this a match?
thresh = 1;


% CODE START %
%%%%%%%%%%%%%%

%process input arguments
if nargin < 4
    error('must pass subimage, image and their distance transforms');
elseif nargin > 4
    process_optional_args(varargin{:});
end

sz = size(sub_img);
res = inf(size(img,1)-sz(1)+1, size(img,2)-sz(2)+1);

sub_on = find(sub_img == 1);

for r=1:size(res,1)
    fprintf('row: %d\r', r);
    for c=1:size(res,2)
        d_img = dist_img(r:r+sz(1)-1, c:c+sz(2)-1);
        img_on = find(img(r:r+sz(1)-1, c:c+sz(2)-1) == 1);
        if isempty(img_on)
            res(r,c) = inf;
        else
            res(r,c) = max(max(d_img(sub_on)), max(dist_sub(img_on)));
        end
    end
end
fprintf('\n');

[row_idx, col_idx] = find(res <= thresh);
idx = sub2ind(size(res), row_idx, col_idx);
score = res(idx);
t_row = row_idx; % - ceil(sz(1)/2) + 1;
l_col = col_idx; % - ceil(sz(2)/2) + 1;
