function [t_row, l_col, score] = subimage_match(sub_img, img, varargin)
% SUBIMAGE_MATCH   Look for copies of the subimage inside the image passed.
%
%   [TOP_IDX, LEFT_IDX, SCORE] = SUBIMAGE_MATCH(SUB_IMG, IMG, [VAR1, VAL1]...)
%
%   SUB_IMG and IMG should be image matrices, with SUB_IMG smaller in size
%   than IMG.  images are assumed to contain values from 0 to 1, with 0
%   representing background pixels.
%
%   TOP_IDX, LEFT_IDX, are vectors containing the top and left pixel
%   co-ordinate offsets in IMG that SUB_IMG matches at.
%
%   SCORE is the score of the match.  It will be the sum of the values of the
%   pixels that mismatch based on the bounding box of SUB_IMG (this will be the
%   Hamming distance for binary images)
%
%   VAR1 and VAL1 are optional, and can be used to override the default values
%   for the LOCAL VARS defined below.  Each VAR1 should be a string giving the
%   exact name of the variable to override.  Each VAL1 should be the updated
%   value to use for that variable.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: subimage_match.m,v 1.1 2006-12-04 23:14:06 scottl Exp $
%
% REVISION HISTORY
% $Log: subimage_match.m,v $
% Revision 1.1  2006-12-04 23:14:06  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%what proportion of the total pixel density can mismatch to still be considered
%a valid match?
thresh = 0;

% CODE START %
%%%%%%%%%%%%%%

%process input arguments
if nargin < 2
    error('must pass a subimage and the image to match against as parameters');
elseif nargin > 2
    process_optional_args(varargin{:});
end

sz = size(sub_img);
area = prod(sz);
sub_sum = sum(sub_img(:));
inv_sub = 1 - sub_img;
inv_sub_sum = sum(inv_sub(:));

pos_corr = abs(filter2(sub_img,img) - sub_sum);
neg_corr = filter2(inv_sub,img);

[row_idx, col_idx] = find(pos_corr + neg_corr <= thresh * area);
idx = sub2ind(size(pos_corr), row_idx, col_idx);
score = (pos_corr(idx) + neg_corr(idx)) / area;
t_row = row_idx - ceil(sz(1)/2) + 1;
l_col = col_idx - ceil(sz(2)/2) + 1;
