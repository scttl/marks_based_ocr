function img = imgcell2mat(img, varargin)
% IMGCELL2MAT   Converts a cell array of images into a padded single image
%
%   IMG = IMGCELL2MAT(IMG, [VAR1, VAL1]...)
%
%   IMG is a cell vector of matrices each of which represents an image.  The
%   output concatenates the elements in the cell vector, padding them so that 
%   the built-in cell2mat can be called without error.
%
%   VAR1 and VAL1 are optional, and can be used to override the default values
%   for the LOCAL VARS defined below.  Each VAR1 should be a string giving the
%   exact name of the variable to override.  Each VAL1 should be the updated
%   value to use for that variable.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: imgcell2mat.m,v 1.1 2006-12-02 16:20:34 scottl Exp $
%
% REVISION HISTORY
% $Log: imgcell2mat.m,v $
% Revision 1.1  2006-12-02 16:20:34  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%how much space to add between rows if creating a column vector image (or cols
%if creating a row vector image)
elem_padding = 4;

%which value to use for padding
pad_val = 0;


% CODE START %
%%%%%%%%%%%%%%
%process input arguments
if nargin < 1
    error('must pass in an image');
elseif nargin > 1
    process_optional_args(varargin{:});
end

%arg checking
if ~iscell(img)
    error('input image must be a cell vector');
elseif length(img) ~= prod(size(img))
    error('input image array must be a vector (not a matrix)');
end

num = length(img);
max_len = 0;
sz = zeros(num,2);
row_major = true;
if size(img,2) > 1
    row_major = false;
end

for ii=1:num
    sz(ii,:) = size(img{ii});
    if row_major
        max_len = max(sz(ii,2), max_len);
    else
        max_len = max(sz(ii,1), max_len);
    end
end

for ii=1:num
    if row_major
        img{ii} = [img{ii}, pad_val + zeros(sz(ii,1),max_len-sz(ii,2))];
        if ii < num && elem_padding > 0
            img{ii} = [img{ii}; pad_val + zeros(elem_padding,max_len)];
        end
    else
        img{ii} = [img{ii}; pad_val + zeros(max_len-sz(ii,1),sz(ii,2))];
        if ii < num && elem_padding > 0
            img{ii} = [img{ii}, pad_val + zeros(max_len,elem_padding)];
        end
    end
end
img = cell2mat(img);
