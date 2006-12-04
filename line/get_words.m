function words = get_words(line_img, space_width, varargin)
% GET_WORDS   Segment a line image into a cell array of word images.
%
%   WORDS = GET_WORDS(LINE_IMG, SPACE_WIDTH, [VAR1, VAL1]...)
%
%   LINE_IMG should be an array depicting an image to be chopped.
%
%   SPACE_WIDTH should be a positive scalar listing how many blank columns must
%   exist in the image for it to be considered for splitting.
%
%   VAR1 and VAL1 are optional, and can be used to override the default values
%   for the LOCAL VARS defined below.  Each VAR1 should be a string giving the
%   exact name of the variable to override.  Each VAL1 should be the updated
%   value to use for that variable.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: get_words.m,v 1.2 2006-12-04 23:27:55 scottl Exp $
%
% REVISION HISTORY
% $Log: get_words.m,v $
% Revision 1.2  2006-12-04 23:27:55  scottl
% small bugfix in creating the cell array
%
% Revision 1.1  2006-12-04 23:15:31  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%which pixel value denotes background?
bg_val = 0;


% CODE START %
%%%%%%%%%%%%%%

%process input arguments
if nargin < 2
    error('incorrect number of arguments');
elseif nargin > 2
    process_optional_args(varargin{:});
end

img_sum = sum(line_img);
bg_cols = find(img_sum == bg_val);
start_idx = [];
end_idx = [];
if bg_cols(1) ~= 1
    start_idx = 1;
end
while ~isempty(bg_cols)
    if length(bg_cols) <= space_width
        %can't have any more valid splits
        new_end = length(img_sum);
        while img_sum(new_end) == 0
            new_end = new_end - 1;
        end
        end_idx = [end_idx, new_end];
        bg_cols = [];
    elseif bg_cols(1)+space_width-1 == bg_cols(space_width)
        %valid space found
        end_idx = [end_idx, bg_cols(1) - 1];
        bg_cols = bg_cols(space_width:end);
        while length(bg_cols) > 1 && bg_cols(1) +1 == bg_cols(2)
            bg_cols = bg_cols(2:end);
        end
        start_idx = [start_idx, bg_cols(1) + 1];
        if length(bg_cols) == 1
            %last transition, we want to make sure we don't miss the last end
            %idx, so prefix a dummy bg_col (it will get removed below)
            bg_cols = [NaN, bg_cols];
        end
    end
    bg_cols = bg_cols(2:end);
end
if length(end_idx) > length(start_idx)
    %this can happen if we start the line at a space gap
    end_idx = end_idx(2:end);
end

words = cell(length(start_idx),1);
for ii=1:length(start_idx)
    words{ii} = line_img(:, start_idx(ii):end_idx(ii));
end
