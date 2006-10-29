function imgs = get_line_imgs(Lines, varargin)
% GET_LINE_IMGS  Use the Lines struct to return images of each line
%
%   IMGS = GET_LINE_IMGS(LINES, ['var1', val1]...)
%
%   The image of each line in the Lines struct is retrieved and returned as an
%   image.
%
%   LINES is a struct meant to hold all line related information, for the lines
%   associated with the componenets already discovered.  See GET_LINES() for
%   more info
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: get_line_imgs.m,v 1.1 2006-10-29 17:13:08 scottl Exp $
%
% REVISION HISTORY
% $Log: get_line_imgs.m,v $
% Revision 1.1  2006-10-29 17:13:08  scottl
% initial check-in
%


% LOCAL VARIABLES %
%%%%%%%%%%%%%%%%%%%
line_idx = [];  %indices of lines to return (leave empty for all lines)


% CODE START %
%%%%%%%%%%%%%%
if nargin < 1
    error('incorrect number of arguments specified!');
elseif nargin > 1
    process_optional_args(varargin{:});
end

if isempty(line_idx)
    line_idx = 1:Lines.num;
end
if size(line_idx,1) > 1
    line_idx = line_idx';
end

imgs = cell(length(line_idx),1);

pgs = unique(Lines.pg(line_idx))';
for pp=pgs
    M = ~imread(Lines.files{pp});
    idx = find(Lines.pg(line_idx) == pp);
    for jj=line_idx(idx)
        pos = Lines.pos(jj,:);
        imgs{jj} = M(pos(2):pos(4),pos(1):pos(3));
    end
end

