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
% $Id: display_lines.m,v 1.3 2007-01-25 18:42:10 scottl Exp $
%
% REVISION HISTORY
% $Log: display_lines.m,v $
% Revision 1.3  2007-01-25 18:42:10  scottl
% added ability to draw component bounding boxes, reverse the output colours,
% and removed an unneeded subfunction.
%
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

%should we draw component bounding boxes? (this requires Comps to be passed)
Comps = [];
draw_comps = false;
exclude_comps = [];  %input component numbers to exclude.  Like spaces
comps_col = reshape([255,0,0],1,1,3);  %this is red

%set save_lines to true to write the line images to disk based on the params
%below it
save_lines = false;
global MOCR_PATH;  %make use of the globally defined MOCR_PATH variable
img_prefix = [MOCR_PATH, '/results/line_image'];
img_format = 'png';

%should we reverse the display?  black ink on white background
reverse_display = false;


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
    if reverse_display
        M{ii} = 1 - M{ii};
    end
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
        if draw_comps && ~isempty(Comps)
            comp_idcs = find(Comps.line == idx(ii));
            if ~isempty(exclude_comps)
                comp_idcs = setdiff(comp_idcs, exclude_comps);
            end
            line_l = Lines.pos(idx(ii),1);
            line_t = Lines.pos(idx(ii),2);
            for jj=comp_idcs'
                l = Comps.pos(jj,1); t = Comps.pos(jj,2);
                r = Comps.pos(jj,3); b = Comps.pos(jj,4);
                l = l-line_l+1; t = t-line_t+1;
                r = r-line_l+1; b = b-line_t+1;
    
                M{ii}(t,l:r,:) = repmat(comps_col,1,r-l+1);
                M{ii}(t:b,l,:) = repmat(comps_col,b-t+1,1);
                M{ii}(t:b,r,:) = repmat(comps_col,b-t+1,1);
                M{ii}(b,l:r,:) = repmat(comps_col,1,r-l+1);
            end
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
