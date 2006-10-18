function [Lines, Comps] = get_lines(Comps, varargin)
% get_lines  Determine each "line" on a page based on neighbours
%
%   [LINES, COMPS] = get_lines(COMPS, ['var1', val1]...)
%
%   The bounding box of each component in COMPS as well as its directional 
%   neighbours is used to extract individual lines of a document (as determined 
%   by the left, top, right, and bottom pixel co-ordinates of the bounding box.
%
%   COMPS should be a struct containing neighbour, page, and position fields.
%   See GET_COMPS for further details.  Note that the COMPS struct will have
%   its found_lines, modal_height, line, scale_factor, descender_off, and 
%   ascender_off fields updated, and this updated struct is returned as the 
%   second output argument
%
%   Optional parameters should be name, value pairs that override the default
%   values for the local variables specified in the code below.
%
%   LINES is a struct meant to hold all line related information, for the lines
%   associated with the componenets already discovered.  It contains the
%   following fields:
%     num - An integer specifying the number of lines found
%     pg - A num length vector specifying upon which page this line belongs.
%     files - A cell array containing the string path to each image file that
%             contains lines.
%     pos - A numx4 matrix specifying the pixel co-ordinates of the left, top,
%           right, and bottom boundaries of the line (the boudning box)
%     baseline - A vector listing how many pixels down from the top line region
%                boundary each baseline is (the baseline is the position at
%                which most characters have their bottom region (exceptions are
%                characters like j, g, p, etc.
%     xheight - A vector listing how many pixels down from the top line region
%               boundary each x-height is (the x-height is the position at
%               which most characters have their top region -- it is the height
%               of lowercase characters like x, e, a, etc.)


% CVS INFO %
%%%%%%%%%%%%
% $Id: get_lines.m,v 1.7 2006-10-18 16:01:20 scottl Exp $
%
% REVISION HISTORY
% $Log: get_lines.m,v $
% Revision 1.7  2006-10-18 16:01:20  scottl
% by default use the full height (instead of x-height) as scale factor
% since this improves results.
%
% Revision 1.6  2006-10-09 16:27:52  scottl
% rewritten, introducing a Lines struct, etc.
%
% Revision 1.5  2006-08-14 01:15:48  scottl
% remove dependence on imview
%
% Revision 1.4  2006/07/05 00:51:40  scottl
% Re-written based on Component and Cluster structure refinements.
%
% Revision 1.3  2006/06/19 21:52:12  scottl
% small change to use length instead of size.
%
% Revision 1.2  2006/06/12 20:59:40  scottl
% use cell array for comps, allow one to specify the number of lines returned
%
% Revision 1.1  2006/06/03 20:56:01  scottl
% Initial check-in.



% LOCAL VARIABLES %
%%%%%%%%%%%%%%%%%%%
max_num = inf;  %maximum # of lines to find.

display_images = false;
save_images = false;
%RGB colour values used to display images if set to true above
out_col = reshape([255,0,0],1,1,3);  % this is red

%the filename and file type to save the image (save_images must be true)
global MOCR_PATH;
img_prefix = [MOCR_PATH, '/results/lines_pg'];
img_format = 'png';

%what percentage of pixels in a row must be 'on' for consideration as a
%baseline
base_thresh = 0.2;

%what percentage of pixels in a row must be 'on' for consideration as an
%x-height line.
xheight_thresh = 0.3;

%by defeault turn-off the override warning display since it will be generated
%for each line we create.
override_display = 'OFF';   %other legal value is 'ON'

%print each line and column as it is processed?
print_lines = false;



% CODE START %
%%%%%%%%%%%%%%
tic;
if nargin < 1
    error('incorrect number of arguments specified!');
elseif nargin > 1
    process_optional_args(varargin{:});
end
warning(override_display, 'MBOCR:override');

num_pages = size(Comps.pg_size,1);

%initialize the Lines struct as well as the new Comps fields
Lines = init_lines(Comps.files);
Comps.found_lines = true;
Comps.line = zeros(Comps.max_comp, 1, 'uint64');
Comps.descender_off = zeros(Comps.max_comp, 1, 'int16');
Comps.ascender_off = zeros(Comps.max_comp, 1, 'int16');
Comps.modal_height = NaN;
Comps.scale_factor = ones(Comps.max_comp, 1);

for pp=1:num_pages

    if max_num == 0
        break;
    end

    fprintf('%.2fs: processing page %d\n', toc, pp);

    %find top-left component on this page
    more_rows = false;
    pg_comps = find(Comps.pg == pp);
    if ~isempty(pg_comps)
        [vc, vc]  = min(Comps.pos(pg_comps,2));
        vc = pg_comps(vc);
        more_rows = true;
    end
    first_comps = [];

    %open the page image (for determining line statistics like baseline, etc.)
    M = ~imread(Comps.files{pp});
    
    %follow neighbours to get bounding box
    while more_rows && max_num > 0
    
        %there may be left neighbours that aren't as high up the page.
        vc = find_left(Comps, vc);

        pos = Comps.pos(vc,:);
        l = pos(1); t = pos(2); r = pos(3); b = pos(4);
        hc = vc;
        more_cols = true;
        more_rows = false;
        next_top = Inf;
        first_comps = [first_comps, hc];
        if print_lines fprintf('->'); end

        Lines.num = Lines.num + 1;
        Comps.line(hc) = Lines.num;
    
        while more_cols
            %should we expand the box dimensions?
            if print_lines fprintf('.'); end
            pos = Comps.pos(hc,:);
            t = min(t, pos(2));
            r = max(r, pos(3));
            b = max(b, pos(4));
    
            %have we found a new closer next row?
            new_vc = Comps.nb(hc,4);
            if new_vc ~= 0
                if Comps.pos(new_vc,2) < next_top
                    new_firstvc = find_left(Comps, new_vc);
                    if all(new_firstvc ~= first_comps)
                        %the first component's bottom neighbour could be on the 
                        %same line, so we scan across to ensure that they have 
                        %different far right neighbours to ensure we don't 
                        %repeat the same line
                        if hc ~= first_comps(end) || ...
                           any(find_right(Comps,hc) ~= find_right(Comps,new_vc))
                            more_rows = true;
                            next_top = Comps.pos(new_vc,2);
                            vc = new_vc;
                        end
                    else
                        %since on the same line, see if we should expand the
                        %right, or bottom dimensions based on those of the
                        %bottom neighbour
                        r = max(r, Comps.pos(new_vc,3));
                        b = max(b, Comps.pos(new_vc,4));
                        Comps.line(new_vc) = Lines.num;
                    end
                end
            end

            %is there a lower top neighbour in the same row?
            new_vc = Comps.nb(hc,2);
            if new_vc ~= 0
                if Comps.pos(new_vc,2) < t
                    newfirstvc = find_left(Comps, new_vc);
                    if newfirstvc == first_comps(end)
                        t = Comps.pos(new_vc,2);
                        Comps.line(new_vc) = Lines.num;
                    end
                end
            end
    
            %see if there are more items to the right in this row
            if Comps.nb(hc,3) ~= 0
                hc = Comps.nb(hc,3);
                Comps.line(hc) = Lines.num;
            else
                more_cols = false;
                if print_lines fprintf('\n'); end
            end
        end

        %since some components may have been missed during line boundary 
        %detection above, we must add their associated line too
        missing = Comps.pg==pp & Comps.line==0 & Comps.pos(:,1)>=l & ...
                  Comps.pos(:,2)>=t & Comps.pos(:,3)<=r & Comps.pos(:,4)<=b;
        Comps.line(missing) = Lines.num;
    
        %set the remaining fields
        Lines.pg(Lines.num,:) = pp;
        Lines.pos(Lines.num,:) = [l, t, r, b];

        %get the baseline and x-height
        [Lines.baseline(Lines.num,:), Lines.xheight(Lines.num,:)] = ...
        get_line_props(M(t:b,l:r), 'base_thresh', base_thresh, ...
                       'xheight_thresh', xheight_thresh);

        max_num = max_num - 1;
    end

    %components that at least partially overlap line boundary edges may still
    %be unassociated with a line, so we must add them to an appropriate line
    %and extend boundaries as required.
    for ii = find(Comps.pg == pp & Comps.line == 0)'
        nbs = Comps.nb(ii,:);
        if nbs(1) == 0 && nbs(3) == 0
            fprintf('new line for comp %d\n', ii);
            %create a new line for this component
            Lines.num = Lines.num + 1;
            match_line = Lines.num;
            Lines.pg(Lines.num) = pp;
            Lines.pos(Lines.num,:) = [Inf, Inf, -Inf, -Inf];
        else
            %try and add to the same line as its left or right neighbour
            if nbs(1) == 0 || (nbs(1) ~= 0 && nbs(3) ~= 0 && ...
                Comps.nb_dist(ii,1) > Comps.nb_dist(ii,3))
                first_dir = 3;  %try right first
                second_dir = 1;
            else
                first_dir = 1;  %try left first
                second_dir = 3;
            end
            match_line = Comps.line(nbs(first_dir));
            while match_line == 0 && Comps.nb(nbs(first_dir),first_dir) ~= 0;
                nbs(first_dir) = Comps.nb(nbs(first_dir),first_dir);
                match_line = Comps.line(nbs(first_dir));
            end
            if match_line == 0
                %try the other direction
                while match_line == 0 && nbs(second_dir) ~= 0;
                    match_line = Comps.line(nbs(second_dir));
                    nbs(second_dir) = Comps.nb(nbs(second_dir),second_dir);
                end
                if match_line == 0
                    %can't find any left or right neighbour belonging to a
                    %line.  Create a new one.
                    Lines.num = Lines.num + 1;
                    match_line = Lines.num;
                    Lines.pg(Lines.num) = pp;
                    Lines.pos(Lines.num,:) = [Inf, Inf, -Inf, -Inf];
                end
            end
        end
        Comps.line(ii) = match_line;
        l = min(Lines.pos(match_line,1), Comps.pos(ii,1));
        t = min(Lines.pos(match_line,2), Comps.pos(ii,2));
        r = max(Lines.pos(match_line,3), Comps.pos(ii,3));
        b = max(Lines.pos(match_line,4), Comps.pos(ii,4));
        Lines.pos(match_line,:) = [l,t,r,b];
        %top or bottom change could have an impact on baselines and x-height
        %so recalculate these
        [Lines.baseline(match_line,:), Lines.xheight(match_line,:)] = ...
        get_line_props(M(t:b,l:r), 'base_thresh', base_thresh, ...
                       'xheight_thresh', xheight_thresh);
    end

    %draw the bounding boxes found
    if display_images || save_images
        M = label2rgb(M, 'white', 'k');
        for ii=find(Lines.pg == pp)';
            l = Lines.pos(ii,1); t = Lines.pos(ii,2); r = Lines.pos(ii,3); 
            b = Lines.pos(ii,4);
            M(t,l:r,:) = repmat(out_col,1,r-l+1);
            M(t:b,l,:) = repmat(out_col,b-t+1,1);
            M(t:b,r,:) = repmat(out_col,b-t+1,1);
            M(b,l:r,:) = repmat(out_col,1,r-l+1);
        end
        if display_images
            imshow(M);
            fprintf('press any key to continue (to next page if it exists)\n');
            pause;
        end
        if save_images
            fprintf('%.2fs: saving page %d to disk\n', toc, pp);
            imwrite(M, [img_prefix, num2str(pp), '.', img_format], img_format);
        end
    end
end

fprintf('%.2fs: calculating Component ascender and descender offsets\n', toc);
line_tops = Lines.pos(Comps.line, 2);
line_bases = Lines.baseline(Comps.line);
line_xheights = Lines.xheight(Comps.line);
Comps.descender_off(:) = Comps.pos(:,4) - (line_tops + line_bases);
Comps.ascender_off(:) = (line_tops + line_xheights) - Comps.pos(:,2);

fprintf('%.2fs: calculating Component scale factors\n', toc);
%comp_heights = double(line_bases - line_xheights + 1);
comp_heights = double(Comps.pos(:,4) - Comps.pos(:,2) + 1);
Comps.modal_height = mode(comp_heights);
%Comps.modal_height = max(comp_heights);
Comps.scale_factor = Comps.modal_height ./ comp_heights;

%if warnings have been turned off, turn them back on
if override_display == 'OFF'
    warning('ON', 'MBOCR:override');
end

fprintf('%.2fs: TOTAL NUMBER OF LINES FOUND: %d\n', toc, Lines.num);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this function creates an empty line structure
function Lines = init_lines(files)
Lines.num = 0;
Lines.pg = uint32([]);
Lines.files = files;
Lines.pos = uint16([]);
Lines.baseline = uint16([]);  %must be same class as pos, so +,- work
Lines.xheight = uint16([]);


% get the left-most component in this line
function lc = find_left(Comps, lc)
while Comps.nb(lc,1) ~= 0
    lc = Comps.nb(lc,1);
end


% get the right-most component in this line
function rc = find_right(Comps, rc)
while Comps.nb(rc,3) ~= 0
    rc = Comps.nb(rc,3);
end
