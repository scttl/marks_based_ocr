function Pos = get_lines(Comps, num)
% get_lines  Determine bounding box of each "line" on a page based on neighbours
%
%   Pos = get_lines(Comps, [num])
%
%   The bounding box of each component in Comps as well as its directional 
%   neighbours is used to extract individual lines of a document (as determined 
%   by the left, top, right, and bottom pixel co-ordinates of the bounding box 
%   created and stored as a row in Pos).
%
%   Comps should be a struct containing neighbour, page, and position fields.
%   See cluster_comps for further details.
%
%   num is optional and if specified determines the maximum number of lines to
%   return (assuming there are that many lines in the files).  If not specified
%   all lines found are returned.

% CVS INFO %
%%%%%%%%%%%%
% $Id: get_lines.m,v 1.4 2006-07-05 00:51:40 scottl Exp $
%
% REVISION HISTORY
% $Log: get_lines.m,v $
% Revision 1.4  2006-07-05 00:51:40  scottl
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
%
%

% LOCAL VARIABLES %
%%%%%%%%%%%%%%%%%%%
max_num = inf;
num_pages = size(Comps.pg_size,1);
Pos = cell(num_pages,1);

display_images = false;
save_images = false;
%RGB colour values used to display images if set to true above
out_col = reshape([255,0,0],1,1,3);  % this is red

%the filename and file type to save the image (display_images must be true)
img_prefix = 'results/nips5_lines';
img_format = 'png';


% CODE START %
%%%%%%%%%%%%%%
if nargin < 1 || nargin > 2
    error('incorrect number of arguments specified!');
elseif nargin == 2
    max_num = num;
end

for pp=1:num_pages

    if max_num == 0
        %don't process further pages, and trim the empty cells
        Pos = Pos(1:pp-1);
        break;
    end

    fprintf('processing page %d\n', pp);

    %find top-left component 
    more_rows = false;
    pg_comps = find(Comps.pg == pp);
    if ~isempty(pg_comps)
        [v_val, vc]  = min(Comps.pos(pg_comps,2));
        vc = pg_comps(vc);
        more_rows = true;
    end
    first_comps = [];
    
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
        fprintf('new row found');
    
        while more_cols
            %should we expand the box dimensions?
            fprintf('.');
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
                    end
                end
            end
    
            %see if there are more items to the right in this row
            if Comps.nb(hc,3) ~= 0
                hc = Comps.nb(hc,3);
            else
                more_cols = false;
                fprintf('\n');
            end
        end
    
        %update the line boundaries
        Pos{pp} = [Pos{pp}; l, t, r, b];

        max_num = max_num - 1;
    end

    %draw the bounding boxes found
    if display_images || save_images
        M = label2rgb(~imread(Comps.files{pp}), 'white', 'k');
        for ii=1:size(Pos{pp},1)
            l = Pos{pp}(ii,1); t = Pos{pp}(ii,2); r = Pos{pp}(ii,3); 
            b = Pos{pp}(ii,4);
            M(t,l:r,:) = repmat(out_col,1,r-l+1);
            M(t:b,l,:) = repmat(out_col,b-t+1,1);
            M(t:b,r,:) = repmat(out_col,b-t+1,1);
            M(b,l:r,:) = repmat(out_col,1,r-l+1);
        end
        if display_images
            imview(M);
            fprintf('press any key to continue (to next page if it exists)\n');
            pause;
        end
        if save_images
            fprintf('saving page %d to disk\n', pp);
            imwrite(M, [img_prefix, num2str(pp), '.', img_format], img_format);
        end
        clear M;
    end
end


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the left-most cluster on this line
function lc = find_left(Comps, lc)
while Comps.nb(lc,1) ~= 0
    lc = Comps.nb(lc,1);
end


% get the right-most cluster on this line
function rc = find_right(Comps, rc)
while Comps.nb(rc,3) ~= 0
    rc = Comps.nb(rc,3);
end
