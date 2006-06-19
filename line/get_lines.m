function Pos = get_lines(Clust, Comps, num)
% get_lines  Determine bounding box of each "line" on a page based on neighbours
%
%   Pos = get_lines(Clust, Comps, [num])
%
%   The bounding box of each component in Clust as well as its directional 
%   neighbours is used to extract individual lines of a document (as determined 
%   by the left, top, right, and bottom pixel co-ordinates of the bounding box 
%   created and stored as a row in Pos).
%
%   Clust should be an array of structs, each of which is assumed to contain 
%   several fields of a particular format.   See cluster_comps for details.
%
%   Comps should be a 2 or 3 dimensional image matrix, where each entry 
%   represents a pixel, and each 'on' pixel is labelled with the component
%   number to which it belongs.
%
%   num is optional and if specified determines the maximum number of lines to
%   return (assuming there are that many lines in the files).  If not specified
%   all lines found are returned.

% CVS INFO %
%%%%%%%%%%%%
% $Id: get_lines.m,v 1.3 2006-06-19 21:52:12 scottl Exp $
%
% REVISION HISTORY
% $Log: get_lines.m,v $
% Revision 1.3  2006-06-19 21:52:12  scottl
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
num_pages = length(Comps);
Pos = cell(num_pages,1);
num_clust = length(Clust);

display_images = false;
save_images = false;
%RGB colour values used to display images if set to true above
out_col = reshape([255,0,0],1,1,3);  % this is red

%the filename and file type to save the image (display_images must be true)
img_prefix = 'results/nips5_lines';
img_format = 'png';


% CODE START %
%%%%%%%%%%%%%%
if nargin < 2 || nargin > 3
    error('incorrect number of arguments specified!');
elseif nargin == 3
    max_num = num;
end

for p=1:num_pages

    fprintf('processing page %d\n', p);

    %find top-left component 
    more_rows = false;
    min_top  = Inf;
    vcl = 0;
    for i=1:num_clust
        idx = find(Clust(i).pg == p);
        if ~isempty(idx)
            [val, off] = min(Clust(i).pos(idx,2));
            if (val < min_top)
                min_top = val;
                vcl = i;
                voff = idx(off);
            end
        end
    end
    
    if min_top ~= Inf
        % at least one row found
        more_rows = true;
    end
    
    %follow neighbours to get bounding box
    while more_rows && max_num > 0
    
        %there may be left neighbours that aren't as high up the page.
        [vcl, voff] = find_left(Clust, vcl, voff);
    
        l = Clust(vcl).pos(voff,1);
        t = Clust(vcl).pos(voff,2);
        r = Clust(vcl).pos(voff,3);
        b = Clust(vcl).pos(voff,4);
        hcl = vcl;
        hoff = voff;
        more_cols = true;
        more_rows = false;
        min_top = Inf;
        first_comp = Clust(hcl).comp(hoff);
        fprintf('new row found');
    
        while more_cols
            %should we expand the box dimensions?
            fprintf('.');
            if Clust(hcl).pos(hoff,2) < t
                t = Clust(hcl).pos(hoff,2);
            end
            if Clust(hcl).pos(hoff,3) > r
                r = Clust(hcl).pos(hoff,3);
            end
            if Clust(hcl).pos(hoff,4) > b
                b = Clust(hcl).pos(hoff,4);
            end
    
            %have we found a new closer next row?
            if Clust(hcl).nb(hoff,4) ~=0
                [newvcl, newvoff] = get_cl_off(Clust, Clust(hcl).nb(hoff,4));
                if Clust(newvcl).pos(newvoff,2) < min_top
                    [newfirstcl, newfirstoff] = find_left(Clust, newvcl, ...
                                                newvoff);
                    if Clust(newfirstcl).comp(newfirstoff) ~= first_comp
                        %the first component's bottom neighbour could be on the 
                        %same line, so we scan across to ensure that they have 
                        %different far right neighbours to ensure we don't 
                        %repeat the same line
                        if Clust(hcl).comp(hoff) ~= first_comp || ...
                           any(find_right(Clust, hcl, hoff) ~= ...
                               find_right(Clust, newvcl, newvoff))
                            more_rows = true;
                            min_top = Clust(newvcl).pos(newvoff,2);
                            vcl  = newvcl;
                            voff = newvoff;
                        end
                    else
                        %since on the same line, see if we should expand the
                        %right, or bottom dimensions based on those of the
                        %bottom neighbour
                        if Clust(newvcl).pos(newvoff,3) > r
                            r = Clust(newvcl).pos(newvoff,3);
                        end
                        if Clust(newvcl).pos(newvoff,4) > b
                            b = Clust(newvcl).pos(newvoff,4);
                        end
                    end
                end
            end

            %is there a lower top neighbour in the same row?
            if Clust(hcl).nb(hoff,2) ~= 0
                [newvcl, newvoff] = get_cl_off(Clust, Clust(hcl).nb(hoff,2));
                if Clust(newvcl).pos(newvoff,2) < t
                    [newfirstcl, newfirstoff] = find_left(Clust, newvcl, ...
                                                newvoff);
                    if Clust(newfirstcl).comp(newfirstoff) == first_comp
                        t = Clust(newvcl).pos(newvoff,2);
                    end
                end
            end
    
            %see if there are more items to the right in this row
            if Clust(hcl).nb(hoff,3) ~= 0
                [hcl, hoff] = get_cl_off(Clust, Clust(hcl).nb(hoff,3));
            else
                more_cols = false;
                fprintf('\n');
            end
        end
    
        %update the line boundaries
        Pos{p} = [Pos{p}; l, t, r, b];

        max_num = max_num - 1;
    end

    %draw the bounding boxes found
    if display_images || save_images
        M = label2rgb(Comps{p}, 'white', 'k');
        for i=1:size(Pos{p},1)
            l = Pos{p}(i,1); t = Pos{p}(i,2); r = Pos{p}(i,3); b = Pos{p}(i,4);
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
            fprintf('saving page %d to disk\n', p);
            imwrite(M, [img_prefix, num2str(p), '.', img_format], img_format);
        end
    end
end


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the left-most cluster on this line
function [cls, off] = find_left(Clust, cls, off)
while Clust(cls).nb(off,1) ~= 0
    [cls, off] = get_cl_off(Clust, Clust(cls).nb(off,1));
end


% get the right-most cluster on this line
function [cls, off] = find_right(Clust, cls, off)
while Clust(cls).nb(off,3) ~= 0
    [cls, off] = get_cl_off(Clust, Clust(cls).nb(off,3));
end
