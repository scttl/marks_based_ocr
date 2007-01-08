function Comps = get_comps(Files, varargin)
% GET_COMPS   Determine connected component blobs in the files passed
%
%   Comps = GET_COMPS(FILES, ['var1', new_val1]...)
%   FILES should either be a string, or a cell array of strings listing the
%   path and name of the file(s) to be processed.  We assume that each file
%   has a binary (black and white) intensity.  One can also override the
%   default value of any of the local variables defined below by passing in
%   pairs of optional arguments.  The first should be a string specifying the
%   name of the variable to override (typed exactly as seen below).  The second
%   should be the value to replace the original with.
%
%   Comps is a struct meant to hold information on all the connected components
%   found while processing FILES.  It contains the following fields:
%     max_comp - integer specifying the largest component number index.
%     files - a cell array of strings listing the path and name of the file(s)
%             processed.
%     pg_size - a 2 column array giving the number of rows and columns of each
%               page processed.
%     clust - a vector listing the index into the Clust array to which 
%             each component belongs.
%     pos - an nx4 numeric matrix listing the pixel co-ordinates of the left,
%           top, right, and bottom edges of the tight bounding box around the
%           component on the page it was found.  There is one such row for each
%           connected component.
%     pg - a vector of length n giving the page number upon which the component
%          was found.  This number can be used to index Comps.files and 
%          determine the filename corresponding to this page.
%     nb - an nx4 numeric matrix listing its left, top, right, and bottom 
%          nearest neighbours by component id.
%     nb_dist - an nx4 numeric matrix listing its left, top, right, and bottom 
%               nearest neighbours distances (in pixels)
%     regions - this mx6 matrix lists the page number, region number, and 4 
%               positional co-ordinates (left,top,right, and bottom) of each 
%               region found if an associated jtag file is being used to crop 
%               regions.  If this isn't the case, this is left empty.
%     found_lines - this is a boolean that will be set to true if the line to
%                   which each component belongs is found, and the parameters
%                   below have been set.  It is initially false, and the
%                   parameters below are initially emtpy.  They can be
%                   calculated with a call to get_lines() though.
%     line - a vector listing the index into the Line struct to which this
%            component belongs
%     modal_height - a scalar value listing the most commonly found height (as
%                    measured from the baseline to the x-height) of the
%                    components
%     scale_factor - a vector listing how much the component must be scaled up
%                    or down so that its height matches the modal height.  A
%                    value > 1 means the component must be increased in size to
%                    match.  A value = 1 means no change, and a value < 1 
%                    requires a decrease in size.
%     descender_off - this vector will hold the relative position of the 
%                     bottom edge of the component as the number of pixels 
%                     below (if positive) or above (if negative) from the 
%                     baseline of the line to which it was found.
%     ascender_off - analogous to the descender_offset, this vector holds
%                    the relative position of the top edge of the component
%                    as the number of pixels above (if positive) or below
%                    (if negative) from the x-height of the line to which it
%                    was found.
%     found_true_labels - this is a boolean that will be set to true if the
%                         ground truth labels for each component have been
%                         caluculated.  This can be done with a call to
%                         ground_truth_label()
%     truth_label - this vector will hold decimal ASCII character values for
%                   each character (or partial character) that has been asigned.
%     model_spaces - this boolean will be set to true if spaces have been
%                    modelled and counts taken.  See add_space_model()


% CVS INFO %
%%%%%%%%%%%%
% $Id: get_comps.m,v 1.9 2007-01-08 22:03:39 scottl Exp $
%
% REVISION HISTORY
% $Log: get_comps.m,v $
% Revision 1.9  2007-01-08 22:03:39  scottl
% added region number to the region information (so they can be associated
% with the corect ground truth file)
%
% Revision 1.8  2007-01-05 17:10:21  scottl
% added region field to Comps structure, useful for selecting and ordering
% text from individual regions (for ground truth comparison etc.)
%
% Revision 1.7  2007-01-02 19:23:12  scottl
% changed cropping so that regions are kept instead of removed.
%
% Revision 1.6  2006-12-17 19:52:26  scottl
% update method of calculating neighbour distance based on boudning box
% distances (instead of ink to ink distances).  Also fix a couple of
% mis-references to the page variable.
%
% Revision 1.5  2006-11-25 20:06:54  scottl
% small addition of model_spaces field.
%
% Revision 1.4  2006-11-07 02:52:06  scottl
% no change.
%
% Revision 1.3  2006/10/18 15:52:10  scottl
% small fix for optional argument passing.
%
% Revision 1.2  2006-10-09 16:30:39  scottl
% added parameters for modal height, scale factor, and groundtruth label.
%
% Revision 1.1  2006-09-22 18:04:51  scottl
% initial check-in.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%the number of directions to consider when looking for conn. components (4 or
%8) are the only valid entries
num_dirs = 8;

%the maximum allowable aspect ratio (difference between width and height).
max_ar_diff = 10.0;

%the minimum size (in number of pixels) of an element for it to be considered
%for reference: dot's on i's are about 3x3 on NIPS papers
min_elem_width = 2;
min_elem_height = 2;

%the maximum size (in number of pixels) of an element for it to be considered
%for reference: typical text characters in NIPS papers are about 10-30x10-30
max_elem_width = 100;
max_elem_height = 100;

%should we remove certain regions if an associated jtag file exits?
crop_regions = true;
%give a cell array list of the regions to be kept.  Possible choices are:
%section_heading, subsection_heading, footer, references, bullet_item, 
%table_label, header, authour_list, code_block, main_title, figure_label, 
%figure, image, text, equation, footnote, figure_caption, decoration, abstract,
%table, graph, eq_number, editor_list, table_caption,  pg_number
keep_region_list = {'section_heading', 'subsection_heading', 'footer', ...
     'references', 'bullet_item', 'table_label', 'header', 'authour_list', ...
     'code_block', 'main_title', 'figure_label', 'text', 'footnote', ...
     'figure_caption', 'abstract', 'editor_list', 'table_caption', 'pg_number'};
jtag_extn = 'jtag';  %file extension for jtag files

%should we display the page image before and after processing?  This uses
%memory and takes longer to run
display_page = false;

%RGB colour values used to display updates if set to true above
bg_col = reshape([0,0,0],1,1,3);  % this is black
out_col = reshape([255,0,0],1,1,3);  % this is red

%conventional pixel value to use for background intensities
bg_val = 0; 

%should we save connected component images?  This uses memory and takes longer
%to run
save_cc_page = false;
global MOCR_PATH;
%prefix of filename to use when saving connected component images.  If multiple
%files are passed, the page numbers are also added to the filename
img_prefix = [MOCR_PATH, '/results/conn_comps'];
%file format to use when saving the image
img_format = 'png';


% CODE START %
%%%%%%%%%%%%%%
tic;

%process arguments
if nargin < 1
    error('must specify at least one file to process!');
elseif nargin > 1
    process_optional_args(varargin{:});
end


%the number of pages of elements to cluster (we assume 1 page per file passed)
num_pgs = size(Files,1);

%initialize the Comps struct
Comps = init_comps(Files);

%iterate through each page adding items to our list of components
for pp=1:num_pgs

    fprintf('Processing page %d\n', pp);

    %begin by loading the page contents into memory.  Since Matlab convention
    %is 1 for background intensities, flip this.
    Lbl_img = imread(Comps.files{pp});
    if bg_val ~= 1
        Lbl_img = ~Lbl_img;
    end
    Comps.pg_size(pp,:) = size(Lbl_img);

    if crop_regions
        %keep only the regions listed based on the associated jtag file (if one
        %exists)
        dot_pos = strfind(Comps.files{pp}, '.');
        if isempty(dot_pos)
            warning('MBOCR:noExtn', ...
                    'file has no extension. Unable to parse regions');
        else
            jtag_file = [Comps.files{pp}(1:dot_pos(end)), jtag_extn];
            [Lbl_img, r_num, r] = crop_jtag_regions(jtag_file,Lbl_img, ...
                                  keep_region_list);
            Comps.regions = [Comps.regions; [repmat(pp,size(r,1),1),r_num,r]];
            fprintf('%.2fs: done cropping unwanted regions from this page\n',...
                    toc);
        end
    end

    %determine the connected components on this page
    [Lbl_img, num_new_comps] = bwlabel(Lbl_img, num_dirs);
    fprintf('%.2fs: %d connected components found on this page\n', toc, ...
            num_new_comps);

    if display_page || save_cc_page
        if num_new_comps > 0
            %this colormap uses the black background with white letter font
            Rgb_img = label2rgb(Lbl_img, 'white', 'k');
        else
            Rgb_img =  repmat(bg_col, size(Lbl_img));
        end
        if display_page
            imshow(Rgb_img);
            drawnow;
            hold on;
            fprintf('%.2fs: Image of connected components rendered\n', toc);
        end
    end

    %determine the positions of the components on this page
    prev_max_comp = Comps.max_comp;
    Comps.max_comp = Comps.max_comp + num_new_comps;
    Comps.pos = [Comps.pos; get_comp_bb(Lbl_img, num_new_comps)];
    fprintf('%.2fs: Position of each component found\n', toc);


    %prune any components that don't meet the size requirements
    start_comp = prev_max_comp + 1;
    comp_idcs = start_comp:Comps.max_comp;
    V_diff = Comps.pos(comp_idcs,4) - Comps.pos(comp_idcs,2) + 1;
    H_diff = Comps.pos(comp_idcs,3) - Comps.pos(comp_idcs,1) + 1;
    VH_ratio = double(H_diff) ./ double(V_diff);
    reject_list = prev_max_comp + find(V_diff < min_elem_height | ...
                   H_diff < min_elem_width | V_diff > max_elem_height | ...
                   H_diff > max_elem_width | VH_ratio < 1/max_ar_diff | ...
                   VH_ratio > max_ar_diff);
    for ii=1:length(reject_list)
        %also remove them from the component image
        pos = Comps.pos(reject_list(ii),:);
        region = Lbl_img(pos(2):pos(4), pos(1):pos(3));
        idx = region == (reject_list(ii) - prev_max_comp);
        region(idx) = bg_val;
        Lbl_img(pos(2):pos(4), pos(1):pos(3)) = region;
    end
    %squeeze out the rejected components, and update the numbers in the label
    %image
    keep_list = setdiff(comp_idcs, reject_list);
    Comps.pos = Comps.pos([1:prev_max_comp, keep_list],:);
    reg_len = length(reject_list);
    num_new_comps = num_new_comps - reg_len;
    Comps.max_comp = Comps.max_comp - reg_len;
    comp_idcs = start_comp:Comps.max_comp;
    for ii=comp_idcs
        %update the label image
        pos = Comps.pos(ii,:);
        region = Lbl_img(pos(2):pos(4), pos(1):pos(3));
        idx = region == (keep_list(ii - prev_max_comp) - prev_max_comp);
        region(idx) = ii;
        Lbl_img(pos(2):pos(4), pos(1):pos(3)) = region;
    end
    fprintf('%.2fs: Rejected Components pruned\n', toc);

    %determine the nearest neighbour for each component in each direction
    Comps.nb = [Comps.nb; zeros(num_new_comps,4)];
    Comps.nb_dist = [Comps.nb_dist; Inf(num_new_comps, 4)];
    %first get the top and bottom distances
    prev_loc  = 0;
    prev_row  = Inf;
    last_col  = NaN;
    [R, C] = find(Lbl_img ~= 0);
    for ii=1:length(R)
        loc = Lbl_img(R(ii), C(ii));
        if last_col == C(ii) && prev_loc ~= loc
            %not first valid component encountered in this column 
            val = R(ii) - prev_row;
            if val < Comps.nb_dist(loc,2) && ...
               Comps.pos(prev_loc,2) < Comps.pos(loc,2)
                %new closest top neighbour for comp.  The second test ensures
                %we cant end up with self-pointing neighbours due to nested 
                %components (ex. a C with a dot in its center), which would have
                %C point to the dot as its closest top neighbour which 
                %points to C as its top neighbour etc.
                Comps.nb_dist(loc,2) = val;
                Comps.nb(loc,2) = prev_loc;
            elseif val==Comps.nb_dist(loc,2) && Comps.nb(loc,2) ~= prev_loc ...
                   && Comps.pos(prev_loc,2) < Comps.pos(loc,2)
                %2nd component the same distance to comp.  Choose the
                %component with larger L-R overlap with comp.
                pn_loc = Comps.nb(loc,2);
                pn_ovr = min(Comps.pos(loc,3), Comps.pos(pn_loc,3)) - ...
                         max(Comps.pos(loc,1), Comps.pos(pn_loc,1));
                pl_ovr = min(Comps.pos(loc,3), Comps.pos(prev_loc,3)) - ...
                         max(Comps.pos(loc,1), Comps.pos(prev_loc,1));
                if pl_ovr > pn_ovr
                    Comps.nb(loc,2) = prev_loc;
                end
            end
            if val < Comps.nb_dist(prev_loc,4) && ...
               Comps.pos(loc,4) > Comps.pos(prev_loc,4)
                %loc is new closest bottom neighbour for prev_loc
                Comps.nb_dist(prev_loc,4) = val;
                Comps.nb(prev_loc,4) = loc;
            elseif val==Comps.nb_dist(prev_loc,4) && ...
                   Comps.nb(prev_loc,4) ~= loc && ...
                   Comps.pos(loc,4) > Comps.pos(prev_loc,4)
                %2nd component the same distance to prev_loc.  Choose the
                %component with larger L-R overlap with prev_loc.
                pn_loc = Comps.nb(prev_loc,4);
                pn_ovr = min(Comps.pos(prev_loc,3), Comps.pos(pn_loc,3)) - ...
                         max(Comps.pos(prev_loc,1), Comps.pos(pn_loc,1));
                l_ovr = min(Comps.pos(prev_loc,3), Comps.pos(loc,3)) - ...
                        max(Comps.pos(prev_loc,1), Comps.pos(loc,1));
                if l_ovr > pn_ovr
                    Comps.nb(prev_loc,4) = loc;
                end
            end
        end
        %update component, row, and column accordingly
        prev_loc  = loc;
        prev_row  = R(ii);
        last_col  = C(ii);
    end

    %now get the left and right distances.  This requires working with the
    %transpose to ensure we get items in L-R order
    [C, R] = find(Lbl_img' ~= 0);
    prev_loc  = 0;
    prev_col  = Inf;
    last_row  = NaN;
    for ii=1:length(C)
        loc = Lbl_img(R(ii), C(ii));
        if last_row == R(ii) && prev_loc ~= loc
            %not first valid component encountered in this row 
            val = C(ii) - prev_col;
            if val < Comps.nb_dist(loc,1) && ...
               Comps.pos(prev_loc,1) < Comps.pos(loc,1)
                %new closest left neighbour for loc
                Comps.nb_dist(loc,1) = val;
                Comps.nb(loc,1) = prev_loc;
            elseif val==Comps.nb_dist(loc,1) && Comps.nb(loc,1) ~= prev_loc ...
                   &&  Comps.pos(prev_loc,1) < Comps.pos(loc,1)
                %2nd component the same distance to loc.  Choose the
                %component with larger T-B overlap with loc.
                pn_loc = Comps.nb(loc,1);
                pn_ovr = min(Comps.pos(loc,4), Comps.pos(pn_loc,4)) - ...
                         max(Comps.pos(loc,2), Comps.pos(pn_loc,2));
                pl_ovr = min(Comps.pos(loc,4), Comps.pos(prev_loc,4)) - ...
                         max(Comps.pos(loc,2), Comps.pos(prev_loc,2));
                if pl_ovr > pn_ovr
                    Comps.nb(loc,1) = prev_loc;
                end
            end
            if val < Comps.nb_dist(prev_loc,3) && ...
               Comps.pos(loc,3) > Comps.pos(prev_loc,3)
                %loc is new closest right neighbour for prev_loc
                Comps.nb_dist(prev_loc,3) = val;
                Comps.nb(prev_loc,3) = loc;
            elseif val==Comps.nb_dist(prev_loc,3) && ...
                   Comps.nb(prev_loc,3) ~= loc && ...
                   Comps.pos(loc,3) > Comps.pos(prev_loc,3)
                %2nd component the same distance to prev_loc.  Choose the
                %component with larger T-B overlap with prev_loc.
                pn_loc = Comps.nb(prev_loc,3);
                pn_ovr = min(Comps.pos(prev_loc,4), Comps.pos(pn_loc,4)) - ...
                         max(Comps.pos(prev_loc,2), Comps.pos(pn_loc,2));
                l_ovr = min(Comps.pos(prev_loc,4), Comps.pos(loc,4)) - ...
                        max(Comps.pos(prev_loc,2), Comps.pos(loc,2));
                if l_ovr > pn_ovr
                    Comps.nb(prev_loc,3) = loc;
                end
            end
        end
        %update loc, column and row accordingly
        prev_loc  = loc;
        prev_col  = C(ii);
        last_row  = R(ii);
    end
    clear C R;
    %update the neighbour distances based on bounding boxes
    for ii = 1:4
        idx = find(Comps.nb(:,ii) ~= 0);
        if ii<=2
            alt = ii+2;
            nbs = Comps.nb(idx,ii);
            Comps.nb_dist(idx,ii) = Comps.pos(idx,ii) - Comps.pos(nbs,alt);
        else
            alt = ii-2;
            nbs = Comps.nb(idx,ii);
            Comps.nb_dist(idx,ii) = Comps.pos(nbs,alt) - Comps.pos(idx,ii);
        end
    end
    fprintf('%.2fs: Neighbours of each component found\n', toc);

    %initialize the remaining Comps fields
    Comps.pg(comp_idcs,1) = pp;

    %show and/or save the appropriate RGB image if required.
    if display_page || save_cc_page
        for ii=comp_idcs
            l = Comps.pos(ii,1); t = Comps.pos(ii,2); r = Comps.pos(ii,3); ...
            b = Comps.pos(ii,4);
    
            Rgb_img(t,l:r,:) = repmat(out_col,1,r-l+1);
            Rgb_img(t:b,l,:) = repmat(out_col,b-t+1,1);
            Rgb_img(t:b,r,:) = repmat(out_col,b-t+1,1);
            Rgb_img(b,l:r,:) = repmat(out_col,1,r-l+1);
        end
        if display_page
            imshow(Rgb_img);
            fprintf('%.2fs: Displaying image. Press a key to continue\n', toc);
            pause;
        end
        if save_cc_page
            fprintf('%.2fs: Writing components of page %d to disk\n', toc, pp);
            imwrite(Rgb_img, [img_prefix, num2str(pp), '.', img_format], ...
                    img_format);
        end
    end
end

fprintf('%.2fs: TOTAL NUMBER OF COMPONENTS FOUND: %d\n', toc, Comps.max_comp);



% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this function creates an empty component stucture
function Comps = init_comps(Files)

Comps.max_comp = 0;
if ~iscell(Files)
    Comps.files = {Files};
else
    Comps.files = Files;
end
Comps.pg_size = uint16([]);
Comps.clust = uint16([]);
Comps.pos = uint16([]);
Comps.pg = uint32([]);
Comps.nb = [];
Comps.nb_dist = uint16([]);
Comps.regions = [];
Comps.found_lines = false;
Comps.line = uint64([]);   %uint32([]);
Comps.modal_height = NaN;
Comps.scale_factor = [];
Comps.descender_off = int16([]);
Comps.ascender_off = int16([]);
Comps.found_true_labels = false;
Comps.truth_label = cell(0);
Comps.model_spaces = false;
