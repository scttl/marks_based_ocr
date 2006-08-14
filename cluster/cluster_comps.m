function [Clust, Comps] = cluster_comps(Files)
% CLUSTER_COMPS   Cluster together connected component blobs in the files passed
%
%   [Clust, Comps] = CLUSTER_COMPS(FILES)
%   FILES should either be a string, or a cell array of strings listing the
%   path and name of the file(s) to be processed.  We assume that each file
%   has a binary (black and white) intensity.
%
%   Clust is a struct meant to hold information on all clusters found while
%   processing FILES.  It contains the following fields:
%     num   - a scalar denoting the number of clusters in the struct.
%     num_comps - a vector giving the number of components belonging to each
%                 cluster.
%     comps - a cell array denoting for each cluster, which components belong
%             to it (listed as a vector of component id's).  Each component id 
%             is specified by its row number in the fields of the struct.
%     avg - a cell array, where each entry is a numeric matrix that gives the 
%           average pixel intensity of each pixel based on the components that 
%           belong to that cluster
%     norm_sq - a vector specifying the L2 norm squared of each avg (used to 
%               speed up distance calculations)
%     refined - a boolean vector used to determine whether a cluster must be
%               further processed
%     changed - a boolean vector used to determine whether a cluster has
%               changed (gained or lost components) during a refinement pass.
%     offset - this holds the most frequently occuring relative position of the
%              bottom edge of the component from the baseline.  Positive values
%              refer to the number of pixels below the baseline, the bottom of
%              the image lies, and negative refers to pixels above the baseline.
%              Note that this is not calculated in this function (thus 
%              found_offsets is false -- see below).
%     found_offsets - this is a boolean that will be set to true if offset 
%                     values for each cluster have been calculated and false
%                     otherwise.
%     num_trans - a scalar that determines the number of left-right component
%                 transitions found.  Note that this is not calculated here
%     bigram - a num x num matrix that denotes (smoothed) transition 
%              probabilities.  Again, this isn't calculated here.
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
%     offset - this vector will hold the relative position of the bottom edge
%              of the component as the number of pixels below (if positive) or
%              above (if negative) from the baseline of the line to which it 
%              was found.  Note that this is not calculated in this function.
%     found_offsets - this is a boolean that will be set to true if offset 
%                     values for each cluster have been calculated and false
%                     otherwise.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: cluster_comps.m,v 1.5 2006-08-14 01:40:53 scottl Exp $
%
% REVISION HISTORY
% $Log: cluster_comps.m,v $
% Revision 1.5  2006-08-14 01:40:53  scottl
% added new Clust.changed field, finished implementing merges.
%
% Revision 1.4  2006/07/05 01:14:30  scottl
% rewritten based on new Cluster and Component structures.
%
% Revision 1.3  2006/06/19 20:56:05  scottl
% bugfix: better handling of neighbours of nested components
%
% Revision 1.2  2006/06/12 20:56:01  scottl
% changed comps to a cell array (indexed by page) to allow different sized pages
%
% Revision 1.1  2006/06/03 20:55:47  scottl
% Initial check-in.
%
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%the number of directions to consider when looking for conn. components (4 or
%8) are the only valid entries
num_dirs = 8;
%num_dirs = 4;

%these parameters control the distance metric used during grouping.  Valid
%choices are: 'euc', 'conv', or 'hausdorff' for Euclidian, Convolutional
%Euclidian, or Hausdorff (with underlying Euclidian) distance measurements
%used.
match_metric='hausdorff';
%match_metric='euc';
scale_metric='euc';
split_metric='euc';

%the following parameters control when things are grouped together.
straight_match_thresh = .009;
match_thresh = 1.3;  %hausdorff euclidian distance thresh
%match_thresh = .009;
scale_thresh = .009;
split_thresh = .013;  %good euc thresh
max_splits = 2;
merge_thresh = 3;  %look for compoonents no more than 3 pixels apart
merge_pct = .85;  %ensure at least 85% of the elements in each cluster match
merge_min_comps = 3; %ensure that the cluster has at least 3 elements 

%the minimum size (in number of pixels) of an element for it to be considered
%for reference: dot's on i's are about 3x3 on NIPS papers
min_elem_width = 4;
min_elem_height = 4;

%the maximum size (in number of pixels) of an element for it to be considered
%for reference: typical text characters in NIPS papers are about 10-30x10-30
max_elem_width = 100;
max_elem_height = 100;

%should we remove certain regions if an associated jtag file exits?
crop_regions = true;
%give a cell array list of the regions to be removed.  Possible choices are:
%section_heading, subsection_heading, footer, references, bullet_item, 
%table_label, header, authour_list, code_block, main_title, figure_label, 
%figure, image, text, equation, footnote, figure_caption, decoration, abstract,
%table, graph, eq_number, editor_list, table_caption,  pg_number
rem_region_list = {'figure', 'figure_label', 'image', 'decoration', 'table', ...
    'code_block', 'table_label', 'graph', 'equation', 'eq_number'};

%should we display the page image before and after processing?  This uses
%memory and takes longer to run
display_page = false;

%RGB colour values used to display updates if set to true above
bg_col = reshape([0,0,0],1,1,3);  % this is black
out_col = reshape([255,0,0],1,1,3);  % this is red
bg_val = 0;  %used in Comps
fg_val = 1;  %used in cluster averages initially

%should we save connected component images?  This uses memory and takes longer
%to run
save_cc_page = false;
%prefix of filename to use when saving connected component images
img_prefix = 'results/nips5_comps';
%file format to use when saving the image
img_format = 'png';


% CODE START %
%%%%%%%%%%%%%%
tic;
%the number of pages of elements to cluster (we assume 1 page per file passed)
num_pgs = size(Files,1);

%initialize the Comps and Clust structs
Clust.num = uint16(0);
Clust.num_comps = [];
Clust.comps = {};
Clust.avg = {};
Clust.norm_sq = [];
Clust.refined = logical([]);
Clust.changed = logical([]);
Clust.offset = int16([]);
Clust.found_offsets = false;
Clust.num_trans = 0;
Clust.bigram = [];

Comps.max_comp = 0;
if num_pgs == 1
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
Comps.offset = int16([]);
Comps.found_offsets = false;


%this is used to determine when to use straight-matching
pass = 1;

%iterate through each page adding items to our cluster list
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
        %remove the regions listed based on the associated jtag file
        %(we assume the jtag file is in the same directory as the image file,
        % and has the same name, but with the extension .jtag instead)
        dot_pos = strfind(Comps.files{pp}, '.');
        if isempty(dot_pos)
            warning('file has no extension, unable to parse regions');
        else
            jtag_file = [Comps.files{pp}(1:dot_pos(end)), 'jtag'];
            pos = jtag_region_finder(jtag_file, rem_region_list);
            for i=1:size(pos,1)
                Lbl_img(pos(i,2):pos(i,4), pos(i,1):pos(i,3)) = bg_val;
            end
        end
        fprintf('%.2fs: done cropping unwanted regions from this page\n', ...
                toc);
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
    Comps.pos = [Comps.pos; zeros(num_new_comps,4)];

    %first we need to fill in the L T R B pixel boundary values for each 
    %component
    [R, C] = find(Lbl_img ~= 0);
    for ii=1:length(R)
        comp = Lbl_img(R(ii), C(ii));
        loc = comp + prev_max_comp;
        %is this the first time seeing this component?
        if Comps.pos(loc,1) == 0
            %set all 4 positions for this component to this location
            Comps.pos(loc,1) = C(ii); Comps.pos(loc,2) = R(ii); 
            Comps.pos(loc,3) = C(ii); Comps.pos(loc,4) = R(ii);
        else
            %update the R position and see if we should update T or B
            Comps.pos(loc,3) = C(ii);
            if R(ii) < Comps.pos(loc,2)
                Comps.pos(loc,2) = R(ii);
            end
            if R(ii) > Comps.pos(loc,4)
                Comps.pos(loc,4) = R(ii);
            end
        end
    end
    fprintf('%.2fs: Position of each component found\n', toc);


    %prune any components that don't meet the size requirements
    start_comp = prev_max_comp + 1;
    comp_idcs = start_comp:Comps.max_comp;
    V_diff = Comps.pos(comp_idcs,4) - Comps.pos(comp_idcs,2);
    H_diff = Comps.pos(comp_idcs,3) - Comps.pos(comp_idcs,1);
    reject_list = prev_max_comp + find(V_diff < min_elem_height | ...
                   H_diff < min_elem_width | V_diff > max_elem_height | ...
                   H_diff > max_elem_width);
    for ii=1:length(reject_list)
        %also remove them from the component image
        pos = Comps.pos(reject_list(ii),:);
        region = Lbl_img(pos(2):pos(4), pos(1):pos(3));
        idx = find(region == (reject_list(ii) - prev_max_comp));
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
        idx = find(region == (keep_list(ii - prev_max_comp) - prev_max_comp));
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
    fprintf('%.2fs: Neighbours of each component found\n', toc);


    %initialize the remaining Comps fields
    Comps.pg(comp_idcs,1) = pp;
    clust_idcs = Clust.num+1:Clust.num+num_new_comps;
    Comps.clust(comp_idcs,1) = clust_idcs;
    Comps.offset(comp_idcs,1) = 0;
    fprintf('%.2fs: Remaining component fields initialized\n', toc);


    %add each component as a new cluster
    Clust.num = Clust.num + num_new_comps;
    Clust.num_comps(clust_idcs,1) = 1;
    Clust.comps(clust_idcs,1) = num2cell(comp_idcs');
    Clust.avg(clust_idcs,1) = cell(num_new_comps,1);
    Clust.norm_sq(clust_idcs,1) = 0;
    Clust.refined(clust_idcs,1) = false;
    Clust.changed(clust_idcs,1) = false;
    Clust.offset(clust_idcs,1) = 0;

    %update the avg and norm_sq fields
    for ii=comp_idcs
        l = Comps.pos(ii,1); t = Comps.pos(ii,2); r = Comps.pos(ii,3); ...
        b = Comps.pos(ii,4);
        cc = Comps.clust(ii);
        Clust.avg{cc} = zeros(b-t+1, r-l+1);
        on_idx = find(Lbl_img(t:b,l:r) == ii);
        %the norm_squared value is the number of 'on' pixels since each gets
        %an initial intensity of 1.
        Clust.norm_sq(cc) = length(on_idx);
        Clust.avg{cc}(on_idx) = fg_val;

        if display_page || save_cc_page
            Rgb_img(t,l:r,:) = repmat(out_col,1,r-l+1);
            Rgb_img(t:b,l,:) = repmat(out_col,b-t+1,1);
            Rgb_img(t:b,r,:) = repmat(out_col,b-t+1,1);
            Rgb_img(b,l:r,:) = repmat(out_col,1,r-l+1);
        end
    end
    if display_page
        imshow(Rgb_img);
        fprintf('Ready to start matching.  Press a key to begin\n');
        pause;
    end
    if save_cc_page
        fprintf('%.2fs: Writing components of page %d to disk\n', toc, p);
        imwrite(Rgb_img, [img_prefix, num2str(p), '.', img_format], img_format);
    end
    fprintf('%.2fs: New components added to individual new clusters\n', toc);

    %perform refinements until the number of clusters remains constant
    while true
        num_clusts = Clust.num;

        %First look for straight matches among the new clusters.  Initially this
        %should drastically reduce the number of clusters.  Just use simple
        %euclidian distance (to keep run-time low).
        if pp ~= num_pgs || pass == 1
            fprintf('\n\n%.2fs: Starting straight-match refine pass\n', toc);
            [Clust, Comps] = match_refine(Clust, Comps, 'euc', ...
                             straight_match_thresh);
            if pp == num_pgs
                pass = 2;
                Clust.refined(:) = false;
                Clust.changed(:) = false;
            end
        else
            %since all pages have been processed with straight matches
            %we now further refine using splits, merges, and matches using
            %the appropriate distance metrics
            if pass == 2
                %first time through we want to attempt to match over all
                %clusters
                Clust.refined(:) = false;
            end
            fprintf('\n\n%.2fs: Starting %s match refine pass\n', toc, ...
                    match_metric);
            [Clust, Comps] = match_refine(Clust, Comps, match_metric, ...
                             match_thresh);

            if pass == 2
                %first time through we want to attempt to merge over all
                %clusters
                Clust.refined(:) = false;
            end
            fprintf('\n\n%.2fs: Starting merge refine pass\n', toc);
            [Clust, Comps] = merge_refine(Clust, Comps,  merge_thresh, ...
                             merge_pct, merge_min_comps);

            if pass == 2
                %first time through we want to attempt to split over all
                %clusters
                Clust.refined(:) = false;
            end
            fprintf('\n%.2fs: Starting horizontal split refine pass\n', toc);
            [Clust, Comps] = split_refine(Clust, Comps, split_metric, ...
                                       max_splits, split_thresh);

            %to speed up subsequent refinements we only need to further process
            %those Clusters that have changed during the last round of
            %refinements
            Clust.refined(Clust.changed == true) = false;
            Clust.changed(:) = false;
            pass = pass + 1;
        end

        if num_clusts == Clust.num;
            fprintf('\n%.2fs: no further reduction in cluster size: %d\n', ...
                    toc, Clust.num);
            break;
        end
    end
end

%perform a final scale refine over the remaining elements
%fprintf('\n%.2fs: Starting scale refine pass\n', toc);
%[Clust, new_list] = scale_refine(Clust, Comps, scale_metric, scale_thresh);

%finally, sort the clusters
[Clust, Comps] = sort_clusters(Clust, Comps);

fprintf('TOTAL NUMBER OF COMPONENTS FOUND: %d\n', Comps.max_comp);
fprintf('TOTAL NUMBER OF CLUSTERS FOUND: %d\n', Clust.num);



% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
