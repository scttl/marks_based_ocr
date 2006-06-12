function [Clust, Comps] = cluster_comps(Files)
% CLUSTER_COMPS   Cluster together connected component blobs in the files passed
%
%   [Clust, Comps] = CLUSTER_COMPS(FILES)
%   FILES should either be a string, or a cell array of strings listing the
%   path and name of the file(s) to be processed.  We assume that each file
%   has a binary (black and white) intensity, and that each file is the
%   same size (number of pixels).
%
%   Clust is a cell array, with each entry containing a vector of all the
%   connected components found to be equivalent after processing
%
%   Comps is a cell array of connected component labels, one page per entry.

% CVS INFO %
%%%%%%%%%%%%
% $Id: cluster_comps.m,v 1.2 2006-06-12 20:56:01 scottl Exp $
%
% REVISION HISTORY
% $Log: cluster_comps.m,v $
% Revision 1.2  2006-06-12 20:56:01  scottl
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

%this parameter controls the distance metric used during grouping.  Valid
%choices are: 'euc', 'conv', or 'hausdorff' for Euclidian, Convolutional
%Euclidian, or Hausdorff (with underlying Euclidian) distance measurements
%used.
metric='euc';

%the following parameters control when things are grouped together.
match_thresh = .009;
scale_thresh = .009;
split_thresh = .013;
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
save_cc_page = true;
%prefix of filename to use when saving connected component images
img_prefix = 'results/nips5_comps';
%file format to use when saving the image
img_format = 'png';

%the number of pages of elements to cluster (we assume 1 page per file passed)
num_pgs = size(Files,1);

%the total number of components found over all pages examined
tot_comps = 0;

%this is used to determine when to use straight-matching
pass = 1;


% CODE START %
%%%%%%%%%%%%%%

%begin by loading the file contents into memory.  Each file is a 2 dimensional 
%matrix of binary pixel intensities.  We create a 3d array Comps with one 
%such page matrix per entry.
tic;
if num_pgs > 1
    for p=1:num_pgs
        %since Matlab convention is 1 for background intensities, flip this and
        %convert to a single precision matrix since we will renumber the 
        %foreground pixels with component numbers
        Comps{p} = single(~ imread(Files{p}));
    end
else
    Comps{1} = single(~ imread(Files));
end

%initialize our cluster groupings
Clust = [];
fprintf('%.2fs: Initialization complete\n', toc);

%iterate through each page adding items to our cluster list
for p=1:num_pgs

    fprintf('Processing page %d\n', p);

    if crop_regions
        %remove the regions listed based on the associated jtag file
        %(we assume the jtag file is in the same directory as the image file,
        % and has the same name, but with the extension .jtag instead)
        if num_pgs > 1
            fl_name = Files{p};
        else
            fl_name = Files;
        end
        dot_pos = strfind(fl_name, '.');
        if length(dot_pos) == 0
            warning('file has no extension, unable to parse regions');
        else
            jtag_file = [fl_name(1:dot_pos(end)), 'jtag'];
            pos = jtag_region_finder(jtag_file, rem_region_list);
            for i=1:size(pos,1)
                Comps{p}(pos(i,2):pos(i,4), pos(i,1):pos(i,3)) = bg_val;
            end
        end
    end

    %to ensure unique component numbers, update the new components by the 
    %previous largest element
    if p ~= 1
        tot_comps = max(max(Comps{p-1}));
    end
    %start by determining the connected components on this page
    [Comps{p}, num_comps] = bwlabel(Comps{p}, num_dirs);
    fprintf('%.2fs: %d connected components found on this page\n', toc, ...
            num_comps);
    
    Comps{p} = (tot_comps .* (Comps{p} ~= 0)) + Comps{p};
    fprintf('%.2fs: Connected component numbers updated\n', toc);
    start_comp = tot_comps + 1;
    tot_comps = tot_comps + num_comps;

    if display_page || save_cc_page
        if num_comps > 0
            %this colormap uses the black background with white letter font
            M = label2rgb(Comps{p}, 'white', 'k');
        else
            M =  repmat(bg_col, size(Comps{p}));
        end
        if display_page
            imshow(M);
            drawnow;
            hold on;
            fprintf('%.2fs: Image of connected components rendered\n', toc);
        end
    end

    %first we need to fill in the L T R B pixel boundary values for each 
    %component
    Pos = zeros(num_comps, 4);
    [RR, CC] = find(Comps{p} ~= 0);
    for i=1:length(RR)
        %is this the first time seeing this component?
        comp = Comps{p}(RR(i), CC(i));
        loc = comp - (start_comp - 1);
        if Pos(loc,1) == 0
            %set all 4 positions for this component to this location
            Pos(loc,1) = CC(i); Pos(loc,2) = RR(i); Pos(loc,3) = CC(i);
            Pos(loc,4) = RR(i);
        else
            %update the R position and see if we should update T or B
            Pos(loc,3) = CC(i);
            if RR(i) < Pos(loc,2)
                Pos(loc,2) = RR(i);
            end
            if RR(i) > Pos(loc,4)
                Pos(loc,4) = RR(i);
            end
        end
    end
    %prune any components that don't meet the size requirements
    reject_list = [];
    for i=1:num_comps
        l = Pos(i,1); t = Pos(i,2); r = Pos(i,3); b = Pos(i,4);
        v_diff = (b - t);
        h_diff = (r - l);
        if (v_diff < min_elem_height) || (h_diff < min_elem_width) || ...
           (v_diff > max_elem_height) || (h_diff > max_elem_width)
            reject_list = [reject_list, i];
            %also remove them from Comps
            region = Comps{p}(t:b, l:r);
            idx = find(region == start_comp - 1 + i);
            region(idx) = bg_val;
            Comps{p}(t:b, l:r) = region;
        end
    end
    keep_list = setdiff(1:num_comps, reject_list);
    fprintf('%.2fs: Position of each component found\n', toc);

    %determine the nearest neighbour for each component in each direction
    Nb  = zeros(num_comps, 4);
    N_dist = Inf + zeros(num_comps, 4);
    %first get the top and bottom distances
    prev_comp = 0;
    prev_loc  = 0;
    prev_row  = Inf;
    last_col  = NaN;
    [RR, CC] = find(Comps{p} ~= 0);
    for i=1:length(RR)
        comp = Comps{p}(RR(i), CC(i));
        loc = comp - (start_comp - 1);
        if last_col == CC(i) && prev_comp ~= comp
            %not first valid component encountered in this column 
            val = RR(i) - prev_row;
            if val < N_dist(loc,2)
                %new closest top neighbour for comp
                N_dist(loc,2) = val;
                Nb(loc,2) = prev_comp;
            elseif val == N_dist(loc,2) && Nb(loc,2) ~= prev_comp
                %2nd component the same distance to comp.  Choose the
                %component with larger L-R overlap with comp.
                pn_loc = Nb(loc,2) - (start_comp - 1);
                pn_ovr = min(Pos(loc,3), Pos(pn_loc,3)) - ...
                          max(Pos(loc,1), Pos(pn_loc,1));
                pc_ovr = min(Pos(loc,3), Pos(prev_loc,3)) - ...
                          max(Pos(loc,1), Pos(prev_loc,1));
                if pc_ovr > pn_ovr
                    Nb(loc,2) = prev_comp;
                end
            end
            if val < N_dist(prev_loc,4)
                %comp is new closest bottom neighbour for prev_comp
                N_dist(prev_loc,4) = val;
                Nb(prev_loc,4) = comp;
            elseif val == N_dist(prev_loc,4) && Nb(prev_loc,4) ~= comp
                %2nd component the same distance to prev_comp.  Choose the
                %component with larger L-R overlap with prev_comp.
                pn_loc = Nb(prev_loc,4) - (start_comp - 1);
                pn_ovr = min(Pos(prev_loc,3), Pos(pn_loc,3)) - ...
                          max(Pos(prev_loc,1), Pos(pn_loc,1));
                c_ovr = min(Pos(prev_loc,3), Pos(loc,3)) - ...
                          max(Pos(prev_loc,1), Pos(loc,1));
                if c_ovr > pn_ovr
                    Nb(prev_loc,4) = comp;
                end
            end
        end
        %update component and row accordingly
        prev_comp = comp;
        prev_row  = RR(i);
        prev_loc  = loc;
        last_col  = CC(i);
    end
    %now get the left and right distances.  This requires working with the
    %transpose to ensure we get items in L-R order
    [CC, RR] = find(Comps{p}' ~= 0);
    prev_comp = 0;
    prev_loc  = 0;
    prev_col  = Inf;
    last_row  = NaN;
    for i=1:length(CC)
        comp = Comps{p}(RR(i), CC(i));
        loc = comp - (start_comp - 1);
        if last_row == RR(i) && prev_comp ~= comp
            %not first valid component encountered in this row 
            val = CC(i) - prev_col;
            if val < N_dist(loc,1)
                %new closest left neighbour for comp
                N_dist(loc,1) = val;
                Nb(loc,1) = prev_comp;
            elseif val == N_dist(loc,1) && Nb(loc,1) ~= prev_comp
                %2nd component the same distance to comp.  Choose the
                %component with larger T-B overlap with comp.
                pn_loc = Nb(loc,1) - (start_comp - 1);
                pn_ovr = min(Pos(loc,4), Pos(pn_loc,4)) - ...
                          max(Pos(loc,2), Pos(pn_loc,2));
                pc_ovr = min(Pos(loc,4), Pos(prev_loc,4)) - ...
                          max(Pos(loc,2), Pos(prev_loc,2));
                if pc_ovr > pn_ovr
                    Nb(loc,1) = prev_comp;
                end
            end
            if val < N_dist(prev_loc,3)
                %comp is new closest right neighbour for prev_comp
                N_dist(prev_loc,3) = val;
                Nb(prev_loc,3) = comp;
            elseif val == N_dist(prev_loc,3) && Nb(prev_loc,3) ~= comp
                %2nd component the same distance to prev_comp.  Choose the
                %component with larger T-B overlap with prev_comp.
                pn_loc = Nb(prev_loc,3) - (start_comp - 1);
                pn_ovr = min(Pos(prev_loc,4), Pos(pn_loc,4)) - ...
                          max(Pos(prev_loc,2), Pos(pn_loc,2));
                c_ovr = min(Pos(prev_loc,4), Pos(loc,4)) - ...
                          max(Pos(prev_loc,2), Pos(loc,2));
                if c_ovr > pn_ovr
                    Nb(prev_loc,3) = comp;
                end
            end
        end
        %update component and row accordingly
        prev_comp = comp;
        prev_col  = CC(i);
        prev_loc  = loc;
        last_row  = RR(i);
    end
    clear CC RR;
    fprintf('%.2fs: Neighbours of each component found\n', toc);
    num_clusts = size(Clust,1);
    for i=keep_list
        l = Pos(i,1); t = Pos(i,2); r = Pos(i,3); b = Pos(i,4);
        Clust = [Clust; struct('comp', start_comp - 1 + i, ...
          'pos', Pos(i,:), 'nb', Nb(i,:), 'pg', p, 'num', 1, ...
          'avg', zeros(b-t+1, r-l+1))]; 
        Clust(end).avg(find(Comps{p}(t:b,l:r) == Clust(end).comp(1))) = fg_val;

        if display_page || save_cc_page
            M(t,l:r,:) = repmat(out_col,1,r-l+1);
            M(t:b,l,:) = repmat(out_col,b-t+1,1);
            M(t:b,r,:) = repmat(out_col,b-t+1,1);
            M(b,l:r,:) = repmat(out_col,1,r-l+1);
        end
    end
    if display_page
        imshow(M);
        fprintf('Ready to start matching.  Press a key to begin\n');
        pause;
    end
    if save_cc_page
        fprintf('%.2fs: Writing components of page %d to disk\n', toc, p);
        imwrite(M, [img_prefix, num2str(p), '.', img_format], img_format);
    end
    fprintf('%.2fs: New components added to individual new clusters\n', toc);

    %perform refinements until the number of clusters remains constant
    refine_list = size(Clust,1):-1:num_clusts+1;
    while true
        num_clusts = size(Clust,1);

        %First look for straight matches among the clusters.  Initially this
        %should drastically reduce the number of clusters.  Just use simple
        %euclidian distance (to keep run-time down).
        if p ~= num_pgs || pass == 1
            fprintf('\n\n%.2fs: Starting straight-match refine pass\n', toc);
            [Clust, refine_list] = match_refine(Clust, refine_list, 'euc', ...
                                   match_thresh);
            if p == num_pgs
                pass = 2;
            end
        else
            %since all pages have been processed with straight matches
            %we now further refine using splits, merges, and matches using
            %the appropriate distance metric
            if pass > 2
                fprintf('\n\n%.2fs: Starting %s match refine pass\n', toc, ...
                        metric);
                [Clust, new_list] = match_refine(Clust, refine_list, metric, ...
                                    match_thresh);
                if length(new_list) ~= 0
                    refine_list = new_list;
                end
            end

            if pass == 2
                %first time through we want to attempt to merge over all
                %clusters
                refine_list = size(Clust,1):-1:1;
            end
            fprintf('\n\n%.2fs: Starting merge refine pass\n', toc);
            [Clust, Comps, new_list] = merge_refine(Clust, Comps, ...
                refine_list, metric, merge_thresh, merge_pct, merge_min_comps);

            if pass == 2
                %first time through we want to attempt to split over all
                %clusters
                refine_list = size(Clust,1):-1:1;
            elseif pass > 2 && length(new_list) ~= 0
                refine_list = new_list;
            end
            fprintf('\n%.2fs: Starting horizontal split refine pass\n', toc);
            [Clust, Comps, new_list] = split_refine(Clust, Comps, ...
                refine_list, metric, max_splits, split_thresh);

            if pass == 2
                refine_list = size(Clust,1):-1:1;
            elseif pass > 2 && length(new_list) ~= 0
                refine_list = new_list;
            end

            pass = pass + 1;
        end

        if num_clusts == size(Clust,1);
            fprintf('\n%.2fs: no further reduction in cluster size: %d\n', ...
                    toc, num_clusts);
            break;
        end
    end

    if display_page
        imshow(M);
        fprintf('Press any key to continue (to next page if one exists)\n');
        pause;
    end
end

%perform a final scale refine over the remaining elements
%fprintf('\n%.2fs: Starting scale refine pass\n', toc);
%[Clust, new_list] = scale_refine(Clust, size(Clust,1):-1:1, scale_thresh);

%finally, sort the clusters
Clust = sort_clusters(Clust);

fprintf('TOTAL NUMBER OF COMPONENTS FOUND: %d\n', tot_comps);
fprintf('TOTAL NUMBER OF CLUSTERS FOUND: %d\n', size(Clust,1));



% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
