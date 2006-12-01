function [Clust, Comps] = cluster_comps(Comps, varargin)
% CLUSTER_COMPS   Cluster together connected component blobs
%
%   [Clust, Comps] = CLUSTER_COMPS(COMPS, [VAR1, VAL1]...)
%   COMPS is a struct like that returned from get_comps().  See get_comps.m for
%   field details.
%
%   VAR1 and VAL1 are optional, and can be used to override the default values
%   for the LOCAL VARS defined below.  Each VAR1 should be a string giving the
%   exact name of the variable to override.  Each VAL1 should be the updated
%   value to use for that variable.
%
%   Clust is a struct meant to hold information on all clusters found while
%   processing FILES.  It contains the following fields:
%     num   - a scalar denoting the number of clusters in the struct.
%     num_comps - a vector giving the number of components belonging to each
%                 cluster.
%     mode_num - a vector giving the value for the mode of the cluster.  Used
%                in re-averaging where the most frequent cluster (based on this
%                mode_num) is used instead of taking the mean of the clusters.
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
%     found_offsets - this is a boolean that will be set to true if offset 
%                     values for each cluster have been calculated and false
%                     otherwise.
%     descender_off - this holds the most frequently occuring relative position 
%                     of the bottom edge of the components from the baseline.  
%                     Positive values refer to the number of pixels below the 
%                     baseline, the bottom of the image lies, and negative 
%                     refers to pixels above the baseline.
%     ascender_off - this holds the most frequently occuring relative position 
%                    of the top edge of the components from the x-height.  
%                     Positive values refer to the number of pixels below the 
%                     x-height, the top of the image lies, and negative 
%                     refers to pixels above the x-height.
%     num_trans - a scalar that determines the number of left-right component
%                 transitions found.  Note that this is not calculated here
%     bigram - a num x num matrix that denotes (smoothed) transition 
%              probabilities.  Again, this isn't calculated here.
%
%     found_true_labels - this is a boolean that will be set to true if the
%                         ground truth labels for each cluster have been
%                         caluculated.  This can be done with a call to
%                         clust_ground_truth_label()
%     truth_label - this cell array will hold strings containing the true
%                   symbols this cluster blob of ink refers to.
%     model_spaces - this boolean will be set to true if spaces have been
%                    modelled and counts taken.  See add_space_model()


% CVS INFO %
%%%%%%%%%%%%
% $Id: cluster_comps.m,v 1.11 2006-12-01 22:57:37 scottl Exp $
%
% REVISION HISTORY
% $Log: cluster_comps.m,v $
% Revision 1.11  2006-12-01 22:57:37  scottl
% added erosion ability, updated some documentation, removed spurious lines
%
% Revision 1.10  2006-11-13 17:56:32  scottl
% small spacing improvements.
%
% Revision 1.9  2006-10-29 17:24:53  scottl
% change to cluster struct, to use descender and ascender offsets, instead
% of a single offset field.
%
% Revision 1.8  2006/10/18 15:51:16  scottl
% Excised Component related work, and moved to a different block of code.
% Updates to better process optional arguments.
%
% Revision 1.7  2006/08/30 17:41:08  scottl
% small fix to not attempt to process files which don't have a corresponding
% jtag file.
%
% Revision 1.6  2006/08/24 21:40:04  scottl
% added ability to use the mode instead of taking the average of cluster
% intensities while refining.
%
% Revision 1.5  2006/08/14 01:40:53  scottl
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


% LOCAL VARS %
%%%%%%%%%%%%%%

%the number of directions to consider when looking for conn. components (4 or
%8) are the only valid entries
num_dirs = 8;

%these parameters control the distance metric used during grouping.  Valid
%choices are: 'euc', 'conv', or 'hausdorff' for Euclidian, Convolutional
%Euclidian, or Hausdorff (with underlying Euclidian) distance measurements
%used.
match_metric='euc';  %'hausdorff';
split_metric='euc';

%the following parameters control when things are grouped together.
straight_match_thresh = .009;  %good for NIPS papers
%straight_match_thresh = .013;  %infinity book setting
match_thresh = .009;  %1.3;  %hausdorff euclidian distance thresh

split_thresh = .013;  %good euc thresh
max_splits = 2;

merge_thresh = 3;  %look for compoonents no more than 3 pixels apart
merge_pct = .85;  %ensure at least 85% of the elements in each cluster match
merge_min_comps = 3; %ensure that the cluster has at least 3 elements 

%these parameters determine whether averaging of combined components will be
%done (or whether the modal image is used instead -- recommended for Hausdorff)
avg_splits = false;
avg_matches = true;

%by default we don't rescale components, unless the boolean below is set to
%true.  This also requires that Comps.scale_factor exists (and thus that
%get_lines() has been run)
resize_imgs = false;
%which method should be used to resize components (see imresize() for choices)
resize_method = 'nearest';  

%by default we don't work with thinned representations of the images
use_thinned_imgs = false;

%by default we don't erode the images (in an effort to encourage splitting
%glued-together components)
erode_imgs = false;
erosion_se = strel('square', 3);  %try help strel for other choices
eroision_iters = 1;

%what value is assigned to 'on' pixels in the cluster averages initially
fg_val = 1;

%by defeault turn-off the override warning display since it will be generated
%each time we call match_refine etc.
override_display = 'OFF';   %other legal value is 'ON'


% CODE START %
%%%%%%%%%%%%%%
tic;

%process input arguments
if nargin < 1
    error('Comps struct must be specified as the first argument!');
elseif nargin > 1
    process_optional_args(varargin{:});
end
warning(override_display, 'MBOCR:override');

%initialize the Clust struct
Clust = init_clust();

%iterate through each page adding its components to our cluster list and
%performing straight match refinements to reduce the number of clusters
for pp=1:length(Comps.files)

    fprintf('Processing page %d\n', pp);
    comp_idcs = find(Comps.pg == pp)';
    num_pg_comps = length(comp_idcs);

    if num_pg_comps == 0
        %skip to the next page
        continue;
    end

    M = ~imread(Comps.files{pp});
    Lbl_img = bwlabel(M, num_dirs);

    clust_idcs = Clust.num+1:Clust.num+num_pg_comps;
    Comps.clust(comp_idcs,1) = clust_idcs;

    %initially add each component on this page as a new cluster
    Clust.num = Clust.num + num_pg_comps;
    Clust.num_comps(clust_idcs,1) = 1;
    Clust.mode_num(clust_idcs,1) = 1;
    Clust.comps(clust_idcs,1) = num2cell(comp_idcs');
    Clust.avg(clust_idcs,1) = cell(num_pg_comps,1);
    Clust.norm_sq(clust_idcs,1) = 0;
    Clust.refined(clust_idcs,1) = false;
    Clust.changed(clust_idcs,1) = false;
    if Comps.found_lines
        Clust.found_offsets = true;
        Clust.descender_off(clust_idcs,1) = Comps.descender_off(comp_idcs);
        Clust.ascender_off(clust_idcs,1) = Comps.ascender_off(comp_idcs);
    end

    %update the avg and norm_sq fields
    [Clust.avg(clust_idcs), Clust.norm_sq(clust_idcs)] = get_comp_imgs(...
                                                         Comps, comp_idcs);

    %perform a single straight match refinement over these new clusters to
    %reduce their number.
    fprintf('%.2fs: Starting straight-match refine pass\n', toc);
    [Clust, Comps] = match_refine(Clust, Comps, 'dist_metric', 'euc', ...
                     'distance_thrsh', straight_match_thresh);
    fprintf('\n');
end

if resize_imgs && Comps.found_lines
    %try rescaling the images
    fprintf('\n\n%.2fs: Rescaling odd-sized clusters\n', toc);

    for ii=1:Clust.num
        co = Clust.comps{ii};
        sc = mode(Comps.scale_factor(co));
        if sc ~= 1
            Clust.avg{ii} = imresize(Clust.avg{ii}, sc, resize_method);
            Clust.norm_sq(ii) = sum(Clust.avg{ii}(:));
        end
    end
end

if use_thinned_imgs
    %try thinning the images
    fprintf('\n\n%.2fs: thinning Cluster averages\n', toc);
    for ii=1:Clust.num
        Clust.avg{ii} = double(bwmorph(Clust.avg{ii}, 'thin', Inf));
        Clust.norm_sq(ii) = sum(Clust.avg{ii}(:));
    end
elseif erode_imgs
    fprintf('\n\n%.2fs: eroding Cluster averages\n', toc);
    for ii=1:erosion_iters
        for jj=1:Clust.num
            Clust.avg{ii} = imerode(Clust.avg{ii}, erosion_se);
            Clust.norm_sq(ii) = sum(Clust.avg{ii}(:));
        end
    end
end

%now repeatedly perform merges, splits, and matches until the number of
%clusters remains constant.  First time through, refine all clusters.
Clust.refined(:) = false;
Clust.changed(:) = false;
first_pass = true;
while true
    num_clusts = Clust.num;

    if first_pass
        %first time through we want to attempt to merge over all
        %clusters
        Clust.refined(:) = false;
    end
    fprintf('\n\n%.2fs: Starting merge refine pass\n', toc);
    [Clust, Comps] = merge_refine(Clust, Comps,  'max_dist', merge_thresh, ...
                     'valid_pct', merge_pct, 'min_match_comp', ...
                     merge_min_comps, 'resize_imgs', resize_imgs, ...
                     'resize_method', resize_method, 'use_thinned_imgs', ...
                     use_thinned_imgs);

    if first_pass
        %first time through we want to attempt to match over all
        %clusters
        Clust.refined(:) = false;
    end
    fprintf('\n\n%.2fs: Starting %s match refine pass\n', toc, ...
            match_metric);
    [Clust, Comps] = match_refine(Clust, Comps, 'dist_metric', match_metric, ...
                     'distance_thresh', match_thresh, 'use_avg', ...
                     avg_matches);

    if first_pass
        %first time through we want to attempt to split over all
        %clusters
        Clust.refined(:) = false;
    end
    fprintf('\n%.2fs: Starting horizontal split refine pass\n', toc);
    [Clust, Comps] = split_refine(Clust, Comps, 'dist_metric', split_metric, ...
                     'max_splits', max_splits, 'dist_thresh', split_thresh, ...
                     'use_avg', avg_splits, 'resize_imgs', resize_imgs, ...
                     'resize_method', resize_method, 'use_thinned_imgs', ...
                     use_thinned_imgs);


    %to speed up subsequent refinements we only need to further process
    %those Clusters that have changed during the last round of
    %refinements
    Clust.refined(Clust.changed == true) = false;
    Clust.changed(:) = false;
    first_pass = false;

    if num_clusts == Clust.num;
        fprintf('\n%.2fs: no further reduction in cluster size: %d\n', ...
                toc, Clust.num);
        break;
    end
end

%finally, sort the clusters
[Clust, Comps] = sort_clusters(Clust, Comps);

%if warnings have been turned off, turn them back on
if strcmp(override_display, 'OFF')
    warning('ON', 'MBOCR:override');
end

fprintf('TOTAL NUMBER OF COMPONENTS NOW: %d\n', Comps.max_comp);
fprintf('TOTAL NUMBER OF CLUSTERS FOUND: %d\n', Clust.num);



% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize an empty cluster struct
function Clust = init_clust()
Clust.num = 0;
Clust.num_comps = [];
Clust.mode_num = [];
Clust.comps = {};
Clust.avg = {};
Clust.norm_sq = [];
Clust.refined = logical([]);
Clust.changed = logical([]);
Clust.found_offsets = false;
Clust.descender_off = int16([]);
Clust.ascender_off = int16([]);
Clust.num_trans = 0;
Clust.bigram = [];
Clust.found_true_labels = false;
Clust.true_labels = {};
Clust.model_spaces = false;
