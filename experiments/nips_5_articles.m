%this script will load 5 articles (each with ~ 8 pages of images) from
%the NIPS 2001 training set, then cluster the blobs of ink found on each page,
%determine the lines of the page, and create training data out of each of them
%(both the lines, and the cluster counts for the "language model")

global MOCR_PATH;  %used to determine where to save results

%set the following line to true to record process
create_diary = true;
if create_diary
    diary_file = [MOCR_PATH, '/results/nips_5_articles.diary'];
    if exist(diary_file)
        delete(diary_file);
    end
    diary(diary_file);
    diary on;
    fprintf('EXPERIMENT STARTED: %s\n', datestr(now));
end

run_comps=true;
run_lines=true;
run_cluster=true;
run_dictionary=true;
run_map=true;
run_ocr_analysis=true;

%if attempting to determine mappings, this should list the file containing the
%Syms corpora struct
syms_struct_file = [MOCR_PATH, '/data/reuters_syms.mat'];

Files = ls('~scottl/research/TRAINING_DATA/nips_2001/AA0[1-7]*.tif');
res_datafile = [MOCR_PATH, '/results/nips_5_articles.mat'];
res_txtfile = [MOCR_PATH, '/results/nips_5_articles.txt'];


%must convert Files from a single long newline delimited string into a cell 
%array
idx = strfind(Files, char(10));
for ii=1:length(idx)
    if ii == 1
        len(ii) = idx(ii);
    else
        len(ii) = idx(ii) - idx(ii-1);
    end
end
Files = strtrim(mat2cell(Files, 1, len)');

%ensure that if we are doing mappings, that the struct file can be loaded
if run_map
    load(syms_struct_file);
end

%get the components
tic;
if run_comps
    Comps = get_comps(Files, 'min_elem_width',2, 'min_elem_height',2, ...
            'max_elem_width',100, 'max_elem_height',100);
    fprintf('components complete: %f\n', toc);
    save(res_datafile, 'Comps');
end

load(res_datafile);

%determine line boundaries
if run_lines
    [Lines, Comps] = get_lines(Comps, 'base_thresh', .2, 'xheight_thresh', .3);
    fprintf('lines complete: %f\n', toc);
    Comps = merge_diacritic_comps(Comps, 'max_dist', 7);
    fprintf('diacritic component merger complete: %f\n', toc);
    save(res_datafile, 'Comps', 'Lines');
end

load(res_datafile);

%now cluster the components
if run_cluster
    [Clust, Comps] = cluster_comps(Comps, 'Lines',Lines, ...
      'refine_clusters',true, 'match_metric','hausdorff', ...
      'split_metric','euc', 'straight_match_thresh',.009, ...
      'match_thresh', 1.5, 'split_thresh',.013, 'max_splits',2, ...
      'merge_thresh',3, 'merge_pct',.85, 'merge_min_comps',3, ...
      'avg_splits',false, 'avg_matches',true, 'resize_imgs',false, ...
      'resize_method','nearest', 'use_thinned_imgs',false);
    [Clust, Comps] = add_space_model(Clust, Comps, 'space_width', [], ...
                     'space_height', []);
    [Clust, Comps] = sort_clusters(Clust, Comps);
    save(res_datafile, 'Clust', 'Comps', 'Lines');
    fprintf('clustering complete: %f\n', toc);
end

load(res_datafile);

%now create a language model from the same dataset
if run_dictionary
    [Clust, Comps] = create_cluster_dictionary(Clust, Comps);

    fprintf('dictionary complete: %f\n', toc);
    save(res_datafile, 'Clust', 'Comps', 'Lines');
end

load(res_datafile);

%now attempt to infer mappings based on positional information
if run_map
    [order, score] = positional_learn_mappings(Clust, Syms, ...
                     'dist_metric', 'manhattan');
    fprintf('mapping complete: %f\n', toc);
    save(res_datafile, 'Clust', 'Comps', 'Lines', 'order', 'score');
end

load(res_datafile);

%print out mapped ground truth to a text file, and determine accuracy stats
if run_ocr_analysis
    print_ocr_text(1:Lines.num, Comps, Syms, order(:,1), ...
                   'display_text', false, 'save_results', true, ...
                   'save_file', res_txtfile);
    %can't analyze these results since no ground truth text available
    fprintf('ocr analysis complete: %f\n', toc);
end

if create_diary
    fprintf('EXPERIMENT ENDED: %s\n', datestr(now));
    diary off;
end
