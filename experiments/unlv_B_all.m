%this script will attempt to cluster and save the results from the UNLV ISRI
%OCR dataset (the 'B' dataset in particular)

global MOCR_PATH;  %used to determine where to save results

%set the following line to true to record process
create_diary = true;
if create_diary
    diary_file = [MOCR_PATH, '/results/unlv_B_all.diary'];
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
run_pos_map=true;
run_word_map=true;
run_ocr_analysis=true;

%if attempting to determine mappings, this should list the file containing the
%Syms corpora struct
syms_struct_file = [MOCR_PATH, '/data/reuters_syms.mat'];

%this file should point at a file containing the list of pages to run
pg_file = [MOCR_PATH, '/data/unlv_ocr/B/PAGES'];

%these determine where to find the pages listed above
pg_prefix = [MOCR_PATH, '/data/unlv_ocr/B/B_B/'];
pg_suffix = '.FF';  %use fine-mode fax to compare with Nagy paper

%these determine where to find the ground truth text files for OCR analysis
gt_prefix = [MOCR_PATH, '/data/unlv_ocr/B/B_GT/'];

%this should give the path to the base part of where results will be kept
res_base = [MOCR_PATH, '/results/B'];
if ~exist(res_base, 'dir')
    [s,w] = unix(['mkdir -p ', res_base]);
    if s~=0
        error('problem creating dir: %s', res_base);
    end
end

%open and read the list of pages from the pg_file
imgs = textread(pg_file, '%s');
docs = unique(regexprep(imgs, '(\w*)\-.*', '$1'));
num_docs = length(docs);
%convert any  '-' in the listed files to '_'
imgs = regexprep(imgs, '-', '_');

%ensure that if we are doing mappings or analysis, the struct file can be loaded
if run_pos_map || run_word_map || run_ocr_analysis
    load(syms_struct_file);
end

tic;
for ii=1:num_docs
    fprintf('%.2f: Processing document: %s\n', toc, docs{ii});
    res_dir = [res_base, '/', docs{ii}];
    if ~exist(res_dir, 'dir')
        fprintf('%.2f: Creating new dir\n', toc);
        [s,w] = unix(['mkdir -p ', res_dir]);
        if s~=0
            warning('MBOCR:NoDir', 'problem creating dir: %s\n', res_dir);
            continue;
        end
    end

    res_datafile = [res_dir, '/data.mat'];
    idx = strmatch(docs{ii}, imgs);
    Files = regexprep(imgs(idx), '(.*)', [pg_prefix, '$1', pg_suffix]);

    %get the components
    if run_comps
        Comps = get_comps(Files, 'min_elem_width',2, 'min_elem_height',2, ...
                'max_elem_width',100, 'max_elem_height',100);
        fprintf('components complete: %f\n', toc);
        save(res_datafile, 'Comps');
    end

    load(res_datafile);

    %determine line boundaries
    if run_lines
        [Lines, Comps] = get_lines(Comps, 'base_thresh', .2, ...
                         'xheight_thresh', .3);
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
      'split_metric','euc', 'straight_match_thresh', .01, ...
      'match_thresh', 1.5, 'split_thresh',.015, 'max_splits',2, ...
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
    if run_pos_map
        [order, score] = positional_learn_mappings(Clust, Syms, ...
                         'dist_metric', 'manhattan');
        fprintf('position based ordering complete: %f\n', toc);
        save(res_datafile, 'Clust', 'Comps', 'Lines', 'order', 'score');
    end

    load(res_datafile);

    if run_word_map
        map = word_lookup_map(Clust, Comps, Syms, 'order', order);
        fprintf('word lookup mapping complete: %f\n', toc);
        save(res_datafile, 'Clust', 'Comps', 'Lines', 'order', 'score', 'map');
    end
    
    load(res_datafile);

    %print out mapped ground truth to a text file, and determine accuracy stats
    if run_ocr_analysis
        pgs = unique(Comps.pg);
        for jj=1:length(pgs)
            lines = find(Lines.pg == pgs(jj));
            rgns = find(Comps.regions(:,1) == pgs(jj));
            this_pg = imgs{idx(jj)};
            for kk=1:length(rgns)
                txt_file = [this_pg, '.Z', sprintf('%02d', ...
                            Comps.regions(rgns(kk),2))];
                res_txtfile = [res_dir, '/', txt_file];
                res_char_rprtfile = [res_txtfile, '.char_rprt'];
                res_word_rprtfile = [res_txtfile, '.word_rprt'];
                print_ocr_text(lines, Comps, Syms, map, ...
                       'display_text', false, 'save_results', true, ...
                       'keep_region', Comps.regions(rgns(kk),3:6), ...
                       'save_file', res_txtfile);
                cmd = ['accuracy ', gt_prefix, txt_file, ' ', ...
                       res_txtfile, ' ', res_char_rprtfile];
                s = unix(cmd);
                if s ~= 0
                    error('prob running accuracy. cmd: %s', cmd);
                end
                cmd = ['wordacc ', gt_prefix, txt_file, ' ', ...
                       res_txtfile, ' ', res_word_rprtfile];
                s = unix(cmd);
                if s ~= 0
                    error('prob running word accuracy. cmd: %s', cmd);
                end
            end
        end
        %combine all the report files in this directory into a single
        %cumulative report
        char_rprts = dir([res_dir, '/*.char_rprt']);
        word_rprts = dir([res_dir, '/*.word_rprt']);
        char_rprt_list = '';
        word_rprt_list = '';
        if length(char_rprts) > 1
            for jj=1:length(char_rprts)
                char_rprt_list = [char_rprt_list, res_dir, '/', ...
                                  char_rprts(jj).name, ' '];
                word_rprt_list = [word_rprt_list, res_dir, '/', ...
                                  word_rprts(jj).name, ' '];
            end
            cmd = ['accsum ', char_rprt_list, ' > ', res_dir, '/', ...
                   docs{ii}, '.chartot_rprt'];
            s = unix(cmd);
            if s ~= 0
                error('prob running accsum. cmd: %s', cmd);
            end
            cmd = ['wordaccsum ', word_rprt_list, ' > ', res_dir, '/', ...
                   docs{ii}, '.wordtot_rprt'];
            s = unix(cmd);
            if s ~= 0
                error('prob running wordaccsum. cmd: %s', cmd);
            end
        else
            %just copy the single file for the total count
            cmd = ['cp ', res_dir, '/', char_rprts(1).name, ' ', res_dir, ...
                   '/', docs{ii}, '.chartot_rprt'];
            s = unix(cmd);
            if s ~= 0
                error('prob running cp. cmd: %s', cmd);
            end
            cmd = ['cp ', res_dir, '/', word_rprts(1).name, ' ', res_dir, ...
                   '/', docs{ii}, '.wordtot_rprt'];
            s = unix(cmd);
            if s ~= 0
                error('prob running cp. cmd: %s', cmd);
            end
        end

        fprintf('ocr analysis complete: %f\n', toc);
    end

end

if create_diary
    fprintf('EXPERIMENT ENDED: %s\n', datestr(now));
    diary off;
end
