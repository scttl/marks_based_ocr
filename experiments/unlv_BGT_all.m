%this script will attempt to use the ground truth text from the UNLV ISRI OCR
%'B' dataset, to create clusters, then attempt to infer the labels, and save
%the accuracy results.

global MOCR_PATH;  %used to determine where to save results

%set the following line to true to record process
create_diary = true;
if create_diary
    diary_file = [MOCR_PATH, '/results/unlv_BGT_all.diary'];
    if exist(diary_file)
        delete(diary_file);
    end
    diary(diary_file);
    diary on;
    fprintf('EXPERIMENT STARTED: %s\n', datestr(now));
end

run_cluster=true;
run_pos_map=true;
run_vote_map=false;
run_word_map=true;
run_ocr_analysis=true;

%if attempting to determine mappings, this should list the file containing the
%Syms corpora struct
syms_struct_file = [MOCR_PATH, '/data/reuters_pos_15_syms.mat'];

%this directory determines where to find the ASCII text pages
pg_dir = [MOCR_PATH, '/data/unlv_ocr/B/B_GT/'];

%this should give the path to the base part of where results will be kept
res_base = [MOCR_PATH, '/results/BGT'];
if ~exist(res_base, 'dir')
    [s,w] = unix(['mkdir -p ', res_base]);
    if s~=0
        error('problem creating dir: %s', res_base);
    end
end

%open and read the list of pages from the pg_file
xx = dir(pg_dir);
imgs = cell(length(xx),1);
[imgs{:}] = deal(xx.name);
imgs = imgs(3:end);  %remove . and ..
docs = unique(regexprep(imgs, '(\w*)\_.*', '$1'));
num_docs = length(docs);

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
    Files = regexprep(imgs(idx), '(.*)', [pg_dir, '$1']);

    %now cluster the components
    if run_cluster
        [Clust, Comps, Lines] = create_text_clusters(Files);
        save(res_datafile, 'Clust', 'Comps', 'Lines');
        fprintf('clustering complete: %f\n', toc);
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
    if run_vote_map
        [order, score] = vote_learn_mappings(Clust, Comps, Syms, ...
                         'limit_to_map', true, 'word_count_weight_pct', 1, ...
                         'add_first_pos_up_let', true, ...
                         'add_last_pos_punct_syms', true, ...
                         'valid_punct_syms', '.,:;!?');
        save(res_datafile, 'Clust', 'Comps', 'Lines', 'order', 'score');
    end

    load(res_datafile);

    if run_word_map
        map = word_lookup_map(Clust, Comps, Syms, 'order', order);
        fprintf('word lookup mapping complete: %f\n', toc);
        save(res_datafile, 'Clust', 'Comps', 'Lines', 'order', 'score', 'map');
    else
        %convert the order to a map
        map = cell2mat(order);
        map = map(:,1);
        save(res_datafile, 'Clust', 'Comps', 'Lines', 'order', 'score', 'map');
    end
    
    load(res_datafile);

    %print out mapped ground truth to a text file, and determine accuracy stats
    if run_ocr_analysis
        pgs = Comps.files;
        for jj=1:length(pgs)
            lines = find(Lines.pg == jj);
            txt_file = imgs{idx(jj)};
            res_txtfile = [res_dir, '/', txt_file];
            res_char_rprtfile = [res_txtfile, '.char_rprt'];
            res_word_rprtfile = [res_txtfile, '.word_rprt'];
            print_ocr_text(lines, Comps, Syms, map, ...
                   'display_text', false, 'save_results', true, ...
                   'save_file', res_txtfile);
            cmd = ['accuracy ', pg_dir, txt_file, ' ', ...
                   res_txtfile, ' ', res_char_rprtfile];
            s = unix(cmd);
            if s ~= 0
                error('prob running accuracy. cmd: %s', cmd);
            end
            cmd = ['wordacc ', pg_dir, txt_file, ' ', ...
                   res_txtfile, ' ', res_word_rprtfile];
            s = unix(cmd);
            if s ~= 0
                error('prob running word accuracy. cmd: %s', cmd);
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
