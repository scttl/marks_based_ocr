%this script will generate the Syms struct (symbol counts and information) 
%taken from the Reuters text corpus, with positional counts taken from
%classifying 1/3 of the ground truth pages of the UNLV 'B' dataset.
%
%This struct forms the basis of the output symbols that can be generated, and
%will be used to determine a mapping from clusters of ink.

global MOCR_PATH;  %used to determine where to save results

%set the following line to true to record process
create_diary = true;
if create_diary
    diary_file = [MOCR_PATH, '/results/gen_reuters_BGT_syms.diary'];
    if exist(diary_file)
        delete(diary_file);
    end
    diary(diary_file);
    diary on;
    fprintf('EXPERIMENT STARTED: %s\n', datestr(now));
end

run_syms=true;
run_text_pos=true;

%this file should point at a file containing the list of possible symbols
sym_in_file = [MOCR_PATH, '/data/input_utf8_syms.txt'];

%this should give the path to where results will be kept
res_symfile = [MOCR_PATH, '/data/reuters_BGT_syms.mat'];

%these point at the files containing the actual text of the corpus
corpus_files = {[MOCR_PATH, '/data/reuters/000.body.txt']};

%this should be replaced with directory containing ground truth text to be
%clustered for positional count statistics.
pos_count_dir = [MOCR_PATH, '/data/unlv_ocr/B/B_GT/'];

%what portion of the files in the directory should be taken
pos_proportion = .33;

%how many words are used in positional counts
max_word_len = 15;

tic;
%get the symbols
if run_syms
    Syms = create_alphabet(sym_in_file, 'corpora_files', corpus_files, ...
           'use_srilm', false);
    save(res_symfile, 'Syms');
    fprintf('symbols complete: %f\n', toc);
end

%update the positional counts
if run_text_pos
    xx = dir(pos_count_dir);
    imgs = cell(length(xx),1);
    [imgs{:}] = deal(xx.name);
    imgs = imgs(3:end);  %remove . and ..
    docs = unique(regexprep(imgs, '(\w*)\_.*', '$1'));
    num_docs = floor(.33 * length(docs));
    Syms.num = 0;
    Syms.val = {};
    Syms.count = [];
    Syms.pos_count = cell(1,max_word_len);
    for ii=1:num_docs
        idx = strmatch(docs{ii}, imgs);
        Files = regexprep(imgs(idx), '(.*)', [pos_count_dir, '$1']);
        [Clust, Comps, Lines] = create_text_clusters(Files, 'max_word_len', ...
                                max_word_len);
        Syms.num = Syms.num + Clust.num;
        Syms.val = [Syms.val(:); Clust.truth_label'];
        Syms.count = [Syms.count; Clust.num_comps];
        for jj=1:max_word_len
            Syms.pos_count{jj} = [Syms.pos_count{jj}; Clust.pos_count{jj}];
        end
    end

    save(res_symfile, 'Syms');
    fprintf('separate positional counts complete: %f\n', toc);
end

if create_diary
    fprintf('EXPERIMENT ENDED: %s\n', datestr(now));
    diary off;
end
