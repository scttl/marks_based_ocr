%this script will generate the Syms struct (symbol counts and information) 
%taken from the Reuters text corpus.
%
%This struct forms the basis of the output symbols that can be generated, and
%will be used to determine a mapping from clusters of ink.

global MOCR_PATH;  %used to determine where to save results

%set the following line to true to record process
create_diary = true;
if create_diary
    diary_file = [MOCR_PATH, '/results/gen_reuters_pos_15_syms.diary'];
    if exist(diary_file)
        delete(diary_file);
    end
    diary(diary_file);
    diary on;
    fprintf('EXPERIMENT STARTED: %s\n', datestr(now));
end

run_syms=true;

%this file should point at a file containing the list of possible symbols
sym_in_file = [MOCR_PATH, '/data/input_utf8_syms.txt'];

%this should give the path to where results will be kept
res_symfile = [MOCR_PATH, '/data/reuters_pos_15_syms.mat'];

%these point at the files containing the actual text of the corpus
corpus_files = {[MOCR_PATH, '/data/reuters/000.body.txt']};
%corpus_files = {[MOCR_PATH, '/data/reuters/000.body.txt'];
%                [MOCR_PATH, '/data/reuters/001.body.txt'];
%                [MOCR_PATH, '/data/reuters/002.body.txt'];
%                [MOCR_PATH, '/data/reuters/003.body.txt'];
%                [MOCR_PATH, '/data/reuters/004.body.txt'];
%                [MOCR_PATH, '/data/reuters/005.body.txt'];
%                [MOCR_PATH, '/data/reuters/006.body.txt'];
%                [MOCR_PATH, '/data/reuters/007.body.txt'];
%                [MOCR_PATH, '/data/reuters/008.body.txt'];
%                [MOCR_PATH, '/data/reuters/009.body.txt'];
%                [MOCR_PATH, '/data/reuters/010.body.txt'];
%                [MOCR_PATH, '/data/reuters/011.body.txt'];
%                [MOCR_PATH, '/data/reuters/012.body.txt'];
%                [MOCR_PATH, '/data/reuters/013.body.txt'];
%                [MOCR_PATH, '/data/reuters/014.body.txt'];
%                [MOCR_PATH, '/data/reuters/015.body.txt'];
%                [MOCR_PATH, '/data/reuters/016.body.txt'];
%                [MOCR_PATH, '/data/reuters/017.body.txt'];
%                [MOCR_PATH, '/data/reuters/018.body.txt'];
%                [MOCR_PATH, '/data/reuters/019.body.txt'];
%                [MOCR_PATH, '/data/reuters/020.body.txt'];
%                [MOCR_PATH, '/data/reuters/021.body.txt']};

tic;
%get the symbols
if run_syms
    Syms = create_alphabet(sym_in_file, 'corpora_files', corpus_files, ...
           'use_srilm', false, 'max_word_len', 15, ...
           'template_font', 'USE_TEX', 'font_ptsize', 36);
    %must fixup the normalization term
    for ii=1:length(Syms.pos_count)
        val = Syms.pos_count{ii} .* Syms.pos_total;
        norms = sum(val,2);
        norms(norms == 0) = 1;  %to prevent dividing by 0
        Syms.pos_norms{ii} = repmat(norms, 1, size(val,2));
        Syms.pos_count{ii} = val ./ Syms.pos_norms{ii};
    end
    fprintf('symbols complete: %f\n', toc);
    save(res_symfile, 'Syms');
end

if create_diary
    fprintf('EXPERIMENT ENDED: %s\n', datestr(now));
    diary off;
end
