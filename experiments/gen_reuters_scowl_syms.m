%this script will generate the Syms struct (symbol counts and information) 
%taken from the Reuters text corpus.
%
%This struct forms the basis of the output symbols that can be generated, and
%will be used to determine a mapping from clusters of ink.

global MOCR_PATH;  %used to determine where to save results

%set the following line to true to record process
create_diary = true;
if create_diary
    diary_file = [MOCR_PATH, '/results/gen_reuters_scowl_syms.diary'];
    if exist(diary_file)
        delete(diary_file);
    end
    diary(diary_file);
    diary on;
    fprintf('EXPERIMENT STARTED: %s\n', datestr(now));
end

run_syms=true;
use_separate_words=true;

%this file should point at a file containing the list of possible symbols
sym_in_file = [MOCR_PATH, '/data/input_utf8_syms.txt'];

%this should give the path to where results will be kept
res_symfile = [MOCR_PATH, '/data/reuters_scowl_syms.mat'];
%res_symfile = [MOCR_PATH, '/results/reuters_scowl_syms.mat'];

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

%we use a separate corpus for word lists
word_files = {[MOCR_PATH, '/data/scowl/english-words.10'];
              [MOCR_PATH, '/data/scowl/english-words.20'];
              [MOCR_PATH, '/data/scowl/english-upper.10'];
              [MOCR_PATH, '/data/scowl/english-contractions.10'];
              [MOCR_PATH, '/data/scowl/english-contractions.35'];
              [MOCR_PATH, '/data/scowl/american-words.10'];
              [MOCR_PATH, '/data/scowl/american-words.20'];
              [MOCR_PATH, '/data/scowl/british-words.10'];
              [MOCR_PATH, '/data/scowl/british-words.20']};

tic;
%get the symbols
if run_syms
    Syms = create_alphabet(sym_in_file, 'corpora_files', corpus_files, ...
           'use_srilm', false);
    save(res_symfile, 'Syms');
    fprintf('symbols complete: %f\n', toc);
end

%update the word list
if use_separate_words
    D = create_word_dictionary(word_files);
    Syms.words = D.word;
    Syms.word_count = D.word_count;
    save(res_symfile, 'Syms');
    fprintf('separate wordlist complete: %f\n', toc);
end

if create_diary
    fprintf('EXPERIMENT ENDED: %s\n', datestr(now));
    diary off;
end
