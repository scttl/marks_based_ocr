function D = create_word_dictionary(Files)
% CREATE_DICTIONARY Parse text corpus creatng list and counts of words and chars
%
%   D = CREATE_DICTIONARY(Files)
%   Files should be a cell array (or single string), and give the full 
%   path and name of the file(s) to be used in creating the dictionary.
%
%   D is a struct with fields: word, word_count, char, char_count, char_bigram
%
%  NOTE: This requires matlab7 since it uses textscan, and assumes the text
%  has been preprocessed to remove all unwanted (punctuation/numeric/control)
%  characters
%   use something like:   tr -c -d '[:alnum:] .?!\n' < infile.txt > outfile.txt
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: create_word_dictionary.m,v 1.1 2006-06-10 21:01:40 scottl Exp $
%
% REVISION HISTORY
% $Log: create_word_dictionary.m,v $
% Revision 1.1  2006-06-10 21:01:40  scottl
% Initial revision.
%

% LOCAL VARS %
%%%%%%%%%%%%%%
num_files = size(Files,1);
D = {};
W = {};
C = [];


% CODE START %
%%%%%%%%%%%%%%
tic;
if num_files == 1
    fid = fopen(Files);
    if fid == -1
        error('unable to open file');
    end
    W = textscan(fid, '%s');
    frewind(fid);
    C = fread(fid);
    fclose(fid);
else
    for i = 1:num_files
        fid = fopen(Files{i});
        if fid == -1
            error('unable to open file');
        end
        W = [W; textscan(fid, '%s')];
        frewind(fid);
        C = [C; fread(fid)];
        fclose(fid);
    end
end
fprintf('%.2fs: finished parsing files\n', toc);

[D.word,i,i] = unique(W{1}(:));
D.word_count = accumarray(i,1);
fprintf('%.2fs: finished creating and counting words\n', toc);

[D.char,i,i] = unique(C);
D.char = char(D.char);
D.char_count = accumarray(i,1);
fprintf('%.2fs: finished creating and counting chars\n', toc);

num_chars = length(D.char);
D.char_bigram = zeros(num_chars);
pos = zeros(1, max(single(D.char)));
for i = 1:num_chars
    pos(single(D.char(i))) = i;
end

for i = 2:length(C)
    D.char_bigram(pos(C(i-1)), pos(C(i))) = ...
    D.char_bigram(pos(C(i-1)), pos(C(i))) + 1;
end
fprintf('%.2fs: finished creating character bigram matrix\n', toc);
