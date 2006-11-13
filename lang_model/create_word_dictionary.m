function D = create_word_dictionary(Files, varargin)
% CREATE_WORD_DICTIONARY Parse text corpus creatng list and counts of words and 
% chars
%
%   D = CREATE_WORD_DICTIONARY(Files, [var1, val1]...)
%   Files should be a cell array (or single string), and give the full 
%   path and name of the file(s) to be used in creating the dictionary.
%
%   D is a struct with fields: word, word_count, char, char_count, char_bigram,
%   char_trigram, first_count, and first_bg_count
%
%  NOTE: This requires matlab7 since it uses textscan, and assumes the text
%  has been preprocessed to remove all unwanted (punctuation/numeric/control)
%  characters.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: create_word_dictionary.m,v 1.4 2006-11-13 17:58:03 scottl Exp $
%
% REVISION HISTORY
% $Log: create_word_dictionary.m,v $
% Revision 1.4  2006-11-13 17:58:03  scottl
% added trigram counts.
%
% Revision 1.3  2006-10-29 17:18:45  scottl
% updated argument handling, ability to strip unwanted characters, etc.
%
% Revision 1.2  2006/10/18 16:03:15  scottl
% add first character of line counts (for use in HMM's)
%
% Revision 1.1  2006/06/10 21:01:40  scottl
% Initial revision.
%

% LOCAL VARS %
%%%%%%%%%%%%%%
num_files = size(Files,1);
D = {};
W = {};
C = [];

eol = 10;  %encoding for the character corresponding to the newline character

%list of character to disclude from out counts.  Typically these will be the
%non-printable chars though it can be customized
strip_chars = [1:31,127];


% CODE START %
%%%%%%%%%%%%%%
tic;
for ii=1:num_files
    fid = fopen(Files{ii});
    if fid == -1
        error('unable to open file');
    end
    W = [W; textscan(fid, '%s')];
    frewind(fid);
    C = [C; fread(fid)];
    fclose(fid);
end
fprintf('%.2fs: finished parsing files\n', toc);

[D.word,ii,ii] = unique(W{1}(:));
D.word_count = accumarray(ii,1);
fprintf('%.2fs: finished creating and counting words\n', toc);

[D.char,ii,ii] = unique(C);
D.char = char(D.char);
D.char_count = accumarray(ii,1);
fprintf('%.2fs: finished creating and counting chars\n', toc);

num_chars = length(D.char);
D.char_bigram = zeros(num_chars);
D.char_trigram = zeros(num_chars, num_chars, num_chars);
D.first_count = zeros(num_chars, 1);
D.first_bg_count = zeros(num_chars);

pos = zeros(1, max(single(D.char)));
for ii=1:num_chars
    pos(single(D.char(ii))) = ii;
end

D.first_count(pos(C(1))) = 1;  %the 1st character is the start of the 1st line
D.first_bg_count(pos(C(2)),pos(C(2))) = 1;  %the 1st 2 chars on 1st line
D.char_bigram(pos(C(1)),pos(C(2))) = 1;  %count the 1st 2 chars separately

for ii=3:length(C)
    prev_prev = pos(C(ii-2));
    prev = pos(C(ii-1));
    curr = pos(C(ii));
    D.char_trigram(prev_prev,prev,curr) = ...
    D.char_trigram(prev_prev,prev,curr) + 1;
    D.char_bigram(prev,curr) = D.char_bigram(prev,curr) + 1;
    if C(ii-2) == eol
        D.first_bg_count(prev,curr) = D.first_bg_count(prev,curr) + 1;
    elseif C(ii-1) == eol
        D.first_count(curr) = D.first_count(curr) + 1;
    end
end
fprintf('%.2fs: finished creating character n-gram matrices\n', toc);

[D.char, idx] = setdiff(D.char, char(strip_chars)');
D.char_count = D.char_count(idx);
D.first_count = D.first_count(idx);
D.first_bg_count = D.first_bg_count(idx,idx);
D.char_bigram = D.char_bigram(idx,idx);
D.char_trigram = D.char_trigram(idx,idx,idx);
fprintf('%.2fs: finished stripping unwanted characters\n', toc);
