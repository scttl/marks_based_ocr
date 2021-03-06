function D = create_word_dictionary(Files, varargin)
% CREATE_WORD_DICTIONARY Parse text corpus creatng list and counts of words and 
% chars
%
%   D = CREATE_WORD_DICTIONARY(Files, [var1, val1]...)
%   Files should be a cell array (or single string), and give the full 
%   path and name of the file(s) to be used in creating the dictionary.
%
%   D is a struct with fields: word, word_count, char, char_count, pos_count,
%   char_bigram, char_trigram, first_count, and first_bg_count
%
%  NOTE: This requires matlab7 since it uses textscan
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: create_word_dictionary.m,v 1.8 2007-01-25 18:50:29 scottl Exp $
%
% REVISION HISTORY
% $Log: create_word_dictionary.m,v $
% Revision 1.8  2007-01-25 18:50:29  scottl
% added ability to generate statistics for some prefix of the data instead of
% the entire set.
%
% Revision 1.7  2007-01-13 18:20:26  scottl
% ensure that words and word counts are handled correctly when
% multiple files are used.  Speed up positional count creation.
%
% Revision 1.6  2007-01-05 17:17:58  scottl
% added additional skip characters (limit to basic ASCII symbols for now).
%
% Revision 1.5  2006-12-19 21:38:19  scottl
% added character positional count information.
%
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

%list of character to disclude from our counts.  Typically these will be the
%non-printable chars though it can be customized
strip_chars = [1:31,127:999];

%when creating the positional counts, up to what word length should we include?
max_word_len = 10;

%should we limit ourselves to just some first portion of the dataset?  Set this
%to inf, to keep the entire dataset.  Otherwise, it will specify the number of 
%(8-bit) characters to keep
keep_num_chars = inf;


% CODE START %
%%%%%%%%%%%%%%
tic;
if nargin < 1
    error('incorrect number of args specified');
elseif nargin > 1
    process_optional_args(varargin{:});
end

if ~iscell(Files)
    Files = {Files};
end
for ii=1:num_files
    if keep_num_chars <= 0
        break;
    end
    fid = fopen(Files{ii});
    if fid == -1
        error('unable to open file');
    end
    [tmp,count] = fread(fid, keep_num_chars);
    fclose(fid);
    keep_num_chars = keep_num_chars - count;
    C = [C; tmp];
    tmp = textscan(char(tmp), '%s');
    W = [W; tmp{:}];
end

fprintf('%.2fs: finished parsing files\n', toc);
[D.word,ii,ii] = unique(W);
D.word_count = accumarray(ii,1);
D.word_len = zeros(size(D.word_count));
for ii=1:length(D.word_len)
    D.word_len(ii) = length(D.word{ii});
end
fprintf('%.2fs: finished creating and counting words\n', toc);

[D.char,ii,ii] = unique(C);
D.char = char(D.char);
D.char_count = accumarray(ii,1);
D.pos_count = cell(1,max_word_len);
for ii=1:max_word_len
    this_count = zeros(length(D.char),ii);
    idx = find(D.word_len == ii);
    if ~isempty(idx)
        words = cell2mat(D.word(idx));
        for jj=1:length(D.char)
            for kk=1:ii
                this_count(jj,kk) = sum(D.word_count(idx(strcmp(words(:,kk), ...
                                        {D.char(jj)}))));
            end
        end
    end
    D.pos_count{ii} = this_count;
end
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
for ii=1:max_word_len
    D.pos_count{ii} = D.pos_count{ii}(idx,:);
end
D.first_count = D.first_count(idx);
D.first_bg_count = D.first_bg_count(idx,idx);
D.char_bigram = D.char_bigram(idx,idx);
D.char_trigram = D.char_trigram(idx,idx,idx);
fprintf('%.2fs: finished stripping unwanted characters\n', toc);
