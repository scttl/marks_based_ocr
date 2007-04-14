function [words, per_clust_words] = get_cluster_word_seq(Comps, delim, varargin)
% GET_CLUSTER_WORD_SEQ  Determine reading order sequence of cluster words
%
%   [WORDS,PER_CLUST_WORDS] = get_cluster_word_seq(COMPS, DELIM, [var1,val1]...)
%   This function uses component and delimiter information to tokenize the
%   reading order sequence of cluster identifiers into ``words''.  Typically
%   the delimiter will be the found space cluster.
%
%   COMPS should be a struct like that returned from get_comps()
%
%   DELIM should be a vector of cluster identifiers used to tokenize the
%   sequence.  There can be more than one delimiter, and each of them will not
%   appear in the final list of words (by default).  Tokens are also delimited
%   when the end of a line is reached.
%
%   WORDS will be a cell array, each entry of which will be a delimited word
%         (words are listed in reading order)
%
%   PER_CLUST_WORDS is also a cell array, with one entry per cluster.  Each
%   will be a vector listing all the words to which this cluster belongs.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: get_cluster_word_seq.m,v 1.2 2007-04-14 20:16:00 scottl Exp $
%
% REVISION HISTORY
% $Log: get_cluster_word_seq.m,v $
% Revision 1.2  2007-04-14 20:16:00  scottl
% smal change to ensure there is one word per row.
%
% Revision 1.1  2007-04-10 19:39:43  scottl
% initial check-in.  Separated from word_lookup_map
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%should we include delimiter characters in the tokenized strings?
include_delims = false;

%should tokens be brokend when the end of a line is reached?
eol_delim = true;


% CODE START %
%%%%%%%%%%%%%%
if nargin < 2
    error('incorrect number of arguments');
elseif nargin > 2
    process_optional_args(varargin{:});
end

seq = get_cluster_seq(Comps, unique(Comps.line));
words = cell(0);
per_clust_words = cell(max(Comps.clust),1);
word_num = 1;
for ii=1:length(seq)
    curr_word = [];
    curr_line = seq{ii};
    while ~isempty(curr_line)
        if any(curr_line(1) == delim)
            %end of a word.  Update counts if not empty (ie ignore repeated
            %delimiters)
            if ~isempty(curr_word)
                if include_delims
                    curr_word = [curr_word, curr_line(1)];
                end
                words{word_num,1} = curr_word;
                curr_word = [];
                word_num = word_num + 1;
            end
        else
            %continue to grow the current word
            curr_word = [curr_word, curr_line(1)];
            if all(per_clust_words{curr_line(1)} ~= word_num)
                %ensure each cluster is only added once for each word
                per_clust_words{curr_line(1)} = ...
                                  [per_clust_words{curr_line(1)},word_num];
            end
        end
        curr_line = curr_line(2:end);
    end
    if eol_delim && ~isempty(curr_word)
        if include_delims
            curr_word = [curr_word, curr_line(1)];
        end
        words{word_num,1} = curr_word;
        curr_word = [];
        word_num = word_num + 1;
    end
end
