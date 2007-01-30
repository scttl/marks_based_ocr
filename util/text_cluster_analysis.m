function text_cluster_analysis(Syms, files, varargin)
% TEXT_CLUSTER_ANALYSIS    Use ASCII text to cluster and evaluate mappings
%
%     text_cluster_analysis(SYMS, FILE_LIST, [VAR1, VAL1]...)
%
% FILE_LIST should be either a listing of page files representing the documents 
% to process and collect statistics from, or it can be a directory to be 
% searched for documents.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: text_cluster_analysis.m,v 1.2 2007-01-30 01:29:35 scottl Exp $
%
% REVISION HISTORY
% $Log: text_cluster_analysis.m,v $
% Revision 1.2  2007-01-30 01:29:35  scottl
% removed old debugging code
%
% Revision 1.1  2007-01-13 18:36:46  scottl
% initial check-in.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%pattern representing filenames to look for if given a directory
text_pattern = '*.Z*';

%pattern to group multiple files together into the same document
doc_group_pattern = '(\d{4,4})\_.*';
doc_group_replace = '$1';

%set the font size of the axis labels and units
ft_size = 18;


% CODE START %
%%%%%%%%%%%%%%

if nargin < 2
    error('incorrect number of arguments passed');
elseif nargin > 2
    process_optional_args(varargin{:});
end

if ~iscell(files)
    files = {files};
end

%collect the list of report files from the files and directories passed.
text_list = {};
for ii=1:length(files)
    if isdir(files{ii});
        cmd = ['find ', files{ii}, ' -name "', text_pattern, '" -print'];
        [s,w] = unix(cmd);
        if s~=0
            error('problems running cmd: %s', cmd);
        end
        new_list = convert_to_cell(w);
        text_list = [text_list, new_list(:)];
    elseif exist(files{ii}, 'file')
        text_list{end+1} = files{ii}
    end
end
docs = unique(regexprep(text_list, doc_group_pattern, doc_group_replace));

%collect statistics by reading each of the report files
num_docs = length(docs);
fprintf('GENERATING STATISTICS FOR %d DOCUMENTS\n', num_docs);
fprintf('--------------------------------------\n');
tot_a = zeros(num_docs,3);
for ii=1:num_docs
    idx = strmatch(docs{ii}, text_list);
    Clust = create_word_dictionary(text_list(idx));
    Clust.num = length(Clust.char);
    total_chars = sum(Clust.char_count);
    for jj=1:length(Clust.pos_count)
        Clust.pos_count{jj} = Clust.pos_count{jj} ./ total_chars;
    end
    [order, score] = positional_learn_mappings(Clust, Syms);
    tot_a(ii,1) = total_chars;
    clust_chars = cellstr(Clust.char);
    clust_chars{Clust.char == ' '} = ' '; %since cellstr munches space chars
    tot_a(ii,2) = sum(Clust.char_count(~strcmp(clust_chars, ...
                      Syms.val(order(:,1)))));
    tot_a(ii,3) = 1 - (tot_a(ii,2) / tot_a(ii,1));
end

%create a scatter plot of character accuracy versus size
scatter(tot_a(:,1), tot_a(:,3)*100, 100, 'filled');
xlabel('number of characters in document', 'FontSize', ft_size);
ylabel('percentage of characters correctly identified', 'FontSize', ft_size);
set(gca, 'FontSize', ft_size);
fprintf('average number of characters per document: %.2f\n', mean(tot_a(:,1)));
fprintf('average num of incorrectly identified chars: %.2f\n',mean(tot_a(:,2)));
fprintf('average character accuracy per document: %.2f\n', mean(tot_a(:,3)));
fprintf('minimum document length (chars): %d\n', min(tot_a(:,1)));
fprintf('maximum document length (chars): %d\n', max(tot_a(:,1)));
fprintf('median document length (chars): %d\n', median(tot_a(:,1)));
fprintf('minimum character accuracy: %.2f\n', min(tot_a(:,3)));
fprintf('maximum character accuracy: %.2f\n', max(tot_a(:,3)));
fprintf('median character accuracy: %.2f\n', median(tot_a(:,3)));


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this small helper converts a long newline delimited string into a cell array
%with newlines used to break each string into the appropriate entries
function cell_str = convert_to_cell(in_str)
    idx = strfind(in_str, char(10));  %look for newlines
    start_pos = 1;
    cell_str = {};
    while ~isempty(idx)
        cell_str{end+1} = in_str(start_pos:idx(1)-1);
        while length(idx) > 1 && idx(2) == idx(1) + 1
            %remove consecutive newlines
            idx = idx(2:end);
        end
        start_pos = idx(1)+1;
        idx = idx(2:end);
    end
    if start_pos < length(in_str)
        cell_str{end+1} = in_str(start_pos:end);
    end
  

