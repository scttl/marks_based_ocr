function tot_a = unlv_word_ocr_analysis(files, varargin)
% UNLV_WORD_OCR_ANALYSIS    Collect statistics and plot average word accuracy
%
%     TOTALS = unlv_word_ocr_analysis(FILE_LIST, [VAR1, VAL1]...)
%
% FILE_LIST should be either a listing of accuracy reports to process and 
% collect statistics from, or it can be a directory to be searched for accuracy 
% reports to process
%
% TOTALS represents a matrix specifying the accuracy for each report.
% The first column gives the total number of words appearing, the second 
% the number incorrectly identified, and the third the accuracy
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: unlv_word_ocr_analysis.m,v 1.5 2007-02-05 22:16:21 scottl Exp $
%
% REVISION HISTORY
% $Log: unlv_word_ocr_analysis.m,v $
% Revision 1.5  2007-02-05 22:16:21  scottl
% fixed handling of multiple directories.
%
% Revision 1.4  2007-02-01 17:58:36  scottl
% implemented ability to limit documents considered for generating
% statistics over, removed plot of dashed lines between points
%
% Revision 1.3  2007-01-30 01:31:33  scottl
% don't read the entire file since this could be large, instead we now
% just read the first portion when looking for accuracy info
%
% Revision 1.2  2007-01-18 19:16:36  scottl
% updates to display more digits of accuracies and other stats
%
% Revision 1.1  2007-01-13 18:36:47  scottl
% initial check-in.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%pattern representing filenames to look for if given a directory
rprt_pattern = '*.wordtot_rprt';

%set the font size of the axis labels and units
ft_size = 18;

%to compare with other results, we can choose to only take a best subset
take_best_pct = 1.0;

%we can limit to processing only some of the reports found.  Only really useful
%if processing dictionaries and you know in advance how many you want to proces
%Leave empty to include all reports found
take_first_num = [];

%what if anything should we display a plot of
plot_total=true;
total_style = '*';


% CODE START %
%%%%%%%%%%%%%%

if nargin < 1
    error('incorrect number of arguments passed');
elseif nargin > 1
    process_optional_args(varargin{:});
end

if ~iscell(files)
    files = {files};
end

%collect the list of report files from the files and directories passed.
rprt_list = {};
for ii=1:length(files)
    if isdir(files{ii});
        cmd = ['find ', files{ii}, ' -name "', rprt_pattern, '" -print'];
        [s,w] = unix(cmd);
        if s~=0
            error('problems running cmd: %s', cmd);
        end
        new_rprts = convert_to_cell(w);
        rprt_list = [rprt_list(:); new_rprts(:)]';
    elseif exist(files{ii}, 'file')
        rprt_list{end+1} = files{ii}
    end
end

%collect statistics by reading each of the report files
num_rprts = length(rprt_list);
if ~isempty(take_first_num) && take_first_num > 0 && take_first_num <= num_rprts
    num_rprts = take_first_num;
end
num_best_rprts = floor(take_best_pct * num_rprts);
fprintf('GENERATING STATISTICS FOR %d WORD REPORTS\n', num_best_rprts);
fprintf('------------------------------------------\n');
tot_a = zeros(num_rprts,3);
for ii=1:num_rprts
    %in the UNLV reports, we only really care about the ASCII count accuracy
    %so we extract them from the file
    fid = fopen(rprt_list{ii}, 'r');
    if fid == -1
        error('problems reading file: %s', rprt_list{ii});
    end
    %only need to process first 175 bytes of file since the accuracy is written
    %near the start of the file
    data = char(fread(fid, 175))';
    fclose(fid);
    data = convert_to_cell(data);

    if length(data) >= 5
        tot_a(ii,1) = strread(data{3}, '%d %*[^\n]');
        tot_a(ii,2) = strread(data{4}, '%d %*[^\n]');
        tot_a(ii,3) = 1 - (tot_a(ii,2) / tot_a(ii,1));
    else
        error('incorrectly formatted data file');
    end
end

%sort them in descending order, and chop those if not generating stats over the
%full subset
if take_best_pct < 1
    [idx,idx] = sort(tot_a(:,3), 'descend');
    tot_a = tot_a(idx,:);
    tot_a = tot_a(1:num_best_rprts,:);
end

if plot_total
    %create a plot of accuracy accuracy versus document size
    plot(tot_a(:,1), tot_a(:,3)*100, total_style);
    xlabel('number of words in document', 'FontSize', ft_size);
    ylabel('% of words correctly identified', 'FontSize', ft_size);
    set(gca, 'FontSize', ft_size);
end

fprintf('average number of words per document: %.4f\n', mean(tot_a(:,1)));
fprintf('average num of incorrectly identified words: %.4f\n',mean(tot_a(:,2)));
fprintf('average word accuracy per document: %.4f\n', mean(tot_a(:,3)));
fprintf('minimum document length (words): %d\n', min(tot_a(:,1)));
fprintf('maximum document length (words): %d\n', max(tot_a(:,1)));
fprintf('median document length (words): %d\n', median(tot_a(:,1)));
fprintf('minimum word accuracy: %.4f\n', min(tot_a(:,3)));
fprintf('maximum word accuracy: %.4f\n', max(tot_a(:,3)));
fprintf('median word accuracy: %.4f\n', median(tot_a(:,3)));


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
  

