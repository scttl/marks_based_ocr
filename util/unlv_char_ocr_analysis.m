function [tot_a low_a upp_a dig_a sym_a spc_a] = ...
         unlv_char_ocr_analysis(files, varargin)
% UNLV_CHAR_OCR_ANALYSIS    Collect statistics and plot average char accuracy
%
%     [TOTALS LOWER UPPER DIGITS SYMBOLS SPACES] = ...
%     unlv_char_ocr_analysis(FILE_LIST, [VAR1, VAL1]...)
%
% FILE_LIST should be either a listing of accuracy reports to process and 
% collect statistics from, or it can be a directory to be searched for accuracy 
% reports to process
%
% TOTALS, LOWER, UPPER, DIGITS, SYMBOLS, SPACES represent matrcies specifying 
% the accuracy for each report broken up by symbol type.  The first column gives
% the total number of characters appearing, the second the number incorrectly 
% identified, and the third the accuracy
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: unlv_char_ocr_analysis.m,v 1.5 2007-02-05 22:15:34 scottl Exp $
%
% REVISION HISTORY
% $Log: unlv_char_ocr_analysis.m,v $
% Revision 1.5  2007-02-05 22:15:34  scottl
% fixed multiple directory handling, added display of other averages.
%
% Revision 1.4  2007-02-01 17:57:57  scottl
% added ability to limit number of docs reports are calculated over,
% removed dashed lines between points when plotting (by default)
%
% Revision 1.3  2007-01-30 01:30:31  scottl
% added display of lower, upper, digit, and other symbol accuracy to the
% overall accuracy plot.
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
rprt_pattern = '*.chartot_rprt';

%set the font size of the axis labels and units
ft_size = 18;

%set the size of the plotted points
marker_size = 12;

%to compare with other results, we can choose to only take a best subset
take_best_pct = 1.0;

%we can limit to processing only some of the reports found.  Only really useful
%if processing dictionaries and you know in advance how many you want to proces
%Leave empty to include all reports found
take_first_num = [];


%what if anything should we display a plot of
plot_total=true;
total_style = '*';
total_legend = 'Combined';
plot_lowlet=true;
lowlet_style = 'x';
lowlet_legend = 'Low';
plot_upplet=true;
upplet_style = '+';
upplet_legend = 'Upper';
plot_digits=true;
digits_style = 'o';
digits_legend = 'Digits';
plot_other=true;
other_style = 's';
other_legend = 'Other';
plot_spaces=false;
spaces_style = 'd';
spaces_legend = 'Spaces';

%by default we plot accuracy relative to the total number of occurences of 
%symbols in the document.  Set to false to plot relative to the number of occurences of that type of symbol.
renorm_by_total = true;


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
fprintf('GENERATING STATISTICS FOR %d CHAR REPORTS\n', num_best_rprts);
fprintf('------------------------------------------\n');
spc_a = zeros(num_rprts,3);
sym_a = zeros(num_rprts,3);
dig_a = zeros(num_rprts,3);
upp_a = zeros(num_rprts,3);
low_a = zeros(num_rprts,3);
tot_a = zeros(num_rprts,3);
for ii=1:num_rprts
    %in the UNLV reports, we only really care about the ASCII count accuracy
    %so we extract them from the file
    fid = fopen(rprt_list{ii}, 'r');
    if fid == -1
        error('problems reading file: %s', rprt_list{ii});
    end
    data = char(fread(fid))';
    fclose(fid);
    data = convert_to_cell(data);
    while ~isempty(data)
        if strfind(data{1}, 'ASCII')
            [x(1), x(2), x(3)] = strread(data{1}, '%d %d %f %*[^\n]');
            if strfind(data{1}, 'Spacing Characters')
                spc_a(ii,:) = x;
            elseif strfind(data{1}, 'Special Symbols')
                sym_a(ii,:) = x;
            elseif strfind(data{1}, 'Digits')
                dig_a(ii,:) = x;
            elseif strfind(data{1}, 'Uppercase Letters')
                upp_a(ii,:) = x;
            elseif strfind(data{1}, 'Lowercase Letters')
                low_a(ii,:) = x;
            end
        end
        data = data(2:end);
    end
    tot_a(ii,1) = sum(spc_a(ii,1) + sym_a(ii,1) + dig_a(ii,1) + upp_a(ii,1) ...
                      + low_a(ii,1));
    tot_a(ii,2) = sum(spc_a(ii,2) + sym_a(ii,2) + dig_a(ii,2) + upp_a(ii,2) ...
                      + low_a(ii,2));
    tot_a(ii,3) = 1 - (tot_a(ii,2) / tot_a(ii,1));
end

%sort them in descending order, and chop those if not generating stats over the
%full set
if take_best_pct < 1
    [idx,idx] = sort(tot_a(:,3), 'descend');
    tot_a = tot_a(idx,:);
    tot_a = tot_a(1:num_best_rprts,:);
    [idx,idx] = sort(low_a(:,3), 'descend');
    low_a = low_a(idx,:);
    low_a = low_a(1:num_best_rprts,:);
    [idx,idx] = sort(upp_a(:,3), 'descend');
    upp_a = upp_a(idx,:);
    upp_a = upp_a(1:num_best_rprts,:);
    [idx,idx] = sort(sym_a(:,3), 'descend');
    sym_a = sym_a(idx,:);
    sym_a = sym_a(1:num_best_rprts,:);
    [idx,idx] = sort(dig_a(:,3), 'descend');
    dig_a = dig_a(idx,:);
    dig_a = dig_a(1:num_best_rprts,:);
    [idx,idx] = sort(spc_a(:,3), 'descend');
    spc_a = spc_a(idx,:);
    spc_a = spc_a(1:num_best_rprts,:);
end

%renormalizing kludge
if renorm_by_total
    low_a(:,1) = tot_a(:,1);
    upp_a(:,1) = tot_a(:,1);
    dig_a(:,1) = tot_a(:,1);
    sym_a(:,1) = tot_a(:,1);
end

%create a plot of character accuracy versus size
plot_str = 'plot(';
legend_str = 'legend(';
if plot_total
    plot_str = [plot_str, 'tot_a(:,1), tot_a(:,3)*100, total_style, '];
    legend_str = [legend_str, 'total_legend, '];
end
if plot_lowlet
    plot_str = [plot_str, 'low_a(:,1), low_a(:,3), lowlet_style, '];
    legend_str = [legend_str, 'lowlet_legend, '];
end
if plot_upplet
    plot_str = [plot_str, 'upp_a(:,1), upp_a(:,3), upplet_style, '];
    legend_str = [legend_str, 'upplet_legend, '];
end
if plot_digits
    plot_str = [plot_str, 'dig_a(:,1), dig_a(:,3), digits_style, '];
    legend_str = [legend_str, 'digits_legend, '];
end
if plot_other
    plot_str = [plot_str, 'sym_a(:,1), sym_a(:,3), other_style, '];
    legend_str = [legend_str, 'other_legend, '];
end
if plot_spaces
    plot_str = [plot_str, 'spc_a(:,1), spc_a(:,3), spaces_style, '];
    legend_str = [legend_str, 'spaces_legend, '];
end

if plot_str(end) == ' '
    %plotting at least one type of data
    plot_str = [plot_str, ' ''MarkerSize'', ', num2str(marker_size), ');'];
    legend_str = [legend_str, '''Location'', ''Best'');'];
    eval(plot_str);
    eval(legend_str);
    xlabel('number of characters in document', 'FontSize', ft_size);
    ylabel('% of characters correctly identified', 'FontSize', ft_size);
    set(gca, 'FontSize', ft_size);
end

fprintf('average number of characters per document: %.4f\n', mean(tot_a(:,1)));
fprintf('average num of incorrectly identified chars: %.4f\n',mean(tot_a(:,2)));
fprintf('average character accuracy per document: %.4f\n', mean(tot_a(:,3)));
fprintf('minimum document length (chars): %d\n', min(tot_a(:,1)));
fprintf('median document length (chars): %d\n', median(tot_a(:,1)));
fprintf('maximum document length (chars): %d\n', max(tot_a(:,1)));
fprintf('minimum character accuracy: %.4f\n', min(tot_a(:,3)));
fprintf('median character accuracy: %.4f\n', median(tot_a(:,3)));
fprintf('maximum character accuracy: %.4f\n', max(tot_a(:,3)));

%normalize the non-total accuracies
spc_a(:,3) = spc_a(:,3) ./ 100;
low_a(:,3) = low_a(:,3) ./ 100;
upp_a(:,3) = upp_a(:,3) ./ 100;
dig_a(:,3) = dig_a(:,3) ./ 100;
sym_a(:,3) = sym_a(:,3) ./ 100;

fprintf('mean lowercase letter accuracy: %.4f\n', mean(low_a(:,3)));
fprintf('minimum lowercase letter accuracy: %.4f\n', min(low_a(:,3)));
fprintf('median lowercase letter accuracy: %.4f\n', median(low_a(:,3)));
fprintf('maximum lowercase letter accuracy: %.4f\n', max(low_a(:,3)));
fprintf('mean upper letter accuracy: %.4f\n', mean(upp_a(:,3)));
fprintf('minimum upper letter accuracy: %.4f\n', min(upp_a(:,3)));
fprintf('median upper letter accuracy: %.4f\n', median(upp_a(:,3)));
fprintf('maximum upper letter accuracy: %.4f\n', max(upp_a(:,3)));
fprintf('mean digit accuracy: %.4f\n', mean(dig_a(:,3)));
fprintf('minimum digit accuracy: %.4f\n', min(dig_a(:,3)));
fprintf('median digit accuracy: %.4f\n', median(dig_a(:,3)));
fprintf('maximum digit accuracy: %.4f\n', max(dig_a(:,3)));
fprintf('mean other symbol accuracy: %.4f\n', mean(sym_a(:,3)));
fprintf('minimum other symbol accuracy: %.4f\n', min(sym_a(:,3)));
fprintf('median other symbol accuracy: %.4f\n', median(sym_a(:,3)));
fprintf('maximum other symbol accuracy: %.4f\n', max(sym_a(:,3)));
fprintf('mean space symbol accuracy: %.4f\n', mean(spc_a(:,3)));
fprintf('median space symbol accuracy: %.4f\n', median(spc_a(:,3)));
fprintf('maximum space symbol accuracy: %.4f\n', max(spc_a(:,3)));


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
  

