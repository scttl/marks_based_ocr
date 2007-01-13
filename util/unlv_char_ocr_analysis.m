function unlv_char_ocr_analysis(files, varargin)
% UNLV_CHAR_OCR_ANALYSIS    Collect statistics and plot average char accuracy
%
%     unlv_char_ocr_analysis(FILE_LIST, [VAR1, VAL1]...)
%
% FILE_LIST should be either a listing of accuracy reports to process and 
% collect statistics from, or it can be a directory to be searched for accuracy 
% reports to process
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: unlv_char_ocr_analysis.m,v 1.1 2007-01-13 18:36:47 scottl Exp $
%
% REVISION HISTORY
% $Log: unlv_char_ocr_analysis.m,v $
% Revision 1.1  2007-01-13 18:36:47  scottl
% initial check-in.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%pattern representing filenames to look for if given a directory
rprt_pattern = '*.chartot_rprt';

%set the font size of the axis labels and units
ft_size = 18;

%to compare with other results, we can choose to only take a best subset
take_best_pct = 1.0;


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
        rprt_list = [rprt_list, new_rprts(:)];
    elseif exist(files{ii}, 'file')
        rprt_list{end+1} = files{ii}
    end
end

%collect statistics by reading each of the report files
num_rprts = length(rprt_list);
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
%full subset
if take_best_pct < 1
    [idx,idx] = sort(tot_a(:,3), 'descend');
    tot_a = tot_a(idx,:);
    tot_a = tot_a(1:num_best_rprts,:);
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
  

