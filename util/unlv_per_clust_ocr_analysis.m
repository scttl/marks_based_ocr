function res = unlv_per_clust_ocr_analysis(files, varargin)
% UNLV_CLUST_RPRT_OCR_ANALYSIS    Plot per cluster OCR accuracy results
%
%     res = unlv_per_clust_ocr_analysis(FILE_LIST, [VAR1, VAL1]...)
%
% FILE_LIST should be either a listing of accuracy reports to process and 
% collect statistics from, or it can be a directory to be searched for accuracy 
% reports to process
%
% RES will be a 3 column matrix with one row for each cluster (ordered by 
% frequency).  The first column will indicate the total number of times it was 
% seen (over all documents), the second column will indicate the number of times
% it was identified correctly (over all documents), and the third column will 
% should the average recognition accuracy of that cluster (when averaged over 
% each document's accuracy)
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: unlv_per_clust_ocr_analysis.m,v 1.1 2007-05-14 23:21:19 scottl Exp $
%
% REVISION HISTORY
% $Log: unlv_per_clust_ocr_analysis.m,v $
% Revision 1.1  2007-05-14 23:21:19  scottl
% initial check-in
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%pattern representing filenames to look for if given a directory
rprt_pattern = '*.clusttot_rprt';

%set the font size of the axis labels and units
ft_size = 18;

%set the size of the plotted points
marker_size = 12;

%we can limit to processing only some of the reports found.  Only really useful
%if processing dictionaries and you know in advance how many you want to proces
%Leave empty to include all reports found
take_first_num = [];


%what if anything should we display a plot of
plot_per_clust_acc=true;
plot_per_freq_acc=false;
plot_num_clust_hist=false;


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

fprintf('GENERATING STATISTICS FOR %d CLUSTER REPORTS\n', num_rprts);
fprintf('--------------------------------------------\n');

stats = cell(0);
num_clusts = zeros(num_rprts,1);
for ii=1:num_rprts
    fid = fopen(rprt_list{ii}, 'r');
    if fid == -1
        error('problems reading file: %s', rprt_list{ii});
    end
    data = textscan(fid, '%f%f%f');
    fclose(fid);
    num_clusts(ii) = length(data{1});
    for jj = data{1}'
        seen = data{2}(jj);
        corr = data{3}(jj);
        if seen > 0
            acc = corr/seen;
        else
            acc = 0;
        end
        if jj>length(stats)
            stats{jj} = [seen, corr, acc];
        else
            stats{jj} = [stats{jj}; [seen, corr, acc]];
        end
    end
end

max_clust = length(stats);
res = zeros(max_clust,3);
for ii=1:max_clust
    res(ii,1) = sum(stats{ii}(:,1));
    res(ii,2) = sum(stats{ii}(:,2));
    res(ii,3) = mean(stats{ii}(:,3));
end

%create the appropriate plots
fig_num=1;
if plot_num_clust_hist
    figure(fig_num); clf;
    hist(num_clusts, min(num_clusts):max(num_clusts));
    axis tight;
    xlabel('number of clusters', 'FontSize', ft_size);
    ylabel('frequency', 'FontSize', ft_size);
    set(gca, 'FontSize', ft_size);
    fig_num = fig_num+1;
end

if plot_per_clust_acc
    figure(fig_num); clf;
    plot(res(:,3));
    axis tight;
    xlabel('cluster identifier', 'FontSize', ft_size);
    ylabel('average recognition accuracy', 'FontSize', ft_size);
    set(gca, 'FontSize', ft_size);
    fig_num = fig_num+1;
end

%print useful statistics
fprintf('average number of clusters per document: %.4f\n', mean(num_clusts));
fprintf('min number of clusters per document: %.4f\n', min(num_clusts));
fprintf('median number of clusters per document: %.4f\n', median(num_clusts));
fprintf('max number of clusters per document: %.4f\n', max(num_clusts));

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
  

