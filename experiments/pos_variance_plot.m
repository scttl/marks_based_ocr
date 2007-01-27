function pos_variance_plot(symbol, varargin)
% POS_VARIANCE_PLOT  Collect positional counts and plot them versus doc length
%    pos_variance_plot(symbol, [VAR1, VAL1]...)
%
%  symbol should be a string representing the symbol to collect positional
%  counts and plot the results for
%
%  Files is an optional cell array listing the input files to collect positional
%  counts of the symbol from, and average over.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: pos_variance_plot.m,v 1.1 2007-01-27 01:06:48 scottl Exp $
%
% REVISION HISTORY
% $Log: pos_variance_plot.m,v $
% Revision 1.1  2007-01-27 01:06:48  scottl
% initial check-in.
%

global MOCR_PATH;  %used to determine where to save results

% LOCAL VARS %
%%%%%%%%%%%%%%

%these point at the files containing the actual text of the corpus
Files = {[MOCR_PATH, '/data/reuters/000.body.txt'];
        [MOCR_PATH, '/data/reuters/001.body.txt'];
        [MOCR_PATH, '/data/reuters/002.body.txt'];
        [MOCR_PATH, '/data/reuters/003.body.txt'];
        [MOCR_PATH, '/data/reuters/004.body.txt'];
        [MOCR_PATH, '/data/reuters/005.body.txt'];
        [MOCR_PATH, '/data/reuters/006.body.txt'];
        [MOCR_PATH, '/data/reuters/007.body.txt'];
        [MOCR_PATH, '/data/reuters/008.body.txt'];
        [MOCR_PATH, '/data/reuters/009.body.txt'];
        [MOCR_PATH, '/data/reuters/010.body.txt'];
        [MOCR_PATH, '/data/reuters/011.body.txt'];
        [MOCR_PATH, '/data/reuters/012.body.txt'];
        [MOCR_PATH, '/data/reuters/013.body.txt'];
        [MOCR_PATH, '/data/reuters/014.body.txt'];
        [MOCR_PATH, '/data/reuters/015.body.txt'];
        [MOCR_PATH, '/data/reuters/016.body.txt'];
        [MOCR_PATH, '/data/reuters/017.body.txt'];
        [MOCR_PATH, '/data/reuters/018.body.txt'];
        [MOCR_PATH, '/data/reuters/019.body.txt'];
        [MOCR_PATH, '/data/reuters/020.body.txt'];
        [MOCR_PATH, '/data/reuters/021.body.txt']};

%up to what length words should we include stats and plots for
max_word_len = 10;

%at what document lengths (number of characters) should we take positional
%counts
doc_len_intervals = 100:1000:20000;

save_figure = false;
figure_driver = '-depsc2'; %help print for other choices
figure_res = '-r300';  %output resolution DPI
figure_file_prefix = [MOCR_PATH, '/results/pos_variance_'];
figure_file_suffix = '.eps';


% CODE START %
%%%%%%%%%%%%%%
if nargin < 1
    error('must specify the symbol to collect stats for');
elseif nargin > 1
    process_optional_args(varargin{:});
end

num_pos = (max_word_len * (max_word_len+1)) / 2;
final_stats = zeros(length(doc_len_intervals), num_pos);
for ii=1:length(doc_len_intervals)
    stats = [];
    for jj=1:length(Files)
        D = create_word_dictionary(Files{jj}, 'keep_num_chars', ...
            doc_len_intervals(ii));
        idx = strfind(D.char', symbol);
        if isempty(idx)
            stats = [stats; zeros(1,num_pos)];
        else
            xx = cell2mat(D.pos_count);
            xx = xx ./ sum(D.char_count);
            stats = [stats; xx(idx,:)];
        end
    end
    final_stats(ii,:) = var(stats);
    fprintf('finished interval %d\n', doc_len_intervals(ii));
end

%create plots showing the resultant variance for each word length
num_rows = floor(sqrt(max_word_len));
num_cols = ceil(max_word_len / num_rows);
start_col = 1;
for ii=1:max_word_len
    subplot(num_rows,num_cols,ii);
    plot(doc_len_intervals, final_stats(:,start_col:start_col+ii-1));
    xlabel('num chars');
    ylabel('var');
    title(['''',symbol,''' pv, wl=', num2str(ii)]);
    %legend show;
    start_col = start_col + ii;
end

if save_figure
    print(gcf, figure_driver, figure_res, ...
          [figure_file_prefix, symbol, figure_file_suffix]);
end
