function pos_feature_plot(counts, varargin)
% POS_FEATURE_PLOT    Plot positional features given
%
%     pos_feature_plot(COUNTS, [VAR1, VAL1]...)
%
% COUNTS should either be a matrix or a cell array.  With each row representing
% a particular character to plot, and the columns representing the positional
% counts to be plotted.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: pos_feature_plot.m,v 1.3 2007-02-03 21:06:11 scottl Exp $
%
% REVISION HISTORY
% $Log: pos_feature_plot.m,v $
% Revision 1.3  2007-02-03 21:06:11  scottl
% default font changes
%
% Revision 1.2  2007-01-30 01:29:13  scottl
% added optional legend display, reset axes to be tight to the data
%
% Revision 1.1  2007-01-26 16:57:45  scottl
% initial check-in.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%how should we set markers, lines, and color?
plot_marker_string = 'o';  %type help plot for choices
plot_marker_color = 'b';
plot_marker_size = 8;
plot_xlabel = 'Position and word length';
plot_ylabel = 'Probability';
%plot_title = 'Plot of probability vs position in words of various length';
plot_title = '';
%no legend is used by default.  If overriding this should be a cell array of 
%strings to name each vector
plot_legend_strings = '';  

%set the font size of the axis labels and units
ft_size = 18;

%should we draw word length delimeter lines?
draw_word_dividers = true;
word_line_style = '--';  %type doc line to see the property vals
word_line_color = 'k';


% CODE START %
%%%%%%%%%%%%%%

if nargin < 1
    error('incorrect number of arguments passed');
elseif nargin > 1
    process_optional_args(varargin{:});
end

if iscell(counts)
    counts = cell2mat(counts);
end

sz = size(counts);
plot(1:sz(2), counts, plot_marker_string, ...
    'MarkerFaceColor', plot_marker_color, ...
    'MarkerEdgeColor', plot_marker_color, ...
    'MarkerSize', plot_marker_size);
xlabel(plot_xlabel, 'FontSize', ft_size);
ylabel(plot_ylabel, 'FontSize', ft_size);
title(plot_title);
if ~isempty(plot_legend_strings)
    legend(plot_legend_strings);
end
axis tight;

if draw_word_dividers
    %positive root of quadratic equation used to caculate number of words
    num_words = floor((sqrt(1 + 8*sz(2)) -1)/2);
    X = repmat(cumsum(1:num_words)+.5,2,1);
    Y = repmat([1;0], 1, num_words);
    line(X,Y, 'LineStyle', word_line_style, 'Color', word_line_color);
end

