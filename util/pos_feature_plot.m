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
% $Id: pos_feature_plot.m,v 1.4 2007-02-08 19:22:20 scottl Exp $
%
% REVISION HISTORY
% $Log: pos_feature_plot.m,v $
% Revision 1.4  2007-02-08 19:22:20  scottl
% added bar graph, relabelled axes according to Sam's suggestions.
%
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

%should we draw the plot half-height?
half_height = true;

%should we draw a bar plot instead? (only works when passing a single vector
draw_bar = true;

%how should we set markers, lines, and color?
plot_marker_string = 'o';  %type help plot for choices
plot_marker_color = 'b';
plot_marker_size = 8;
plot_xlabel = 'Word Length';
plot_ylabel = 'Positional Probability';
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

if half_height
    figure(1); clf; set(gcf, 'PaperPosition', [0 0 7 2.5]);
    hold on;
end

sz = size(counts);
if sz(1) == 1 && draw_bar
    bar(counts, 'FaceColor', plot_marker_color)
else
    plot(1:sz(2), counts, plot_marker_string, ...
         'MarkerFaceColor', plot_marker_color, ...
         'MarkerEdgeColor', plot_marker_color, ...
         'MarkerSize', plot_marker_size);
end
xlabel(plot_xlabel, 'FontSize', ft_size);
ylabel(plot_ylabel, 'FontSize', ft_size);
title(plot_title);
if ~isempty(plot_legend_strings)
    legend(plot_legend_strings);
end
axis([-2 sz(2)+1 0 1]);
set(gca, 'YTick', [0 .2 .4 .6 .8 1.0]);

if draw_word_dividers
    %positive root of quadratic equation used to caculate number of words
    num_words = floor((sqrt(1 + 8*sz(2)) -1)/2);
    dividers=[0,cumsum(1:num_words)+.5];
    centres=dividers(1:num_words)+diff(dividers)/2;
    line([1;1]*dividers,[0;1]*ones(1,16), 'LineWidth',1, 'LineStyle', ...
         word_line_style, 'Color', word_line_color);
    set(gca, 'XTick', centres);
    set(gca, 'XTickLabel', num2str((1:num_words)'));
end

