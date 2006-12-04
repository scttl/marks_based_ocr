function plot_clust_sizes(Clust, varargin)
%  PLOT_CLUST_SIZES  Plot the number of components belonging to each cluster
%
%   PLOT_CLUST_SIZES(CLUST, [VAR1, VAL1]...)
%
%   CLUST should be a struct like that returned in cluster_comps()
%
%   optional LOCAL VARS values below can be overriden specifying the name and
%   new value for the variable to be overwritten as additinoal parameters.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: plot_clust_sizes.m,v 1.1 2006-12-04 19:21:08 scottl Exp $
%
% REVISION HISTORY
% $Log: plot_clust_sizes.m,v $
% Revision 1.1  2006-12-04 19:21:08  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%strings controlling plot labels
title_str = 'Plot showing # of components per cluster';
xaxis_str = 'cluster';
yaxis_str = 'number of components';

%plot symbols see plot()
sym_str = 'b.-';

%this can be used to display counts for a subset of the clusters (leave
%empty to display all)
plot_range = [];


%set save_plot to true to write the plot images to disk based on the params
%below it
save_plot = false;
global MOCR_PATH;  %make use of the globally defined MOCR_PATH variable
img_prefix = [MOCR_PATH, '/results/clust_size_plot'];
img_format = '-dpng'; %other choices: -deps, -depsc2, etc. see print()


% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 1
    error('incorrect number of arguments specified!');
elseif nargin > 1
    process_optional_args(varargin{:});
end

if isempty(plot_range)
    plot(Clust.num_comps, sym_str);
else
    plot(plot_range, Clust.num_comps(plot_range), sym_str);
end
title(title_str);
xlabel(xaxis_str);
ylabel(yaxis_str);

%save the plot to disk if required.
if save_plot
    fprintf('%.2fs: writing plot image to disk\n', toc);
    print(gcf, img_format, img_prefix);
end

fprintf('Elapsed time: %f\n', toc);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
