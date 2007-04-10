function plot_space_model(l1, l2, c, S, varargin)
%  PLOT_SPACE_MODEL  Plot space counts and estimated 2-poisson mix model
%
%   PLOT_SPACE_MODEL(L1, L2, C, S, [VAR1, VAL1]...)
%
%   L1, L2, C are the 2-Poisson mixture model paramters (lambda1, lambda2, 
%   and c), as described in the Huang2006 paper.
%
%   S should be a vector of space widths
%
%   optional LOCAL VARS values below can be overriden specifying the name and
%   new value for the variable to be overwritten as additinoal parameters.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: plot_space_model.m,v 1.1 2007-04-10 18:53:33 scottl Exp $
%
% REVISION HISTORY
% $Log: plot_space_model.m,v $
% Revision 1.1  2007-04-10 18:53:33  scottl
% initial check-in
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%strings controlling plot labels
title_str = 'Histogram plot of space widths and the estimated model';
xaxis_str = 'space width';
yaxis_str = 'count';

%plot symbols see plot()
p1_str = 'g.--';  %first poisson dist'n
p2_str = 'r.--';  %2nd poisson dist'n

%up to how wide should spaces be?
max_space_width = 75;

%set save_plot to true to write the plot images to disk based on the params
%below it
save_plot = false;
global MOCR_PATH;  %make use of the globally defined MOCR_PATH variable
img_prefix = [MOCR_PATH, '/results/space_width_plot'];
img_format = '-dpng'; %other choices: -deps, -depsc2, etc. see print()


% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 4
    error('incorrect number of arguments specified!');
elseif nargin > 4
    process_optional_args(varargin{:});
end

%create a histogram of space widths
clf;
num = length(S);
S = round(S);
S(S < 1) = 1;
S(S > max_space_width) = max_space_width;
hist(S, 1:max_space_width);
[maxfreq, maxfreq] = mode(S);
hold on;

%plot the Poisson estimates given the current lambda values over the entire
%range of spaces
x = 1:max_space_width;
facx = factorial(x);
num_p1 = sum(S < c);
num_p2 = num - num_p1;
if num_p2 > num_p1 || l1 > l2
    %swap amount assigned to each peak (by default we assume the first
    %interchar spacing peak should contain more mass)
    tmp = num_p1;
    num_p1 = num_p2;
    num_p2 = tmp;
end
p1 = num_p1 * ((exp(-l1) * l1.^x) ./ facx);
p2 = num_p2 * ((exp(-l2) * l2.^x) ./ facx);
plot(x,p1,p1_str);
plot(x,p2,p2_str);

%add a vertical line at c
line([c;c], [maxfreq;0]);

title(title_str);
xlabel(xaxis_str);
ylabel(yaxis_str);
hold off;

%save the plot to disk if required.
if save_plot
    fprintf('%.2fs: writing plot image to disk\n', toc);
    print(gcf, img_format, img_prefix);
end

fprintf('Elapsed time: %f\n', toc);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
