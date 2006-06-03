%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% FILE: startup.m
%%
%% DESCRIPTION: Sets initial Matlab defaults for use in the marks-based 
%%              OCR project
%%
%% CVS:
%% $Id: marks_ocr_init.m,v 1.1 2006-06-03 20:52:45 scottl Exp $
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('loading matlab defaults file\n');
ver = version;
if str2num(ver(1)) < 7
    warning('Using Matlab version < 7, some functions may not work correctly');
end
more off;  %turn off pagination
whitebg;   %inverse the figure backgrounds
close(gcf);

%add the necessary paths
addpath "cluster" "util" "display" "ocr"
