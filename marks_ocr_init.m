%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% FILE: marks_ocr_init.m
%%
%% DESCRIPTION: Sets initial Matlab defaults for use in the marks-based 
%%              OCR project
%%
%% NOTE: Be sure to set the global MOCR_PATH variable in your startup.m file
%%       to point at the path where this file resides, then executre this
%%       file.
%% CVS:
%% $Id: marks_ocr_init.m,v 1.4 2006-09-22 18:02:11 scottl Exp $
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Initializing Marks-Based OCR\n');
ver = version;
if str2num(ver(1)) < 7
    warning('Using Matlab version < 7, some functions may not work correctly');
end

more off;  %turn off pagination
whitebg;   %inverse the figure backgrounds
close(gcf);

%add the necessary paths
global MOCR_PATH;
if isempty(MOCR_PATH)
    error('Must set global MOCR_PATH to the install directory.  See README');
end
addpath([MOCR_PATH, '/cluster']);
addpath([MOCR_PATH, '/component']);
addpath([MOCR_PATH, '/line']);
addpath([MOCR_PATH, '/util']);
addpath([MOCR_PATH, '/display']);
addpath([MOCR_PATH, '/ocr']);
addpath([MOCR_PATH, '/lang_model']);
