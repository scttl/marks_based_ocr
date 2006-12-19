function pos = jtag_region_finder(jtag_file, rgn_list)
% JTAG_REGION_FINDER    Returns pixel boundaries of regions of the type passed
%
%   pos = jtag_region_finder(jtag_file, region_list)
%
%   This function attempts to open the jtag file passed (a jtag file is an ascii
%   flat text file with a particular format that gives position co-ordinates of
%   hand labelled regions), and extracts the positions of all regions that
%   match any of the types passed.
%
%   jtag_file should give the path and full filename of the jtag file to be 
%   opened.  An error is returned if the file doesn't exist
%
%   region_list should be a cell array of strings which should specify the types
%   of regions to be extracted.  Possible valid choices are:
%   section_heading, subsection_heading, footer, references, bullet_item, 
%   table_label, header, authour_list, code_block, main_title, figure_label, 
%   figure, image, text, equation, footnote, figure_caption, decoration, 
%   abstract, table, graph, eq_number, editor_list, table_caption,  pg_number
%
%   pos is an nx4 matrix containing the pixel offsets (left, top, right, and
%   bottom) of each of the n regions found that correspond to one of the 
%   region types listed in region_list
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: jtag_region_finder.m,v 1.2 2006-12-19 21:41:21 scottl Exp $
%
% REVISION HISTORY
% $Log: jtag_region_finder.m,v $
% Revision 1.2  2006-12-19 21:41:21  scottl
% no change.
%
% Revision 1.1  2006/06/03 20:56:01  scottl
% Initial check-in.
%
%

% LOCAL VARS %
%%%%%%%%%%%%%%
class_str = 'class = ';
pos_str =  'pos = ';
pos = [];


% CODE START %
%%%%%%%%%%%%%%

%try and open the file passed, storing the sequence of class names and 
%associated positions (ignoring lines with '%' characters
data = textread(jtag_file, '%s', 'delimiter', '\n', 'commentstyle', 'matlab');
class_idx = strmatch(class_str, data);
pos_idx = strmatch(pos_str, data);
classes = data(class_idx);
for i=1:size(classes,1)
    classes{i} = classes{i}(length(class_str)+1:end);
end
str_pos = data(pos_idx);
for i=1:size(str_pos,1)
    P = sscanf(str_pos{i}, [pos_str, '%d %d %d %d']);
    pos = [pos; P(1), P(2), P(3), P(4)];
end

keep_idx = [];
for i=1:size(rgn_list,2)
    keep_idx = [keep_idx; strmatch(rgn_list{i}, classes)];
end
pos = pos(keep_idx,:);
