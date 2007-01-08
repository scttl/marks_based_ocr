function [pos,keep_idx] = jtag_region_finder(jtag_file, rgn_list)
% JTAG_REGION_FINDER    Returns pixel boundaries of regions of the type passed
%
%   [pos,reg_num] = jtag_region_finder(jtag_file, region_list)
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
%   reg_num is a vector of length n, and lists what region number each of the
%   kept regions is.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: jtag_region_finder.m,v 1.5 2007-01-08 22:10:37 scottl Exp $
%
% REVISION HISTORY
% $Log: jtag_region_finder.m,v $
% Revision 1.5  2007-01-08 22:10:37  scottl
% add the region number to the list of output parameters.
%
% Revision 1.4  2007-01-05 17:05:21  scottl
% small fix to ensure regions are returned in their original order
%
% Revision 1.3  2007-01-02 19:21:52  scottl
% made class and position matching more robust, removed depenedence on
% iteration over lines.
%
% Revision 1.2  2006-12-19 21:41:21  scottl
% no change.
%
% Revision 1.1  2006/06/03 20:56:01  scottl
% Initial check-in.
%
%

% LOCAL VARS %
%%%%%%%%%%%%%%
class_str = 'class ';
class_pat = 'class\s*=\s*(\w*)\s*';
class_rep = '$1';
pos_str =  'pos ';
pos_pat = 'pos\s*=\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*';
pos_rep = '$1 $2 $3 $4';
pos = [];


% CODE START %
%%%%%%%%%%%%%%

%try and open the file passed, storing the sequence of class names and 
%associated positions (ignoring lines with '%' characters
data = textread(jtag_file, '%s', 'delimiter', '\n', 'commentstyle', 'matlab');
class_idx = strmatch(class_str, data);
pos_idx = strmatch(pos_str, data);
classes = data(class_idx);
classes = regexprep(classes, class_pat, class_rep);
str_pos = data(pos_idx);
str_pos = regexprep(str_pos, pos_pat, pos_rep);
pos = str2num(char(str_pos));

keep_idx = [];
for i=1:size(rgn_list,2)
    keep_idx = [keep_idx; strmatch(rgn_list{i}, classes)];
end
%ensure we return the positions in the same order they were given
keep_idx = sort(keep_idx);
pos = pos(keep_idx,:);
