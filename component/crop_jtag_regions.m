function [new_img,reg_num,pos] = crop_jtag_regions(jtag_file, img, ...
                                 keep_region_list, varargin)
% CROP_JTAG_REGIONS  Use a JTAG file to crop components due to their region type
%
% [CROPPED_IMG, REG_NUM, REG_POS] = crop_jtag_regions(FILENAME, IN_IMG,
%                                   JTAG_CROP_REGIONS, [VAR1, VAL1]...)
%
% This function looks for a JTAG named FILENAME, and if it finds one, it 
% augments the components found in IN_IMG, removing those that are positioned 
% outside regions belonging to the JTAG_CROP_REGIONS list.
%
% The updated component image is returned in CROPPED_IMG, and the region
% numbers in REG_NUM, and associated region positions in REG_POS
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: crop_jtag_regions.m,v 1.5 2007-02-01 18:01:23 scottl Exp $
%
% REVISION HISTORY
% $Log: crop_jtag_regions.m,v $
% Revision 1.5  2007-02-01 18:01:23  scottl
% fix to handle the case where region info is selected for use, but no
% jtag file can be found.
%
% Revision 1.4  2007-01-08 22:04:20  scottl
% added region number to the list of returned parameters.
%
% Revision 1.3  2007-01-05 17:09:21  scottl
% return region information also.
%
% Revision 1.2  2007-01-02 19:23:38  scottl
% modified code so that regions are kept instead of removed (thus background
% is clean)
%
% Revision 1.1  2006-09-20 21:45:33  scottl
% initial check-in.
%

% LOCAL VARS %
%%%%%%%%%%%%%%
bg_val = 0;


% CODE START %
%%%%%%%%%%%%%%
if nargin < 3
    error('incorrect number of arguments specified!');
elseif nargin > 3
    process_optional_args(varargin);
end

%(we assume the jtag file is in the same directory as the image file,
% and has the same name, but with the extension .jtag instead)
new_img = img;
if exist(jtag_file, 'file');
    [pos,reg_num] = jtag_region_finder(jtag_file, keep_region_list);
    new_img(:) = bg_val;
    for ii=1:size(pos,1)
        new_img(pos(ii,2):pos(ii,4), pos(ii,1):pos(ii,3)) = ...
                img(pos(ii,2):pos(ii,4), pos(ii,1):pos(ii,3));
    end
else
    reg_num = 1;
    pos = [1 1 size(img)];
end
