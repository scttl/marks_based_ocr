function img = crop_jtag_regions(jtag_file, img, region_list, varargin)
% CROP_JTAG_REGIONS  Use a JTAG file to crop components due to their region type
%
% CROPPED_IMG = crop_jtag_regions(FILENAME, IN_IMG, JTAG_CROP_REGIONS)
%
% This function looks for a JTAG named FILENAME, and if it finds one, it 
% augments the components found in IN_IMG, removing those that are positioned 
% inside regions belonging to the JTAG_CROP_REGIONS list.
%
% The updated component image is returned in CROPPED_IMG
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: crop_jtag_regions.m,v 1.1 2006-09-20 21:45:33 scottl Exp $
%
% REVISION HISTORY
% $Log: crop_jtag_regions.m,v $
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
if exist(jtag_file, 'file');
    pos = jtag_region_finder(jtag_file, region_list);
    for ii=1:size(pos,1)
        img(pos(ii,2):pos(ii,4), pos(ii,1):pos(ii,3)) = bg_val;
    end
end