function display_components(Comps, varargin)
%  DISPLAY_COMPONENTS  Display the connected component bounding boxes 
%
%   display_lines(COMPS, [VAR1, VAL1]...)
%
%   COMPS should be a struct like that returned in get_comps()
%
%   optional LOCAL VARS values below can be overriden specifying the name and
%   new value for the variable to be overwritten as additinoal parameters.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: display_components.m,v 1.1 2007-05-14 23:15:34 scottl Exp $
%
% REVISION HISTORY
% $Log: display_components.m,v $
% Revision 1.1  2007-05-14 23:15:34  scottl
% initial check-in
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%this will hold the list of components to draw (leave blank to draw all)
comps_idx = [];  


%set this to true to invert the original foreground and background intensities
reverse_display = false;

%show we not display the image? Useful if you just want to save it
display_cc_page = true;

%should we save connected component images?  This uses memory and takes longer
%to run
save_cc_page = false;

global MOCR_PATH;
%prefix of filename to use when saving connected component images.  If multiple
%files are passed, the page numbers are also added to the filename
img_prefix = [MOCR_PATH, '/results/conn_comps'];
%file format to use when saving the image
img_format = 'png';

comps_col = reshape([255,0,0],1,1,3);  %this is red


% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 1
    error('incorrect number of arguments specified!');
elseif nargin > 1
    process_optional_args(varargin{:});
end

if isempty(comps_idx)
    comps_idx = 1:Comps.max_comp;
end
if size(comps_idx,1) > 1
    comps_idx = comps_idx';
end

for pp = unique(Comps.pg(comps_idx))';

    M = imread(Comps.files{pp});
    %@@ideally we should examine colortype of image and act appropriately, but
    %since I've only ever dealt with binary images, this can be put off until 
    %later
    if reverse_display
        %this colormap uses the white background with black letter font
        M = label2rgb(M, zeros(Comps.max_comp,3), 'w');
    else
        %this colormap uses the black background with white letter font
        M = label2rgb(M, 'white', 'k');
    end

    idx = find(Comps.pg(comps_idx) == pp);
    idx = comps_idx(idx);
    for jj=idx
        l = Comps.pos(jj,1); t = Comps.pos(jj,2);
        r = Comps.pos(jj,3); b = Comps.pos(jj,4);
    
        M(t,l:r,:) = repmat(comps_col,1,r-l+1);
        M(t:b,l,:) = repmat(comps_col,b-t+1,1);
        M(t:b,r,:) = repmat(comps_col,b-t+1,1);
        M(b,l:r,:) = repmat(comps_col,1,r-l+1);
    end

    if display_cc_page
        imshow(M);
    end

    %save the image to disk if required.
    if save_cc_page
        fprintf('writing component page image %d to disk\n', pp);
        imwrite(M, [img_prefix, num2str(pp), '.', img_format], img_format);
    end
end

fprintf('Elapsed time: %f\n', toc);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
