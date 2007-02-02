function d = assign_density(imgs, varargin)
% ASSIGN_DENSITY   Computer normalized densities for the images passed
%
%   D = ASSIGN_DENSITY(IMGS, [VAR1, VAL1]...)
%
%   IMGS should be a cell array of pixel intensity images to assign a density
%   to.  Each density computed is normalized (relative to others) so that the
%   entire collection has mean density of 0, and standard deviation 1.
%
%   IMGS could also be a char array which represents the symbols to draw.  A
%   default font and encoding are used, and a tight bounding box is created for
%   each one.
%
%   VAR1 and VAL1 are optional, and can be used to override the default values
%   for the LOCAL VARS defined below.  Each VAR1 should be a string giving the
%   exact name of the variable to override.  Each VAL1 should be the updated
%   value to use for that variable.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: assign_density.m,v 1.1 2007-02-02 05:51:56 scottl Exp $
%
% REVISION HISTORY
% $Log: assign_density.m,v $
% Revision 1.1  2007-02-02 05:51:56  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%parameters that control how we generate symbol images of a string is passed
img_font = 'helvetica';
img_font_sz = '36';

% CODE START %
%%%%%%%%%%%%%%

%process input arguments
if nargin < 1
    error('must pass list of images');
elseif nargin > 1
    process_optional_args(varargin{:});
end

num = length(imgs);

if ~iscell(imgs)
    %must convert to images with bounding boxes etc.  This assumes
    %@@@ImageMagick's convert is installed.  
    tmp = cell(num,1);
    for ii=1:num
        fprintf('creating image of symbol %s\n', imgs(ii));
        if imgs(ii) == 32 || imgs(ii) == 39 || imgs(ii) == 34
            %have to handle spaces and ' separately since convert doesn't like 
            %them
            tmp{ii} = zeros(12);
            continue;
        end
        fid = fopen('/tmp/tmp.txt','w');
        if fid == -1
            error('problem opening file');
        end
        fwrite(fid, imgs(ii), 'char');
        fclose(fid);
        cmd = ['convert -type bilevel +matte -font ', img_font, ...
            ' -pointsize ', img_font_sz, ' label:@/tmp/tmp.txt /tmp/tmp.tif'];
        [s,w] = unix(cmd);
        if s~= 0
            error('problem running %s', cmd);
        end
        tmp{ii} = ~imread('/tmp/tmp.tif'); 
        %crop the bounding box to be tight
        p = get_comp_bb(tmp{ii},1);
        tmp{ii} = tmp{ii}(p(2):p(4), p(1):p(3));
    end
    delete '/tmp/tmp.txt';
    delete '/tmp/tmp.tif';
    imgs = tmp;
end


d = zeros(length(imgs),1);
for ii=1:num
    d(ii) = double(sum(imgs{ii}(:))) / prod(size(imgs{ii}));
end

%now normalize the images
d = (d - mean(d)) ./ std(d);
