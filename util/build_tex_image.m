function img = build_tex_image(str, varargin)
% BUILD_TEX_IMAGE Use TeX to create an image representation of the string passed
%
% IMG = BUILD_TEX_IAMGE(STRING, [VAR1, VAL1]...)
%
% STRING should be a single character array that is a list of characters and
% possible TeX formatting macros that will be passed to the TeX interpreter,
% and compiled into a dvi file.  Using the dvipng utility, an image
% representation of this file is created then loaded back into matlab where it
% is cropped and returned.
%
% IMG is the returned image. 
%
% NOTE: currently cannot change font etc.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: build_tex_image.m,v 1.2 2006-12-04 19:19:27 scottl Exp $
%
% REVISION HISTORY
% $Log: build_tex_image.m,v $
% Revision 1.2  2006-12-04 19:19:27  scottl
% update to fix optional arg processing, added new dvipng parameters,
% use external crop image utility function.
%
% Revision 1.1  2006-11-29 16:40:33  scottl
% split build_tex_image into its own file (so it can be called from other
% functions)
%


% LOCAL VARS %
%%%%%%%%%%%%%%
tex_file='/tmp/tmp_ocr.tex';
log_file='/tmp/tmp_ocr.log';
dvi_file='/tmp/tmp_ocr.dvi';
png_file='/tmp/tmp_ocr.png';

%how big should the image be?
png_dpi = 190;

%at what point should anti-aliased pixel values be considered on?
aa_thresh = 3;

%what to use for foreground pixels?
fg_val = 1;


% CODE START %
%%%%%%%%%%%%%%

if nargin < 1
    error('incorrect number of arguments passed');
elseif nargin > 1
    process_optional_args(varargin{:});
end

img = [];

fid = fopen(tex_file, 'w');
if fid == -1
    error('problems creating file: %s', tex_file);
end
fprintf(fid, '\\nopagenumbers\n%s\n\\end', str);
fclose(fid);

[s,w] = unix(['tex -output-directory=/tmp ', tex_file]);
if s ~= 0
    unix(['rm -f ' tex_file, ' ', dvi_file, ' ', log_file]);
    error('problem running TeX: %s', w);
end
[s,w] = unix(['dvipng -o ', png_file, ' ', ' -D ', num2str(png_dpi), ...
              ' -bg "rgb 0.0 0.0 0.0" -fg "rgb 1.0 1.0 1.0" ' dvi_file]);
if s ~= 0
    unix(['rm -f ' tex_file, ' ', dvi_file, ' ', log_file, ' ', png_file]);
    error('problem running dvipng: %s', w);
end

img = imread(png_file) > aa_thresh;
if fg_val == 0
    img = ~img;
end

%crop the image to have a tight bounding box
img = crop_image(img);

unix(['rm -f ' tex_file, ' ', dvi_file, ' ', log_file, ' ', png_file]);
