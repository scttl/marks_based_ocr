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
% $Id: build_tex_image.m,v 1.1 2006-11-29 16:40:33 scottl Exp $
%
% REVISION HISTORY
% $Log: build_tex_image.m,v $
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
png_dpi='300';


% CODE START %
%%%%%%%%%%%%%%

if nargin < 1
    error('incorrect number of arguments passed');
elseif nargin > 1
    process_option_args(varargin{:});
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

[s,w] = unix(['dvipng -o ', png_file, ' ', ' -D ', png_dpi, ' ' dvi_file]);
if s ~= 0
    unix(['rm -f ' tex_file, ' ', dvi_file, ' ', log_file, ' ', png_file]);
    error('problem running dvipng: %s', w);
end

img = (imread(png_file) > 0);
[r,c] = find(img == 1, 1, 'first');
img = img(:,c:end);
[r,c] = find(img == 1, 1, 'last');
img = img(:,1:c);
unix(['rm -f ' tex_file, ' ', dvi_file, ' ', log_file, ' ', png_file]);
