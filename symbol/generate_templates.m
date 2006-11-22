function bitmaps = generate_templates(fontname, sn, varargin)
% GENERATE_TEMPLATES  Create bitmaps of characters using the fontname passed
%
%   bitmaps = generate_templates(fontname, symnames, [var1, val1]...)
%   This procedure makes use of ImageMagick's convert utility to generate
%   bitmap images of a particular font in a given size.
%
%   fontname should be the name of a font that ImageMagick knows about.
%   If ImageMagick isn't present, or doesn't know about the fontname, an error
%   is returned.
%
%   symnames should be a cell array of symbols to generate bitmaps of.
%
%   NOTE: this program makes use of the ImageMagick 'convert' utility 
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: generate_templates.m,v 1.6 2006-11-22 17:00:56 scottl Exp $
%
% REVISION HISTORY
% $Log: generate_templates.m,v $
% Revision 1.6  2006-11-22 17:00:56  scottl
% updates to fix convert's output of specific characters.
%
% Revision 1.5  2006-11-14 22:51:14  scottl
% changed to a cell array of strings instead of chars.
%
% Revision 1.4  2006-10-29 17:05:49  scottl
% moved and renamed dofont.m to generate_templates.m.  Rewritten to make
% use of imageMagick's convert utility.
%


% LOCAL VARIABLES %
%%%%%%%%%%%%%%%%%%%
ptsize = 32; %default pointsize

%where should the temporary image be generated (will be deleted)
tmp_img = '/tmp/tmp_template.png';

%should we remove the dots above lowercase i, and j characters?
strip_dots = true;

%this temp file is used to write the value of symbols like '@' that convert 
%has problems displaying from the command line.
tmp_sym_file = '/tmp/tmp_sym.txt';


% CODE START %
%%%%%%%%%%%%%%
if nargin < 2
    error('incorrect number of arguments specified!');
elseif nargin > 2
    process_optional_args(varargin{:});
end

numsyms=length(sn);

for ii=1:numsyms
    fprintf(1,'Symbol %s (%d/%d)\n',sn{ii},ii,numsyms);

    if strcmp(sn{ii}, ' ')
        %space character.  Can't render a bitmap of this using imagemagick
        %so we'll just create one by hand, using ptsize as a loose guide
        dims = floor(ptsize/2);
        bitmaps{ii,:} = zeros(dims,dims);
        continue;
    elseif strcmp(sn{ii}, '''')
        %single quotes have to be quoted differently to prevent shell probs
        cmd = ['convert -trim -font ', fontname, ' -pointsize ', ...
               num2str(ptsize), ' label:"\''" ', tmp_img];
        %@@note: for CMR10 this produces the " character.
    elseif strcmp(sn{ii}, '"')
        %double quotes need the '\' as a prefix for CMR10
        cmd = ['convert -trim -font ', fontname, ' -pointsize ', ...
               num2str(ptsize), ' label:''\"'' ', tmp_img];
    elseif ~isempty(strfind('@{}', sn{ii}))
        %the '@' symbol is special to convert
        fid = fopen(tmp_sym_file, 'w');
        if fid == -1
            error('problems creating temporary sym file: %s', tmp_sym_file);
        end
        fprintf(fid, '%s\n', sn{ii});
        fclose(fid);
        cmd = ['convert -trim -font ', fontname, ' -pointsize ', ...
                     num2str(ptsize), ' label:@', tmp_sym_file, ' ', tmp_img];
    else
        cmd = ['convert -trim -font ', fontname, ' -pointsize ', ...
                     num2str(ptsize), ' label:''', sn{ii}, ''' ', tmp_img];
    end
    [s,w] = unix(cmd);

    if s ~= 0
        unix(['rm -f ' tmp_img]);
        error('problem running ImageMagick: %s', w);
    end

    bitmaps{ii,:} = ~im2bw(imread(tmp_img),0.5);  %convert to binary

    unix(['rm -f ' tmp_img, ' ', tmp_sym_file]);

    if strip_dots && (strcmp(sn{ii}, 'i') || strcmp(sn{ii}, 'j'))
        %strip the dots from these lowercase letters to match what we do during
        %connected components finding
        [lbl_img, count] = bwlabel(bitmaps{ii,:});
        if count == 2
            min_label = 1;
            min_count = sum(lbl_img == 1);
            if sum(lbl_img == 2) <= min_count
                min_label = 2;
            end
            bitmaps{ii}(lbl_img == min_label) = 0;
        elseif count > 2
            warning('MBOCR:TooManyComp', 'More than 2 comps found!\n');
        end
    end
end
