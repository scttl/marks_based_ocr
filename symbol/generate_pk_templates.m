function [bitmaps, offsets, charnames] = generate_pk_templates(font, varargin)
% GENERATE_PK_TEMPLATES  Create bitmaps of characters using the pk font passed
%
%   [bitmaps, offsets, charnames] = generate_pk_templates(font, [var1, val1]...)
%   This procedure by default uses ASCII encodings for upper and lower case
%   letters, period and space, to create bitmap images (logical arrays) of 
%   each character using the fontfile passed in font (should include the path 
%   as well).  Note that the characters can be overridden (as well as other
%   LOCAL VARS) by specifying the name of the parameter to be overridden (as a
%   string), and then its value.
%
%   font should give the full path and name of a .pk (packed font file).
%   If the font doesn't exist, an error is returned
%
%   NOTE: This makes use of the UNIX pk2bm utility, so we're limited to the
%   first 128 ASCII characters only (regardless of the machine's encoding),
%   thus any ligatures, or multiple character symbols will not be generated
%   correctly (only the first character is taken)
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: generate_pk_templates.m,v 1.4 2007-04-16 03:30:02 scottl Exp $
%
% REVISION HISTORY
% $Log: generate_pk_templates.m,v $
% Revision 1.4  2007-04-16 03:30:02  scottl
% small change to declaration of logical array
%
% Revision 1.3  2006-11-22 17:01:46  scottl
% small update to warn when trying to recognize multi-char symbols
%
% Revision 1.2  2006-11-07 02:49:06  scottl
% bugfix for handling single quote character.
%
% Revision 1.1  2006-10-29 17:06:21  scottl
% initial revision
%
% Revision 1.3  2006/07/28 22:19:03  scottl
% create tight bounding boxes, and calculate vertical offsets instead of
% resizing.
%
% Revision 1.2  2006/07/05 01:01:54  scottl
% small spelling fixups.
%
% Revision 1.1  2006/06/10 21:01:36  scottl
% Initial revision.
%

% LOCAL VARIABLES %
%%%%%%%%%%%%%%%%%%%

%default characters to generate templates of
charnames=mat2cell(char([32,46,65:90,97:122])');  %A-Z, a-z, ., ' '

pat = ['.* height : (\d+)\n', ... %get the value of the height
       '.* width : (\d+)\n', ... %the value of the width
       '.* xoff : (\d+|-\d+)\n', ... %x offset (can be +ve or -ve)
       '.* yoff : (\d+|-\d+)\n\n', ... %y offset (again +ve or -ve)
       '(.*)'];   %this will initially hold the character bitmap


% CODE START %
%%%%%%%%%%%%%%
if nargin < 1
    error('incorrect number of arguments specified!');
elseif nargin > 1
    process_optional_args(varargin{:});
end

[fid, msg] = fopen(font);
if fid == -1
    error('unable to find/open the font file: %s.\nReason: %s', font, msg);
end
fclose(fid);

numchars=length(charnames);

for ii=1:numchars
    thischar=charnames{ii};
    if length(thischar) > 1
        warning('MBOCR:MulticharSym', 'Multi-char templates not supported');
        thischar = thischar(1);
    end
    fprintf(1,'Char %c (%d/%d)\n',thischar,ii,numchars);

    if strcmp(thischar, '''')
        %since the character is a single quote, we must escape it differently
        %for the shell
        cmd = ['pk2bm -b -c"''" ', font];
    else
        cmd = ['pk2bm -b -c''', thischar, ''' ', font];
    end
    [s,w] = unix(cmd);
    if s ~= 0
        error('problem running pk2bm: %s', w);
    end
    % if run successfully, pk2bm output is something like the following:
    %    
    %character : 97 (a)
    %   height : 18
    %    width : 18
    %     xoff : -2
    %     yoff : 17
    %
    %  ...*******........
    %  ..**.....***......
    %  .****.....***.....
    %  .****......***....
    %  ..**.......***....
    %  ...........***....
    %  ...........***....
    %  .......*******....
    %  ....****...***....
    %  ..***......***....
    %  .***.......***....
    %  .**........***....
    %  ***........***...*
    %  ***........***...*
    %  ***........***...*
    %  .**.......****...*
    %  ..**.....*..***.*.
    %  ...******....***..

    res = regexp(w, pat, 'tokens');

    vals(ii,1) = str2num(res{1}{3});
    vals(ii,2) = str2num(res{1}{4});
    
    vals(ii,3) = str2num(res{1}{1});
    vals(ii,4) = str2num(res{1}{2});

    bitmaps{ii} = zeros(vals(ii,3), vals(ii,4));

    %we convert the string into a matrix, then trim the 2 leading spaces, and
    %trailing single newline character from the matrix, then convert this to
    %a numeric array
    on_idx = find(strtrim(reshape(res{1}{5}, vals(ii,4)+3, ...
             vals(ii,3))') == '*');
    bitmaps{ii}(on_idx) = 1;
end

%the baseline offsets are calculated as the difference between the height of
%the character and its corresponding yoff.  Note that we must subtract 1 from
%this value because of the way the offsets are done in metafont. 
offsets = vals(:,3) - vals(:,2) - 1;

%must manually fixup the space character (if it appears) since it isn't
%completely blank
space_idx = find(charnames == ' ');
if space_idx
    bitmaps{space_idx} = false(size(bitmaps{space_idx}));
    offsets(space_idx) = 0;
end
