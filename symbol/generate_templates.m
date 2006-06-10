function [bitmaps,charnames]=dofont(fontname, cn)
% DOFONT  Create bitmaps of characters using the pk font passed
%
%   [bitmaps,charnames]=dofont(fontname, [character_list])
%   This procedure takes the list of ASCII characters passed in via 
%   character_list, or uses upper and lower case, letters, period and space if
%   not give, to create bitmap images (logical arrays) of each character using
%   the fontfile passed in fontname (should include the path as well)
%
%   fontname should give the full path and name of a .pk (packed font file).
%   If the font doesn't exist, an error is returned
%
%   character_list is optional and if specified should be a character array
%   listing the individal ASCII characters to generate bitmaps off.  If not 
%   given, we generate bitmaps of: A-Z, a-z, ' ', '.'
%
%   See Also: pk2bm UNIX utility
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: generate_templates.m,v 1.1 2006-06-10 21:01:36 scottl Exp $
%
% REVISION HISTORY
% $Log: generate_templates.m,v $
% Revision 1.1  2006-06-10 21:01:36  scottl
% Initial revision.
%

% LOCAL VARIABLES %
%%%%%%%%%%%%%%%%%%%
charnames=char([32,46,65:90,97:122]);

pat = ['.* height : (\d+)\n', ... %get the value of the height
       '.* width : (\d+)\n', ... %the value of the width
       '.* xoff : (\d+|-\d+)\n', ... %x offset (can be +ve or -ve)
       '.* yoff : (\d+|-\d+)\n\n', ... %y offset (again +ve or -ve)
       '(.*)'];   %this will initially hold the character bitmap

% CODE START %
%%%%%%%%%%%%%%
if nargin < 1 || nargin > 2
    error('incorrect number of arguments specified!');
elseif nargin == 2
    charnames = cn;
end

[fid, msg] = fopen(fontname);
if fid == -1
    error('unable to find/open the font file: %s.\nReason: %s', fontname, msg);
end
fclose(fid);

numchars=length(charnames);

for ii=1:numchars
    thischar=charnames(ii);
    fprintf(1,'Char %c (%d/%d)\n',thischar,ii,numchars);

    [s,w] = unix(['pk2bm -b -c''', thischar, ''' ', fontname]);
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

    offsets(ii,1) = str2num(res{1}{3});
    offsets(ii,2) = str2num(res{1}{4});
    
    offsets(ii,3) = str2num(res{1}{1});
    offsets(ii,4) = str2num(res{1}{2});

    bm{ii} = zeros(offsets(ii,3), offsets(ii,4));

    %we convert the string into a matrix, then trim the 2 leading spaces, and
    %training single newline character from the matrix, then convert this to
    %a numeric array
    on_idx = find(strtrim(reshape(res{1}{5}, offsets(ii,4)+3, ...
             offsets(ii,3))') == '*');
    bm{ii}(on_idx) = 1;

end

tops=offsets(:,2)+1;
bots=tops-offsets(:,3);
maxtop=max(tops);
minbot=min(bots);

for ii=1:numchars
    ww=size(bm{ii},2);
    bitmaps{ii}=logical([zeros(maxtop-tops(ii),ww); ...
                         bm{ii}; ...
                         zeros(bots(ii)-minbot,ww)]);
end

%must manually fixup the space character (if it appears) since it isn't
%completely blank
space_idx = find(charnames == ' ');
if space_idx
    bitmaps{space_idx}=logical(zeros(size(bitmaps{space_idx})));
end
