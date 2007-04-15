function out = numerize(in)
% NUMERIZE  Calculate a numerization of each of the symbol sequences passed
%
%   OUT = NUMERIZE(IN)
%   This function returns a numerized representation of each string or numeric
%   sequence passed.  Each occurence of the same symbol is given the same 
%   numerical value, and the values are numbered starting at 1 (no two 
%   different symbols may receive the same value).  ex. 'sassy' gets 
%   [1 2 1 1 3], and [4 1 10 1 9] gets [1 2 3 2 4].  Inputs should list one 
%   entry per row (or cell item).  The output will be of the same type and 
%   size as the input (with strings converted to doubles), but with values 
%   changed to vectors of numerization values.  If the input is not a cell 
%   array, each column must be the same length


% CVS INFO %
%%%%%%%%%%%%
% $Id: numerize.m,v 1.1 2007-04-15 21:29:51 scottl Exp $
%
% REVISION HISTORY
% $Log: numerize.m,v $
% Revision 1.1  2007-04-15 21:29:51  scottl
% initial check-in
%


% GLOBAL VARS %
%%%%%%%%%%%%%%%


% LOCAL VARS %
%%%%%%%%%%%%%%


% CODE START %
%%%%%%%%%%%%%%
if nargin ~= 1
    error('incorrect number of arguments specified!');
end

if iscell(in)
    out = cell(size(in));
    for ii=1:length(in)
        str = in{ii};
        nm = zeros(size(str));
        curr = 1;
        for jj=1:length(nm)
            pos = find(str(1:jj-1) == str(jj),1);
            if isempty(pos)
                nm(jj) = curr;
                curr = curr+1;
            else
                nm(jj) = nm(pos);
            end
        end
        out{ii} = nm;
    end
else
    out = zeros(size(in));
    for ii=1:size(in,1)
        curr = 1;
        for jj=1:size(in,2)
            pos = find(in(ii,1:jj-1) == in(ii,jj),1);
            if isempty(pos)
                out(ii,jj) = curr;
                curr = curr+1;
            else
                out(ii,jj) = out(ii,pos);
            end
        end
    end
end
