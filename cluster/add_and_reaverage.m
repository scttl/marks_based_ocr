function C = add_and_reaverage(C, C2, M, M2)
%  ADD_AND_REAVERAGE  Add the 2nd cluster contents to the first
%
%  C = add_and_reaverage(C, C2, M, M2)
%
%  C and C2 should be cluster items (see cluster_components.m for more info),
%
%  M and M2 should be identically sized cluster averages (like that returned by
%  euc_match etc.

% CVS INFO %
%%%%%%%%%%%%
% $Id: add_and_reaverage.m,v 1.1 2006-06-03 20:55:47 scottl Exp $
%
% REVISION HISTORY
% $Log: add_and_reaverage.m,v $
% Revision 1.1  2006-06-03 20:55:47  scottl
% Initial check-in.
%
%

% CODE START %
%%%%%%%%%%%%%%
C.comp = [C.comp; C2.comp];
C.pos = [C.pos; C2.pos];
C.nb = [C.nb; C2.nb];
C.pg = [C.pg; C2.pg];
old_num = C.num;
C.num = C.num + C2.num;

% note that the averages may be different sizes, so we have that taken into
% account during the distance call
[cr, cc] = size(C.avg);
[c2r, c2c] = size(C2.avg);
C.avg = (old_num/C.num .* M) + (C2.num/C.num .* M2);

