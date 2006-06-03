function [match, MM1, MM2, dist] = euc_match(C1, C2, euc_thresh)
%  EUC_MATCH  Determine the Euclidean distance between two cluster averages
%
% [match, MM1, MM2, dist] = euc_match(C1, C2, euc_thresh)
%
% C1 and C2 are image intensity representations (typically the average of
% all items in that cluster).
%
% euc_thresh is the maximum length of the Euclidean distance allowed 
% between the average pixel intensities to be considered as a match

% CVS INFO %
%%%%%%%%%%%%
% $Id: euc_dist.m,v 1.1 2006-06-03 20:55:47 scottl Exp $
%
% REVISION HISTORY
% $Log: euc_dist.m,v $
% Revision 1.1  2006-06-03 20:55:47  scottl
% Initial check-in.
%
%

match = false;
S1 = size(C1);
S2 = size(C2);
MM1 = zeros(max(S1(1), S2(1)), max(S1(2), S2(2)));
MM2 = MM1;
MM1(1:S1(1), 1:S1(2)) = C1;
MM2(1:S2(1), 1:S2(2)) = C2;

dist = sqrt(sum(sum((MM1 - MM2).^2))) / min((S1(1) * S1(2)), (S2(1) * S2(2)));

if dist <= euc_thresh
    match = true;
end
