function [match, MA, MB, dist] = hausdorff_match(A, B, haus_thresh)
%   HAUSDORFF_MATCH  Calculate the Hausdorff distance between two cluster avgs.
%
% [match, MA, MB, dist] = hausdorff_match(A, B, haus_thresh)
%
% A and B are iamge intensity representations (typically the average of
% all items in that cluster).
%
% haus_thresh is the maximum length of the euclidian distance allowed 
% between the average pixel intensities to be considered as a match

% CVS INFO %
%%%%%%%%%%%%
% $Id: hausdorff_dist.m,v 1.1 2006-06-03 20:55:48 scottl Exp $
%
% REVISION HISTORY
% $Log: hausdorff_dist.m,v $
% Revision 1.1  2006-06-03 20:55:48  scottl
% Initial check-in.
%
%


% LOCAL VARS %
%%%%%%%%%%%%%%
%percent = .95;  %used for 'soft' distance calculation (not yet implemented)
match = false;


% CODE START %
%%%%%%%%%%%%%%

%first center A and B
%this can be done by calling conv_euc_match (though its not efficient)
[match, MA, MB, dist] = conv_euc_match(A, B, thresh);
MMA = MA >= 0.5;  %binarize these matrices
MMB = MB >= 0.5;

%get the minimal euclidian distance from each point in A to each point in B
a_idx = find(MMA ~= 0);
BD = bwdist(MMB);
atob_dist = BD(a_idx);
%now take the smallest maximum that includes percent of the points in B
%@@ to complete.  Currently taking 100%
a_delta = max(atob_dist);

%repeat from each point in B to each point in A
b_idx = find(MMB ~= 0);
AD = bwdist(MMA);
btoa_dist = AD(b_idx);
%now take the smallest maximum that includes percent of the points in A
%@@ to complete.  Currently taking 100%
b_delta = max(btoa_dist);

%take the average of the two deltas as the distance measure
dist = (a_delta + b_delta)/2;

if dist <= haus_thresh
    match = true;
end
