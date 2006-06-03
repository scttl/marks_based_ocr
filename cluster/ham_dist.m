function [match, MM1, MM2, dist] = ham_match(C1, C2, ham_thresh)
% HAM_MATCH  Calculate Hamming distance between the two cluster averages passed
%
% C1 and C2 are representations (pixel matrices) of each cluster, typically the
% representation of a particular element in the cluster.
%
% ham_thresh (if passed) is the proportion of incorrect pixels that will be 
% allowed in Hamming distance calculations (defaults to .01 if not passed).

% CVS INFO %
%%%%%%%%%%%%
% $Id: ham_dist.m,v 1.1 2006-06-03 20:55:48 scottl Exp $
%
% REVISION HISTORY
% $Log: ham_dist.m,v $
% Revision 1.1  2006-06-03 20:55:48  scottl
% Initial check-in.
%
%

% LOCAL VARS %
%%%%%%%%%%%%%%


% CODE START %
%%%%%%%%%%%%%%
if nargin >= 3
    thresh = ham_thresh;
else
    thresh = .01;  %allow up to 1% of the pixels to not match by default
end

match = false;
dist = inf;
swap = false;

S1 = size(C1);
S2 = size(C2);

MM1 = zeros(max(S1(1), S2(1)), max(S1(2), S2(2)));
MM2 = MM1;

%try various translations by shifting the smaller of the two objects around
%there are 3 cases to consider:
%1. small item fits completely inside larger item i.e.
%   small width <= large width && small height <= large height
%2. small width << large width && small height > large height
%3. small width > large width && small height << large height
if sum(S1) >= sum(S2)
    %second component smaller
    hdiff = floor((max(S1(1), S2(1)) - S1(1)) / 2);
    wdiff = floor((max(S1(2), S2(2)) - S1(2)) / 2);
    MM1(hdiff + (1:S1(1)), wdiff + (1:S1(2))) = C1;
    SM = C2;
    SH = S2(1); SW = S2(2); LH = S1(1); LW = S1(2);
else
    %first component smaller
    hdiff = floor((max(S1(1), S2(1)) - S2(1)) / 2);
    wdiff = floor((max(S1(2), S2(2)) - S2(2)) / 2);
    MM1(hdiff + (1:S2(1)), wdiff + (1:S2(2))) = C2;
    SM = C1;
    swap = true;
    SH = S1(1); SW = S1(2); LH = S2(1); LW = S2(2);
end

%just brute force it by putting starting with the top-left-corners matching
%iterating until the top-right corners match, then go to next row, and
%repeat, stopping when the bottom-right corner matches the bottom-right corner.
%This isn't the best we could do, so really this could should be improved.

%we assume that two components don't match if one is twice as big as the
%other
if (LW * LH) > 2 * (SH * SW)
    return;
end

offset = [0,0];
for i=1:(size(MM1,1) - SH + 1)
    for j=1:(size(MM1,2) - SW + 1)
        MM2 = zeros(size(MM1));
        MM2(i + (0:(SH-1)), j + (0:(SW-1))) = SM;

        %count num of non-matching pixels, see if its fewer than the threshold
        ndist = sum(sum(xor(MM1, MM2))) / (size(MM1,1) * size(MM1,2));
        if ndist < dist
            dist = ndist;
        end
        if ndist <= thresh
            match = true;
            if swap
                T = MM1;
                MM1 = MM2;
                MM2 = T;
            end
            return;
        end
    end
end
