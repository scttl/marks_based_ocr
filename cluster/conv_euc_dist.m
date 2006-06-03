function [match, c1_avg, c2_avg, dist] = conv_euc_match(C1, C2, euc_thresh)
%  CONV_EUC_MATCH  Determine the Euclidean distance by convolving the clusters
%
% [match, MM1, MM2, dist] = conv_euc_match(C1, C2, euc_thresh)
%
% C1 and C2 are iamge intensity representations (typically the average of
% all items in that cluster).
%
% euc_thresh is the maximum length of the euclidian distance allowed 
% between the average pixel intensities to be considered as a match

% CVS INFO %
%%%%%%%%%%%%
% $Id: conv_euc_dist.m,v 1.1 2006-06-03 20:55:47 scottl Exp $
%
% REVISION HISTORY
% $Log: conv_euc_dist.m,v $
% Revision 1.1  2006-06-03 20:55:47  scottl
% Initial check-in.
%
%

match = false;
S1 = sum(size(C1));
S2 = sum(size(C2));

%we want to pass the smaller of the two cluster averages across the larger
if S1 >= S2
    %C1 larger than (or equal) to C2
    sm_avg = C2;
    lr_avg = C1;
else
    %C2 larger than C1
    sm_avg = C1;
    lr_avg = C2;
end

%determine the best overlapping point
lr_sum = sum(sum(lr_avg));
[min_row_sc, min_row_pos] = min(abs(filter2(single(lr_avg), ...
                            single(sm_avg)) - lr_sum));
[score, min_col_pos] = min(min_row_sc);

%0 pad the overlapping clusters so that they are the same size
[lr_rows, lr_cols] = size(lr_avg); nlr_rows = lr_rows; nlr_cols = lr_cols;
[sm_rows, sm_cols] = size(sm_avg); nsm_rows = sm_rows; nsm_cols = sm_cols;

ldiff = min_col_pos - ceil(lr_cols/2);
if ldiff > 0
    % add left columns to lr
    lr_avg = [zeros(nlr_rows, ldiff), lr_avg];
    nlr_cols = nlr_cols + ldiff;
elseif ldiff < 0
    % add left columns to sm
    sm_avg = [zeros(nsm_rows, -ldiff), sm_avg];
    nsm_cols = nsm_cols - ldiff;
end
tdiff = min_row_pos(min_col_pos) - ceil(lr_rows/2);
if tdiff > 0
    %add top rows to lr
    lr_avg = [zeros(tdiff, nlr_cols); lr_avg];
    nlr_rows = nlr_rows + tdiff;
elseif tdiff < 0
    %add top rows to sm
    sm_avg = [zeros(- tdiff, nsm_cols); sm_avg];
    nsm_rows = nsm_rows - tdiff;
end
rdiff = sm_cols - (min_col_pos + floor(lr_cols/2));
if rdiff > 0
    %add lright columns to lr
    lr_avg = [lr_avg, zeros(nlr_rows, rdiff)];
    nlr_cols = nlr_cols + rdiff;
elseif rdiff < 0
     %add lright columns to sm
     sm_avg = [sm_avg, zeros(nsm_rows, -rdiff)];
     nsm_cols = nsm_cols - rdiff;
end
bdiff = sm_rows - (min_row_pos(min_col_pos) + floor(lr_rows/2));
if bdiff > 0
     %add bottom rows to lr
     lr_avg = [lr_avg; zeros(bdiff, nlr_cols)];
     nlr_rows = nlr_rows + bdiff;
elseif bdiff < 0
     %add bottom rows to sm
     sm_avg = [sm_avg; zeros(-bdiff, nsm_cols)];
     nsm_rows = nsm_rows - bdiff;
end

dist = sqrt(sum(sum((lr_avg - sm_avg).^2))) / (sm_rows * sm_cols);

if dist <= euc_thresh
    match = true;
end

if S1 >= S2
    c1_avg = lr_avg;
    c2_avg = sm_avg;
else
    c1_avg = sm_avg;
    c2_avg = lr_avg;
end
