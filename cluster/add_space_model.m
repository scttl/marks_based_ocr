function [Clust, Comps] = add_space_model(Clust, Comps, varargin)
% ADD_SPACE_MODEL Attempt to determine transitions for the space character
%
%   [Clust, Comps] = ADD_SPACE_MODEL(Clust, Comps, [VAR1, VAL1]...)
%
%   This will create a single new cluster representing 'space' character
%   transitions, and create many new components representing these empty
%   components.  The width of a space, is (by default) inferred manually, by
%   looking at the histogram of distances between right neighbours.
%   Essentially this should be bimodal, with the first peak representing
%   typical inter-character spacing lengths, and the second (smaller) peak
%   representing typical inter-word spacing lengths.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: add_space_model.m,v 1.10 2007-05-14 23:13:19 scottl Exp $
%
% REVISION HISTORY
% $Log: add_space_model.m,v $
% Revision 1.10  2007-05-14 23:13:19  scottl
% add space distances to the clust struct
%
% Revision 1.9  2007-04-10 15:46:20  scottl
% working implementation of Hunag space model implemented.
%
% Revision 1.8  2007-02-05 20:34:19  scottl
% added density field estimation.
%
% Revision 1.7  2007-02-01 18:10:12  scottl
% added new class field that is assigned based on offset information
%
% Revision 1.6  2007-01-08 22:02:48  scottl
% small fix for pages that don't contain enough spaces to estimate
% modes properly.
%
% Revision 1.5  2006-12-17 20:11:03  scottl
% calculate spaces as gap between neighbouring components who lie at least
% partially between the baseline and x-height.
%
% Revision 1.4  2006-11-22 17:02:54  scottl
% Updated to reflect changes in truth_label field
%
% Revision 1.3  2006-11-13 17:56:12  scottl
% small fix to remove dependence on timers.
%
% Revision 1.2  2006-11-07 02:51:05  scottl
% estimate space width as where the 2nd mode start to rise (instead of the
% peak).  This gives better estimates in practice.  Also add the truth
% label if required.
%
% Revision 1.1  2006-10-29 17:26:23  scottl
% initial check-in.
%


% LOCAL VARS %
%%%%%%%%%%%%%%
bg_val = 0;

%if left blank, this width of a space character will attempt to be guessed.  To
%use a specific value, specify that instead.
space_width = [];

%if left blank, the height of a space character will match the width.  Can be
%overridden by specifying a value for the height below
space_height = [];

%if estimating space width, should we use the Poisson mixture model as described
%in the Huang paper?  If set to false, we will use a weak mode hunting method.
use_poisson_mix_model = true;

%if using the Poisson mixture model, we need initial guesses of the parameters
%that we attempt to infer
init_l1 = log(4);  %lambda for the first Poisson
init_l2 = log(15); %lambda for the second Poisson
%note, initial space width is estimated below.
init_params = [init_l1 init_l2 init_l2]';

%if using the Poisson mixture and inferring space width via minimize(), how
%many iterations (if -ve), or line searches (if +ve) should we allow
num_iters = -1000;

%if using the Poisson mixture, how many restarts should we attempt before 
%giving up?
max_restarts = 10;

%when attempting to infer the width of a space character, what is the maximum
%possible number of pixels wide we should check or allow
max_wordspace = 70;


% CODE START %
%%%%%%%%%%%%%%

if nargin < 2
    error('incorrect number of arguments passed');
elseif nargin > 2
    process_optional_args(varargin{:});
end

%get the distances between left-right neighbours
idx = find(Comps.nb(:,3) ~= 0);
trans_dist = double(Comps.nb_dist(idx,3));
min_char_space = mode(trans_dist);

    
if isempty(space_width)
    %calculate space_width manually.
    %coarsely estimate spaces by mode hunting in the space width histogram.
    %First peak should be interchar spacing but second should be interword 
    %(approximately)
    space_counts = hist(trans_dist, 1:max_wordspace);
    peaks_found = 0;
    increasing = false;
    for ii=2:max_wordspace
        if space_counts(ii) > space_counts(ii-1)
            if ~increasing
                increasing = true;
                start_pt = ii;
            end
        elseif increasing
            increasing = false;
            peaks_found = peaks_found + 1;
        end
        if peaks_found == 2
            space_width = ii-1 - ceil((ii-1 - start_pt)/1);
            break;
        end
    end
    if isempty(space_width)
        %this can happen if there aren't enough spaces to create 2 peaks (ie.
        %only 1 word of components.  Just estimate space width at some value 
        %larger than the largest component
        space_width = max(trans_dist) + 1;
    end
    if use_poisson_mix_model
        %determine space width threshold using a mixture of 2 Poisson's model
        %must trim spaces that are smaller are than 1 or are very large (since 
        %factorial taken)
        num_restarts = 0;
        init_params(3) = log(space_width);
        [T,ft,ii] = minimize(init_params, 'poisson_mix', num_iters, trans_dist);
        est_width = round(exp(T(3)));
        while num_restarts < max_restarts && (est_width < min_char_space || ...
              est_width > max_wordspace)
            fprintf('restart %d, l1= %.2f l2= %.2f c = %f\n', num_restarts, ...
                    exp(T(1)), exp(T(2)), est_width);
            num_restarts = num_restarts+1;
            if rem(est_width,2)
                init_params(1) = init_l1 + num_restarts*rand;
                init_params(2) = init_l2 + 2*num_restarts*rand;
                init_params(3) = log(space_width+num_restarts);
            else
                init_params(1) = init_l1 - num_restarts*rand;
                init_params(2) = init_l2 - 2*num_restarts*rand;
                if init_params(1) < 0
                    init_paras(1) = -init_params(1);
                end
                if init_params(2) < 0
                    init_paras(2) = -init_params(2);
                end
                init_params(3) = max(log(space_width-num_restarts),1);
            end
            [T,ft,ii] = minimize(init_params, 'poisson_mix', num_iters, ...
                        trans_dist);
            est_width = round(exp(T(3)));
        end
        if num_restarts >= max_restarts
            fprintf('unable to find reasonable space width\n');
            est_width = space_width;
        end
        Clust.space_lambda1 = exp(T(1));
        Clust.space_lambda2 = exp(T(2));
        space_width = est_width;
        Clust.space_width = space_width;
        fprintf('minimum found after %d iters of the %d restart.\n', ii, ...
                num_restarts);
    end
    fprintf('estimated space width at %d pixels\n', space_width);
end

if isempty(space_height)
    %set the height equal to the width (by default)
    space_height = space_width;
end

%we will create spaces in the large gaps between neighbouring components that
%lie at least partially between the x-height line and the baseline 
%(assuming line information has been found).
num_spaces = floor(trans_dist / space_width);
space_idx = find(num_spaces >= 1);
if Comps.found_lines && ~isempty(space_idx)
    heights = Comps.pos(idx(space_idx),4) - Comps.pos(idx(space_idx),2) + 1;
    asc = Comps.ascender_off(idx(space_idx));
    desc = Comps.descender_off(idx(space_idx));
    valid_idx = find((asc<=0 | asc < heights) & (desc<=0 | desc < heights));
    space_idx = space_idx(valid_idx);
    invalid_idx = setdiff(1:length(num_spaces), space_idx);
    num_spaces(invalid_idx) = 0;
end


tot_num_spaces = sum(num_spaces);
curr_pos = Comps.max_comp + 1;

Comps.space_width = space_width;
Comps.space_height = space_height;
Comps.max_comp = Comps.max_comp + tot_num_spaces;
space_comps = curr_pos:Comps.max_comp;
Comps.clust = [Comps.clust; Clust.num + 1 + zeros(tot_num_spaces,1,'uint16')];
Comps.pos = [Comps.pos; zeros(tot_num_spaces,4,'uint16')];
Comps.pg = [Comps.pg; zeros(tot_num_spaces,1,'uint32')];
Comps.nb = [Comps.nb; zeros(tot_num_spaces,4)];
Comps.nb_dist = [Comps.nb_dist; zeros(tot_num_spaces,4,'uint16')];
Comps.scale_factor = [Comps.scale_factor; ones(tot_num_spaces,1)];

if Comps.found_lines
    Comps.line = [Comps.line; zeros(tot_num_spaces,1, 'uint64')];
    Comps.descender_off = [Comps.descender_off;zeros(tot_num_spaces,1,'int16')];
    Comps.ascender_off = [Comps.ascender_off;zeros(tot_num_spaces,1,'int16')];
end
if Comps.found_true_labels
    Comps.truth_label = [Comps.truth_label; uint8(' ')+zeros(tot_num_spaces,1)];
end


for ii=space_idx'
    ns = num_spaces(ii);
    lnb = idx(ii);
    rows = curr_pos:curr_pos+double(ns)-1;
    old_nb = Comps.nb(lnb,3);
    old_dist = Comps.nb_dist(lnb,3);

    %first fixup neighbour info for the left neighbour
    Comps.nb(lnb,3) = curr_pos;
    Comps.nb(curr_pos,1) = lnb;
    Comps.nb_dist(lnb,3) = 1;
    Comps.nb_dist(curr_pos,1) = 1;

    %add the positions of the new components
    new_pos = Comps.pos(lnb,:);
    new_pos(1) = new_pos(3) + 1;
    new_pos(3) = new_pos(1) + space_width - 1;
    new_pos(4) = uint16(int16(new_pos(4)) - Comps.descender_off(lnb));
    new_pos(2) = new_pos(4) - space_height + 1;
    Comps.pos(rows,:) = repmat(new_pos, ns, 1);

    if ns > 1
        curr_nb = rows(2);
        while curr_nb <= rows(end)
            %fixup l-r positions and neighbours
            Comps.pos(curr_nb,1) = Comps.pos(curr_nb-1,3) + 1;
            Comps.pos(curr_nb,3) = Comps.pos(curr_nb,1) + space_width - 1;
            Comps.nb(curr_nb,1) = curr_nb - 1;
            Comps.nb_dist(curr_nb,1) = 1;
            Comps.nb(curr_nb-1,3) = curr_nb;
            old_dist = old_dist - space_width;
            curr_nb = curr_nb + 1;
        end
    end
    Comps.nb(rows(end),3) = old_nb;
    Comps.nb_dist(rows(end),3) = old_dist;
    Comps.nb(old_nb,1) = rows(end);
    Comps.nb_dist(old_nb,1) = old_dist;

    %other fields
    Comps.pg(rows) = Comps.pg(lnb);
    if Comps.found_lines
        Comps.line(rows) = Comps.line(lnb);
    end

    curr_pos = rows(end)+1;
end

Clust.num = Clust.num + 1;
Clust.num_comps(Clust.num) = tot_num_spaces;
Clust.mode_num(Clust.num) = Clust.num_comps(Clust.num);
Clust.comps{Clust.num} = space_comps;
Clust.avg{Clust.num} = bg_val + zeros(space_height, space_width);
Clust.density(Clust.num) = 0;
Clust.norm_sq(Clust.num) = 0;
Clust.refined(Clust.num) = 1;
Clust.changed(Clust.num) = 0;
if Clust.found_offsets
    Clust.descender_off(Clust.num) = 0;
    Clust.ascender_off(Clust.num) = 0;
    Clust.class(Clust.num) = assign_class({' '});
end
if isfield(Clust, 'truth_label')
    Clust.truth_label{Clust.num} = ' ';
end


Comps.model_spaces = true;
Clust.model_spaces = true;
Clust.space_dists = trans_dist;

fprintf('finished creating space character model\n');
