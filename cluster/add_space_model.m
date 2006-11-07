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
% $Id: add_space_model.m,v 1.2 2006-11-07 02:51:05 scottl Exp $
%
% REVISION HISTORY
% $Log: add_space_model.m,v $
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

%when attempting to infer the width of a space character, what is the maximum
%possible number of pixels wide we should check.
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
    
if isempty(space_width)
    %calculate space_width manually.  First peak should be interchar spacing but
    %second should be interword (approximately)
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
    fprintf('%.2fs: estimated space width at %d pixels\n', toc, space_width);
end

if isempty(space_height)
    %set the height equal to the width (by default)
    space_height = space_width;
end

%we will create spaces in the large gaps between neighbouring components that
%stretch over the baseline (assuming line information has been found).
num_spaces = floor(trans_dist / space_width);
space_idx = find(num_spaces >= 1);
if Comps.found_lines
    bot_pos = Comps.pos(idx(space_idx),4);
    %top_pos = Comps.pos(idx(space_idx),2);
    bl_pos = uint16(int16(bot_pos) - Comps.descender_off(idx(space_idx)));
    valid_idx = find(bl_pos <= bot_pos); %& (bl_pos - top_pos) >= space_height);
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
Clust.norm_sq(Clust.num) = 0;
Clust.refined(Clust.num) = 1;
Clust.changed(Clust.num) = 0;
if Clust.found_offsets
    Clust.descender_off(Clust.num) = 0;
    Clust.ascender_off(Clust.num) = 0;
end
if isfield(Clust, 'truth_label')
    Clust.truth_label(Clust.num) = ' ';
end


Comps.model_spaces = true;
Clust.model_spaces = true;

%update the bigram too?
%Clust.num_trans = size(Trans,1);
%Clust.bigram = zeros(Clust.num);
%for ii=1:size(Trans,1)
%    fr = Trans(ii,1); to = Trans(ii,2);
%    Clust.bigram(fr,to) = Clust.bigram(fr,to) + 1;
%end

fprintf('%.2fs: finished creating space character model\n', toc);
