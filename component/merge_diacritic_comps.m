function Comps = merge_diacritic_comps(Comps, varargin)
% MERGE_DIACRITIC_COMPS  Merge where one comp is a diacritical mark of the other
%
% COMPS = merge_diacritic_comps(COMPS)
%
% This function takes a component structure (like that returned from get_comps)
% after it has been run through get_lines() to determine line information.  It
% then uses this line information to determine when two vertically aligned 
% components should be merged into a single component like the dot above and 
% i, or accents above characters (like the circonflexe '^' above an 'e' in 
% French)
%
% The updated component struct is returned
%
% The criteria for determining when two vertically aligned components are
% merged is simple:
%   * both belong to the same line
%   * they are within some vertical threshold of pixels (set below) or one is
%     vertically contained within the other.
%   * if one component has left and right boundaries l1 and r1, and the other
%     l2 and r2, then either l1 <= l2 and r1 >= r2, or l1 >= l2 and r1 <= r2
%     (ie. one's width is contained in the other)
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: merge_diacritic_comps.m,v 1.1 2006-12-17 19:52:47 scottl Exp $
%
% REVISION HISTORY
% $Log: merge_diacritic_comps.m,v $
% Revision 1.1  2006-12-17 19:52:47  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%set this to true to display characters that are merged on screen.
display_matches = false;

%how many pixels apart are they allowed to be, to be considered a match?
max_dist = 7;


% CODE START %
%%%%%%%%%%%%%%
if nargin < 1
    error('incorrect number of arguments specified!');
elseif nargin > 1
    process_optional_args(varargin{:})
end

if ~Comps.found_lines
    error('Comps must contain line information');
end

ii = 1;
while ii<=Comps.max_comp
    %line matches
    idx = find(Comps.line(ii) == Comps.line);
    idx = idx(idx ~= ii);
    %position matches
    m_idx = (Comps.pos(ii,1) <= Comps.pos(idx,1) & ...
                  Comps.pos(ii,3) >= Comps.pos(idx,3)) | ...
                 (Comps.pos(ii,1) >= Comps.pos(idx,1) & ...
                  Comps.pos(ii,3) <= Comps.pos(idx,3));
    idx = idx(m_idx);
    %distance matches
    m_idx = abs(single(Comps.pos(ii,2))-single(Comps.pos(idx,4)))<=max_dist ...
          | abs(single(Comps.pos(idx,2))-single(Comps.pos(ii,4)))<=max_dist ...
          | (Comps.pos(ii,2)<=Comps.pos(idx,2) & Comps.pos(ii,4) >= ...
            Comps.pos(idx,4)) ...
          | (Comps.pos(ii,2)>=Comps.pos(idx,2) & Comps.pos(ii,4) <= ...
            Comps.pos(idx,4));
    idx = idx(m_idx);

    if ~isempty(idx)
        %match found
        idx_in_ii = false;
        ii_in_idx = false;
        if length(idx) > 1
            warning('MBOCR:MultDiac', 'found more than one diacritic match');
            idx = idx(1);
        end
        if Comps.pos(ii,2) <= Comps.pos(idx,2)
            top = ii;
            if Comps.pos(ii,4) >= Comps.pos(idx,4)
                %idx contained wholly within ii
                bot = ii;
                idx_in_ii = true;
            else
                bot = idx;
            end
        else
            top = idx;
            if Comps.pos(ii,4) <= Comps.pos(idx,4)
                %ii wholly contained within idx
                ii_in_idx = true;
                bot = idx;
            else
                bot = ii;
            end
        end

        if display_matches
            fprintf('match between %d and %d\n', top, bot);
            imgs = get_comp_imgs(Comps, [top, bot]);
            subplot(2,2,1), imshow(imgs{1}), axis equal;
            subplot(2,2,3), imshow(imgs{2}), axis equal;
        end

        %update component information fields
        Comps.pos(ii,1) = min(Comps.pos(top,1), Comps.pos(bot,1));
        Comps.pos(ii,2) = Comps.pos(top,2);
        Comps.pos(ii,3) = max(Comps.pos(top,3), Comps.pos(bot,3));
        Comps.pos(ii,4) = Comps.pos(bot,4);

        Comps.nb(ii,2) = Comps.nb(top,2);
        Comps.nb_dist(ii,2) = Comps.nb_dist(top,2);
        Comps.nb(ii,4) = Comps.nb(bot,4);
        Comps.nb_dist(ii,4) = Comps.nb_dist(bot,4);
        if ii_in_idx || ~idx_in_ii && Comps.nb_dist(ii,1) > Comps.nb_dist(idx,1)
            Comps.nb(ii,1) = Comps.nb(idx,1);
            Comps.nb_dist(ii,1) = Comps.nb_dist(idx,1);
        end
        if ii_in_idx || ~idx_in_ii && Comps.nb_dist(ii,3) > Comps.nb_dist(idx,3)
            Comps.nb(ii,3) = Comps.nb(idx,3);
            Comps.nb_dist(ii,3) = Comps.nb_dist(idx,3);
        end
        nb_idx = find(Comps.nb(:,1) == ii | Comps.nb(:,1) == idx);
        Comps.nb(nb_idx,1) = ii;
        Comps.nb_dist(nb_idx,1) = Comps.pos(nb_idx,1) - Comps.pos(ii,3);
        nb_idx = find(Comps.nb(:,2) == ii | Comps.nb(:,2) == idx);
        Comps.nb(nb_idx,2) = ii;
        Comps.nb_dist(nb_idx,2) = Comps.pos(nb_idx,2) - Comps.pos(ii,4);
        nb_idx = find(Comps.nb(:,3) == ii | Comps.nb(:,3) == idx);
        Comps.nb(nb_idx,3) = ii;
        Comps.nb_dist(nb_idx,3) = Comps.pos(ii,1) - Comps.pos(nb_idx,3);
        nb_idx = find(Comps.nb(:,4) == ii | Comps.nb(:,4) == idx);
        Comps.nb(nb_idx,4) = ii;
        Comps.nb_dist(nb_idx,4) = Comps.pos(ii,2) - Comps.pos(nb_idx,4);

        if ~isempty(Comps.descender_off) && bot ~= ii
            Comps.descender_off(ii) = Comps.descender_off(bot);
        end
        if ~isempty(Comps.ascender_off) && top ~= ii
            Comps.ascender_off(ii) = Comps.ascender_off(top);
        end
        if ~isempty(Comps.scale_factor)
            Comps.scale_factor(ii) = Comps.modal_height / ...
                              double(Comps.pos(ii,4) - Comps.pos(ii,2) + 1);
        end

        %@@ what should we do about cluster information??.  For now just leave
        %as ii's original cluster

        if display_matches
            img = get_comp_imgs(Comps, ii);
            subplot(2,2,2), imshow(img{1}), axis equal;
            pause(.1);
        end

        %remove the merged component information
        keep_idx = [1:idx-1, idx+1:Comps.max_comp];
        Comps.max_comp = Comps.max_comp - 1;
        Comps.pos = Comps.pos(keep_idx,:);
        Comps.pg = Comps.pg(keep_idx);
        Comps.nb = Comps.nb(keep_idx,:);
        lrg_idx = Comps.nb > idx;
        Comps.nb(lrg_idx) = Comps.nb(lrg_idx) - 1;
        Comps.nb_dist = Comps.nb_dist(keep_idx,:);
        Comps.line = Comps.line(keep_idx);
        if ~isempty(Comps.descender_off)
            Comps.descender_off = Comps.descender_off(keep_idx);
        end
        if ~isempty(Comps.ascender_off)
            Comps.ascender_off = Comps.ascender_off(keep_idx);
        end
        if ~isempty(Comps.scale_factor)
            Comps.scale_factor = Comps.scale_factor(keep_idx);
        end
        if ~isempty(Comps.clust)
            Comps.clust = Comps.clust(keep_idx);
        end


        %ensure that the index moves to the next valid component
        if idx <= ii
            ii = ii - 1;
        end
    end

    ii = ii + 1;
end

