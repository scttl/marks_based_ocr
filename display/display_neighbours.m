function  display_neighbours(Comps, comp)
% display_neighbours    Draws components and immediate surrounding neighbours
%
%   display_neighbours(Comps, comp)
%
%   Comps should be a struct containing various fields (see cluster_comps for
%   specifics)
%
%   comp should be either an individual valid component number, or a vector of
%   components, each of which will be shown on screen for a configurable time
%   period before the next is shown.

% CVS INFO %
%%%%%%%%%%%%
% $Id: display_neighbours.m,v 1.4 2006-10-09 16:31:40 scottl Exp $
%
% REVISION HISTORY
% $Log: display_neighbours.m,v $
% Revision 1.4  2006-10-09 16:31:40  scottl
% small bugfix to ensure dimensions selected are correct.
%
% Revision 1.3  2006/07/05 01:06:49  scottl
% rewritten based on new cluster and component structures.
%
% Revision 1.2  2006/06/19 20:59:04  scottl
% changes to reflect that Comps is now a cell array.
%
% Revision 1.1  2006/06/03 20:55:54  scottl
% Initial check-in.
%
%

% LOCAL VARS %
%%%%%%%%%%%%%%
display_time = 2;  %display time in seconds, if multiple comps passed
comp_col = reshape([0,0,255],1,1,3);  %this is blue
nb_col = reshape([255,0,0],1,1,3);  %this is red

% CODE START %
%%%%%%%%%%%%%%
if nargin ~= 2
    error('incorrect number of arguments specified!');
end

while length(comp) > 0
    %determine the neighbours of this component
    if comp(1) < 1 || comp(1) > Comps.max_comp
        %this component isn't found, skip on to the next one
        comp = comp(2:end);
        continue;
    end

    M = imread(Comps.files{Comps.pg(comp(1))});
    x = Comps.pos(comp(1),:);
    nb = Comps.nb(comp(1),:);
    l = x(1); t = x(2); r = x(3); b = x(4);
    if nb(1) ~= 0
        %left neighbour exists
        l = min(l, Comps.pos(nb(1),1));
        t = min(t, Comps.pos(nb(1),2));
        b = max(b, Comps.pos(nb(1),4));
    end
    if nb(2) ~= 0
        %top neighbour exists
        t = min(t, Comps.pos(nb(2),2));
        l = min(l, Comps.pos(nb(2),1));
        r = max(r, Comps.pos(nb(2),3));
    end
    if nb(3) ~= 0
        %right neighbour exists
        r = max(r, Comps.pos(nb(3),3));
        t = min(t, Comps.pos(nb(3),2));
        b = max(b, Comps.pos(nb(3),4));
    end
    if nb(4) ~= 0
        %bottom neighbour exists
        b = max(b, Comps.pos(nb(4),4));
        l = min(l, Comps.pos(nb(4),1));
        r = max(r, Comps.pos(nb(4),3));
    end
    M = label2rgb(~M(t:b,l:r), 'white', 'k');
    titlestr = {sprintf('Component %d: [%d,%d,%d,%d]', comp(1),x)};

    M(x(2)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(comp_col,1,x(3)-x(1)+1);
    M(x(2)-t+1:x(4)-t+1,x(1)-l+1,:) = repmat(comp_col,x(4)-x(2)+1,1);
    M(x(2)-t+1:x(4)-t+1,x(3)-l+1,:) = repmat(comp_col,x(4)-x(2)+1,1);
    M(x(4)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(comp_col,1,x(3)-x(1)+1);
    if nb(1) ~= 0
        x = Comps.pos(nb(1),:);
        M(x(2)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(nb_col,1,x(3)-x(1)+1);
        M(x(2)-t+1:x(4)-t+1,x(1)-l+1,:) = repmat(nb_col,x(4)-x(2)+1,1);
        M(x(2)-t+1:x(4)-t+1,x(3)-l+1,:) = repmat(nb_col,x(4)-x(2)+1,1);
        M(x(4)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(nb_col,1,x(3)-x(1)+1);
        titlestr = {titlestr{:}, sprintf('lnb [%d,%d,%d,%d]', x)};
    end
    if nb(2) ~= 0
        x = Comps.pos(nb(2),:);
        M(x(2)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(nb_col,1,x(3)-x(1)+1);
        M(x(2)-t+1:x(4)-t+1,x(1)-l+1,:) = repmat(nb_col,x(4)-x(2)+1,1);
        M(x(2)-t+1:x(4)-t+1,x(3)-l+1,:) = repmat(nb_col,x(4)-x(2)+1,1);
        M(x(4)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(nb_col,1,x(3)-x(1)+1);
        titlestr = {titlestr{:}, sprintf('tnb [%d,%d,%d,%d]', x)};
    end
    if nb(3) ~= 0
        x = Comps.pos(nb(3),:);
        M(x(2)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(nb_col,1,x(3)-x(1)+1);
        M(x(2)-t+1:x(4)-t+1,x(1)-l+1,:) = repmat(nb_col,x(4)-x(2)+1,1);
        M(x(2)-t+1:x(4)-t+1,x(3)-l+1,:) = repmat(nb_col,x(4)-x(2)+1,1);
        M(x(4)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(nb_col,1,x(3)-x(1)+1);
        titlestr = {titlestr{:}, sprintf('rnb [%d,%d,%d,%d]', x)};
    end
    if nb(4) ~= 0
        x = Comps.pos(nb(4),:);
        M(x(2)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(nb_col,1,x(3)-x(1)+1);
        M(x(2)-t+1:x(4)-t+1,x(1)-l+1,:) = repmat(nb_col,x(4)-x(2)+1,1);
        M(x(2)-t+1:x(4)-t+1,x(3)-l+1,:) = repmat(nb_col,x(4)-x(2)+1,1);
        M(x(4)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(nb_col,1,x(3)-x(1)+1);
        titlestr = {titlestr{:}, sprintf('bnb [%d,%d,%d,%d]', x)};
    end
    imshow(M);
    drawnow;
    xlabel(titlestr');
    if length(comp) > 1
        pause(display_time);
    end
    comp = comp(2:end);
end
