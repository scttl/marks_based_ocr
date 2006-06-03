function  display_neighbours(Clust, Comps, comp)
% display_neighbours    Draws components and immediate surrounding neighbours
%
%   display_neighbours(Clust, Comps, comp)
%
%   Clust should be an array of structs, each of which is assumed to contain 
%   a avg field which is a matrix giving the average pixel intensity
%   corresponding to the elements of that cluster.
%
%   Comps should be a 2 or 3 dimensional image matrix, where each entry
%   represents a pixel, and each 'on' pixel is labelled with the component 
%   number to which it belongs.
%
%   comp should be either an individual valid component number, or a vector of
%   components, each of which will be shown on screen for a configurable time
%   period before the next is shown.

% CVS INFO %
%%%%%%%%%%%%
% $Id: display_neighbours.m,v 1.1 2006-06-03 20:55:54 scottl Exp $
%
% REVISION HISTORY
% $Log: display_neighbours.m,v $
% Revision 1.1  2006-06-03 20:55:54  scottl
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
if nargin ~= 3
    error('incorrect number of arguments specified!');
end

while length(comp) > 0
    %determine the neighbours of this component
    [cl, off] = get_cl_off(Clust, comp(1));
    if isnan(cl)
        %this component isn't found, skip on to the next one
        comp = comp(2:end);
        continue;
    end
    pg = Clust(cl).pg(off);
    l = Clust(cl).pos(off,1); t = Clust(cl).pos(off,2);
    r = Clust(cl).pos(off,3); b = Clust(cl).pos(off,4);
    nb = Clust(cl).nb(off,:);
    if nb(1) ~= 0
        %left neighbour exists
        [lcl, loff] = get_cl_off(Clust, nb(1));
        l = Clust(lcl).pos(loff,1);
        if Clust(lcl).pos(loff,2) < t
            t = Clust(lcl).pos(loff,2);
        end
        if Clust(lcl).pos(loff,4) > b
            b = Clust(lcl).pos(loff,4);
        end
    end
    if nb(2) ~= 0
        %top neighbour exists
        [tcl, toff] = get_cl_off(Clust, nb(2));
        t = Clust(tcl).pos(toff,2);
        if Clust(tcl).pos(toff,1) < l
            l = Clust(tcl).pos(toff,1);
        end
        if Clust(tcl).pos(toff,3) > r
            r = Clust(tcl).pos(toff,3);
        end
    end
    if nb(3) ~= 0
        %right neighbour exists
        [rcl, roff] = get_cl_off(Clust, nb(3));
        r = Clust(rcl).pos(roff,3);
        if Clust(rcl).pos(roff,2) < t
            t = Clust(rcl).pos(roff,2);
        end
        if Clust(rcl).pos(roff,4) > b
            b = Clust(rcl).pos(roff,4);
        end
    end
    if nb(4) ~= 0
        %bottom neighbour exists
        [bcl, boff] = get_cl_off(Clust, nb(4));
        b = Clust(bcl).pos(boff,4);
        if Clust(bcl).pos(boff,1) < l
            l = Clust(bcl).pos(boff,1);
        end
        if Clust(bcl).pos(boff,3) > r
            r = Clust(bcl).pos(boff,3);
        end
    end
    M = label2rgb(Comps(t:b,l:r, pg), 'white', 'k');
    x = Clust(cl).pos(off,:);
    titlestr = {sprintf('Component %d: [%d,%d,%d,%d]', comp(1),x)};
    M(x(2)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(comp_col,1,x(3)-x(1)+1);
    M(x(2)-t+1:x(4)-t+1,x(1)-l+1,:) = repmat(comp_col,x(4)-x(2)+1,1);
    M(x(2)-t+1:x(4)-t+1,x(3)-l+1,:) = repmat(comp_col,x(4)-x(2)+1,1);
    M(x(4)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(comp_col,1,x(3)-x(1)+1);
    if nb(1) ~= 0
        x = Clust(lcl).pos(loff,:);
        M(x(2)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(nb_col,1,x(3)-x(1)+1);
        M(x(2)-t+1:x(4)-t+1,x(1)-l+1,:) = repmat(nb_col,x(4)-x(2)+1,1);
        M(x(2)-t+1:x(4)-t+1,x(3)-l+1,:) = repmat(nb_col,x(4)-x(2)+1,1);
        M(x(4)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(nb_col,1,x(3)-x(1)+1);
        titlestr = {titlestr{:}, sprintf('lnb [%d,%d,%d,%d]', x)};
    end
    if nb(2) ~= 0
        x = Clust(tcl).pos(toff,:);
        M(x(2)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(nb_col,1,x(3)-x(1)+1);
        M(x(2)-t+1:x(4)-t+1,x(1)-l+1,:) = repmat(nb_col,x(4)-x(2)+1,1);
        M(x(2)-t+1:x(4)-t+1,x(3)-l+1,:) = repmat(nb_col,x(4)-x(2)+1,1);
        M(x(4)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(nb_col,1,x(3)-x(1)+1);
        titlestr = {titlestr{:}, sprintf('tnb [%d,%d,%d,%d]', x)};
    end
    if nb(3) ~= 0
        x = Clust(rcl).pos(roff,:);
        M(x(2)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(nb_col,1,x(3)-x(1)+1);
        M(x(2)-t+1:x(4)-t+1,x(1)-l+1,:) = repmat(nb_col,x(4)-x(2)+1,1);
        M(x(2)-t+1:x(4)-t+1,x(3)-l+1,:) = repmat(nb_col,x(4)-x(2)+1,1);
        M(x(4)-t+1,x(1)-l+1:x(3)-l+1,:) = repmat(nb_col,1,x(3)-x(1)+1);
        titlestr = {titlestr{:}, sprintf('rnb [%d,%d,%d,%d]', x)};
    end
    if nb(4) ~= 0
        x = Clust(bcl).pos(boff,:);
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
