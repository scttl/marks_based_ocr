function clust = word_subimage_cluster(Lines, varargin)
% WORD_SUBIMAGE_CLUSTER   Segment lines into components using subimages
%
%   [Clust, Comps] = WORD_SUBIMAGE_CLUSTER(LINES, [VAR1, VAL1]...)
%   LINES is a struct like that returned from get_lines()
%
%   VAR1 and VAL1 are optional, and can be used to override the default values
%   for the LOCAL VARS defined below.  Each VAR1 should be a string giving the
%   exact name of the variable to override.  Each VAL1 should be the updated
%   value to use for that variable.
%
%   This algorithm is taken from "Character segmentation using visual
%   inter-word constraints in text page" paper by Hong and Hull (SPIE'95)
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: word_subimage_cluster.m,v 1.2 2006-12-19 21:43:25 scottl Exp $
%
% REVISION HISTORY
% $Log: word_subimage_cluster.m,v $
% Revision 1.2  2006-12-19 21:43:25  scottl
% check in of partially implemented word subimage clustering.  Since
% this technique does not perform particularly well, I did not bother
% with the rest of the implementation.
%
% Revision 1.1  2006-12-05 15:58:28  scottl
% initial incomplete revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%this controls when a word is split into multiple parts (should be a percentage
%of bounding-box area that can mismatch
split_thresh = 0.01;

%this determines when identical words are initially grouped together 
euc_match_thresh = 0.005;
haus_match_thresh = 1.0;

%prevent showing warning each time we call match_refine etc.
override_display = 'OFF';   %other legal value is 'ON'

%this determines how close to an edge a match must be to still be considerd a
%match at an edge (value is specified in pixels)
edge_tolerance = 3;

%how many pixels must be on in an image for it to be considered valid (and not
%noise)
min_area = 20;

%how many pixels wide must an image be for it to be considered valid
min_width = 7;

%should we display the matches that we find to the screen?
display_matches = false;

%these variables can be given queue img and cluster imgs to use (instead of
%estimating them from the Lines struct
passed_queue = [];
passed_clust = [];

% CODE START %
%%%%%%%%%%%%%%
tic;

%process input arguments
if nargin < 1
    error('Comps struct must be specified as the first argument!');
elseif nargin > 1
    process_optional_args(varargin{:});
end
prev_comp_warn = warning(override_display, 'MBOCR:noCompImg');
prev_blank_warn = warning(override_display, 'MBOCR:blankImg');
prev_over_warn = warning(override_display, 'MBOCR:override');

if ~isempty(passed_queue) && ~isempty(passed_clust)
    fprintf('using passed queue and cluster images\n');
    queue.img = passed_queue;
    clust.img = passed_clust;
else
    %determine space width threshold @@TO DO@@
    space_width = 9;
    
    queue.sum = [];
    queue.img = cell(0);
    clust.sum = [];
    clust.img = cell(0);
    
    %break the lines into words
    curr_pg = Lines.pg(1);
    curr_img = ~imread(Lines.files{curr_pg});
    for ii=1:Lines.num
        if Lines.pg(ii) ~= curr_pg
            curr_pg = Lines.pg(ii);
            curr_img = ~imread(Lines.files(curr_pg));
        end
        line_img = curr_img(Lines.pos(ii,2):Lines.pos(ii,4), ...
                            Lines.pos(ii,1):Lines.pos(ii,3));
    
        queue.img = [queue.img; get_words(line_img, space_width)];
    end
    
    %match identical words (average them) @@TO DO@@
    queue.norm = zeros(length(queue.img),1);
    for ii=1:length(queue.img)
        queue.norm(ii) = sum(queue.img{ii}(:));
    end
    idx = 1;
    fprintf('queue length: %d\n', length(queue.img));
    while (idx <= length(queue.img))
        fprintf('item: %d ', idx);
        curr_img = queue.img{idx};
        curr_norm = queue.norm(idx);
        d = euc_dist(curr_img, queue.img, curr_norm, queue.norm);
        %d = ham_dist(curr_img, queue.img);
        match = find(d <= euc_match_thresh);
        match = setdiff(match, idx);  %index will always match with 0 distance
        if ~isempty(match)
            if display_matches
                subplot(1,2,1), imshow(curr_img);
                subplot(1,2,2), imshow(imgcell2mat(queue.img(match)));
                d(match)
                pause(1);
            end
            queue.img = queue.img(setdiff(1:length(queue.img), match));
            queue.norm = queue.norm(setdiff(1:length(queue.norm), match));
        end
        fprintf(' %d matches\r', length(match));
        idx = idx+1;
    end
end
fprintf('\nqueue length: %d\n', length(queue.img));


%iteratively split items in the queue, adding them to the cluster
while ~isempty(queue.img)
    %subplot(1,2,1), imshow(imgcell2mat(queue.img));
    %subplot(1,2,2), imshow(imgcell2mat(clust.img));
    %pause(1);
    fprintf('queue length %d, cluster length %d\n', length(queue.img), ...
            length(clust.img));
    curr_img = queue.img{1};
    queue.img = queue.img(2:end);
    if sum(curr_img(:)) < min_area || size(curr_img,2) < min_width;
        continue;
    end
    %attempt to match with a cluster image
    d = hausdorff_dist(curr_img, clust.img);
    [d,idx] = sort(d);
    if ~isempty(d) && d(1) <= haus_match_thresh
        %add this image to the cluster @@TO DO@@
        fprintf('image matches cluster %d\n', idx(1));
        if display_matches
            subplot(1,2,1), imshow(curr_img);
            subplot(1,2,2), imshow(clust.img{idx(1)});
            pause(1);
        end
    else
        cus = size(curr_img);
        rem_idx = [];
        add_cluster = true;
        for ii=1:length(clust.img)
            cls = size(clust.img{ii});
            if cus(2) <= cls(2)
                %check if this image is a subimage the current cluster image
                sub_img = curr_img;
                sub_sz = cus;
                img = clust.img{ii};
                img_sz = cls;
            else
                %check if current cluster image is a subimage of this image
                sub_img = clust.img{ii};
                sub_sz = cls;
                img = curr_img;
                img_sz = cus;
            end
            [t,l,score] = subimage_match(sub_img, img, 'thresh', split_thresh);
            if ~isempty(l)
                if display_matches
                    display_subimage_match(sub_img, img, t, l);
                    pause(1);
                end
                r = l+sub_sz(2)-1;
                if l(1) <= edge_tolerance
                    l = l(2:end);
                else
                    %include the piece before the first match
                    r = [0; r];
                end
                if ~isempty(l) && l(end)+sub_sz(2)-1 >= img_sz(2)-edge_tolerance
                    r = r(1:end-1);
                else
                    %include the piece after the last match
                    l = [l; img_sz(2)+1];
                end
                fprintf('found %d matches %d\n', length(r), length(l));
                new_items = cell(length(r),1);
                for jj=1:length(r)
                    new_items{jj} = crop_image(img(:,r(jj)+1:l(jj)-1));
                end
                queue.img = [queue.img; new_items];
                if cus(2) <= cls(2)
                    fprintf('this image a subimage of cluster %d\n', ii);
                    %mark this cluster for removal
                    rem_idx = [rem_idx, ii];
                else
                    fprintf('cluster %d a subimage of this image\n', ii);
                    add_cluster = false;
                end
            end
        end

%                while ~isempty(l)
%                    %split the cluster image and add the pieces to the queue
%                    if l(1) <= edge_tolerance
%                        %left-edge match - 2pieces, but left is in cluster
%                        queue.img = [queue.img; ...
%                                     {crop_image(cl_img(:,l(1)+cus(2):end))}];
%                        fprintf('left edge match for image in cluster %d\n',ii);
%                    elseif l(1)+cus(2)-1 >= cls(2) - edge_tolerance
%                        %right-edge match - 2pieces, but right is in cluster
%                        queue.img = [queue.img; ...
%                                     {crop_image(cl_img(:,1:l(1)-1))}];
%                        fprintf('rght edge match for image in cluster %d\n',ii);
%                    else
%                        %middle match - 3 pieces (but 2nd is current image
%                        %so it isn't re-added)
%                        queue.img = [queue.img; ...
%                                     {crop_image(cl_img(:,1:l(1)-1))}; ...
%                                     {crop_image(cl_img(:,l(1)+cus(2):end))}];
%                        fprintf('central match for image in cluster %d\n',ii);
%                    end
%                    if length(l) == 1
%                        %mark this cluster for removal
%                        rem_idx = [rem_idx, ii];
%                    end
%                    l = l(2:end);
%                end
%            else
%                %check if the current cluster image is a subimage of this image
%                [t,l,score] = subimage_match(cl_img, curr_img, 'thresh', ...
%                              split_thresh);
%                if ~isempty(l)
%                    %split the current image and add the pieces to the queue
%                    if l(1) <= edge_tolerance
%                        %left-edge match - 2pieces
%                        queue.img = [queue.img; ...
%                                     {crop_image(curr_img(:,l(1)+cls(2):end))}];
%                        fprintf('left edge match for cluster %d in image\n',ii);
%                    elseif l(1)+cls(2)-1 >= cus(2) - edge_tolerance
%                        %right-edge match - 2pieces (but 2nd is current cluster
%                        queue.img = [queue.img; ...
%                                     {crop_image(curr_img(:,1:l(1)-1))}];
%                        fprintf('rght edge match for cluster %d in image\n',ii);
%                    else
%                        %middle match - 3 pieces (but 2nd is current cluster
%                        %so it isn't re-added)
%                        queue.img = [queue.img; ...
%                                     {crop_image(curr_img(:,1:l(1)-1))}; ...
%                                     {crop_image(curr_img(:,l(1)+cls(2):end))}];
%                        fprintf('central match for cluster %d in image\n',ii);
%                    end
%                    add_cluster = false;
%                    l = l(2:end);
%                end
%            end
%        end
        %remove split clusters
        if ~isempty(rem_idx)
            clust.img = clust.img(setdiff(1:length(clust.img), rem_idx));
        end
        if add_cluster
            %create a new cluster for the current image
            clust.img = [clust.img; {crop_image(curr_img)}];
        end
    end
end

if strcmp(override_display, 'OFF')
    warning(prev_comp_warn.state, 'MBOCR:noCompImg');
    warning(prev_blank_warn.state, 'MBOCR:blankImg');
    warning(prev_over_warn.state, 'MBOCR:override');
end

fprintf('%.2f: TOTAL NUMBER OF CLUSTERS FOUND: %d\n', toc, length(clust.img));



% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
