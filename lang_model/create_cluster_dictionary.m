function [Clust, Comps] = create_cluster_dictionary(Clust, Comps, varargin)
% CREATE_CLUSTER_DICTIONARY Create a list and counts of blobs (chars)
%
%   [Clust, Comps] = CREATE_CLUSTER_DICTIONARY(Clust, Comps, [var1, val1]...)
%
%   This function uses already clustered page data to come up with counts of 
%   cluster transitions for use in constructing a dictionary.
%
%   We assume that the Comps struct passed contains line related information
%   like ascender and descender offsets (see get_lines() for more info).
%   Similarily, the Clust struct should also contain this information
%
%   Any of the variables below can be overridden by passing in its name and new
%   value as an additional pair of arguments
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: create_cluster_dictionary.m,v 1.11 2007-01-13 18:16:55 scottl Exp $
%
% REVISION HISTORY
% $Log: create_cluster_dictionary.m,v $
% Revision 1.11  2007-01-13 18:16:55  scottl
% just check if the line field exists, instead of checking for related fields.
%
% Revision 1.10  2007-01-08 22:08:55  scottl
% fixed potential bug in repeated neighbour following.  Normalized
% position counts.
%
% Revision 1.9  2006-12-19 21:40:46  scottl
% assume that line offset information is already stored,
% rewrite how space models get added.  Added cluster positional count
% information.
%
% Revision 1.8  2006-09-22 18:01:41  scottl
% added MSGID to warning message
%
% Revision 1.7  2006/08/24 21:40:07  scottl
% added ability to use the mode instead of taking the average of cluster
% intensities while refining.
%
% Revision 1.6  2006/08/14 01:26:42  scottl
% removed Dummy variable
%
% Revision 1.5  2006/08/08 03:18:11  scottl
% fixed a couple of bugs with number of parameters, and space refined count.
%
% Revision 1.4  2006/07/05 01:00:04  scottl
% re-written after changing cluster and component structures.  Character counts,
% and bigram model is now stored in Clust.
%
% Revision 1.3  2006/06/21 21:47:12  scottl
% use if test and iterate over all items instead of using i_offs
%
% Revision 1.2  2006/06/19 21:37:57  scottl
% implemented baseline offset calculation, addition of a 'space' cluster.
%
% Revision 1.1  2006/06/12 20:57:50  scottl
% initial revision
%


% LOCAL VARS %
%%%%%%%%%%%%%%
model_spaces = true;  %should we add a space character bitmap and bigram count?

smoothing_counts = 1;  %add plus-one smoothing to ensure all transitions are
                       %possible.

%up to what length non-space sequence of blobs (i.e. word) should we include in
%our positional counts
max_word_len = 10;


% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 2
    error('incorrect number of arguments passed');
elseif nargin > 2
    process_optional_args(varargin{:});
end

%ensure line information present in clusters and components
if ~isfield(Comps, 'line') || isempty(Comps.line)
    error('we require components to contain line information');
end

num_pgs = size(Comps.pg_size,1);

%add the 'space character' Cluster if specified and it doesn't already exist
if model_spaces && ~Clust.model_spaces
    [Clust, Comps] = add_space_model(Clust, Comps);
    fprintf('%.2fs: finished adding space model\n', toc);
end

%sort the clusters to ensure they are ordered by number of elements
[Clust, Comps] = sort_clusters(Clust, Comps);

if model_spaces
    %attempt to infer which Cluster belongs to the space character.  In any 
    %reasonable document this will almost certainly be the first cluster.  Use 
    %our assigned truth value as a sanity check
    space_idx = 1;
    if isfield(Clust ,'truth_label') && ~isempty(Clust.truth_label)
        if ~strcmp(Clust.truth_label{1}, ' ')
            %attempt to find it using true labels
            space_idx = find(strcmp(Clust.truth_label, ' '));
            if isempty(space_idx)
                warning('MBOCR:spaceNotFound', ...
                'Space character not found.  Will be ignored in counts.\n');
                model_spaces = false;
            end
        end
    else
        warning('MBOCR:spaceNotFound', ...
        'Space character nor truth labels found.  Spaces ignored in counts\n');
        model_spaces = false;
    end
end

%now create the bigram counts
idx = find(Comps.nb(:,3) ~= 0);
Trans = double([Comps.clust(idx), Comps.clust(Comps.nb(idx,3))]);
Clust.num_trans = size(Trans,1);
Clust.bigram = zeros(Clust.num);
for ii=1:size(Trans,1)
    fr = Trans(ii,1); to = Trans(ii,2);
    Clust.bigram(fr,to) = Clust.bigram(fr,to) + 1;
end

%add smothing counts to the totals;
Clust.bigram = Clust.bigram + smoothing_counts;
Z = sum(Clust.bigram,2);
if any(Z == 0)
    %this can happen if a component never transitions to another component
    %i.e. only appeaars on the far right edge of pages.  Smoothing should
    %normally take care of this, but to prevent dividing by 0, we augment
    %the sum
    zero_rows = find(Z == 0);
    warning('MBOCR:noTrans', 'No transitions seen from cluster %d\n',zero_rows);
    Z(zero_rows) = 1;
end
Clust.bigram = Clust.bigram ./ repmat(Z,1,Clust.num);
fprintf('%.2fs: finished creating character bigram matrix\n', toc);

%now determine component positional counts
if ~model_spaces
    warning('MBOCR:spaceNotFound', ...
            'no positional counts can be taken.  Spaces not found');
end
Clust.pos_count = cell(1,max_word_len);
Clust.pos_total = Comps.max_comp;
for ii=unique(Comps.line)'
    m = find(Comps.line == ii);
    [idx,idx] = min(Comps.pos(m,1));  %find left-most component on this line
    nbs = m(idx);
    while Comps.nb(nbs(end),3) ~= 0
        if any(nbs == Comps.nb(nbs(end),3))
            %repeated neighbours, skip to the next neighbour listing this one
            next = find(Comps.nb(:,1) == Comps.nb(nbs(end),3));
            while ~isempty(next) && any(nbs == next(1))
                next = next(2:end);
            end
            if isempty(next)
                break;
            else
                nbs = [nbs; next(1)];
            end
        else
            nbs = [nbs; Comps.nb(nbs(end),3)];
        end
    end
    spc_pos = find(Comps.clust(nbs) == space_idx);
    pos = 0;
    while ~isempty(spc_pos)
        len = spc_pos(1)-1 - pos;
        if len > 0 && len <= max_word_len
            Clust.pos_count{len} = [Clust.pos_count{len}; ...
                                    Comps.clust(nbs(pos+1:spc_pos(1)-1))'];
        end
        pos = spc_pos(1);
        spc_pos = spc_pos(2:end);
    end
    %must also include word between the last space and the end of the line
    len = length(nbs) - pos;
    if len > 0 && len <= max_word_len
        Clust.pos_count{len} = [Clust.pos_count{len}; ...
                                Comps.clust(nbs(pos+1:end))'];
    end
end
fprintf('%.2fs: finished creating cluster word lists for pos counts\n', toc);
%now go through the cluster "word" lists created above, and update the counts
for ii=1:max_word_len
    w = Clust.pos_count{ii};
    Clust.pos_count{ii} = zeros(Clust.num,ii);
    for jj=1:Clust.num
        Clust.pos_count{ii}(jj,:) = sum(w == jj) ./ Clust.pos_total;
    end
end
fprintf('%.2fs: finished counting positions of components\n', toc);

