function [Clust, chg_list] = scale_refine(Clust, rl, st)
% SCALE_REFINE  Attempt to match scaled versions of the same cluster
%
%   [Clust, chg_list] = scale_refine(Clust, refine_list, distance_thresh)
%
%   Clust should be an array of structs, each of which is assumed to contain 
%   several fields of a particular format.  See cluster_comps for details.
%
%   refine_list is optional and if specified, determines which items are
%   to be refined (ie searched for a match).  It should be a vector of cluster
%   indices and will be processed in the order given.  If not specified, this 
%   method will attempt to refine all clusters, starting from the highest 
%   numbered one.
%
%   distance_thresh is optional and if specified, determines the tolerance
%
%   The refined list of clusters is returned in Clust, as is the indices of 
%   those clusters that changed in chg_list

% CVS INFO %
%%%%%%%%%%%%
% $Id: scale_refine.m,v 1.2 2006-06-12 20:56:01 scottl Exp $
%
% REVISION HISTORY
% $Log: scale_refine.m,v $
% Revision 1.2  2006-06-12 20:56:01  scottl
% changed comps to a cell array (indexed by page) to allow different sized pages
%
% Revision 1.1  2006/06/03 20:55:48  scottl
% Initial check-in.
%
%

% LOCAL VARS %
%%%%%%%%%%%%%%
distance_thresh = 0.01;          %default distance theshold
scale_method = 'bicubic';   %other choices are bilinear or nearest


display_images = false;


% CODE START %
%%%%%%%%%%%%%%
if nargin < 1 || nargin > 3
    error('incorrect number of arguments passed!');
elseif nargin >= 2
    refine_list = rl;
    if nargin == 3
        distance_thresh = st;
    end
end

num_clusts = size(Clust, 1);
chg_list = [];

%initially put all elements in the refine if not passed above
%smallest)
if nargin < 2
    refine_list = num_clusts:-1:1;
end
keep_list = 1:num_clusts;

while ~ isempty(refine_list)

    %determine if the first item in the list can be refined by up/down
    %sampling it to the size of the item to check
    r = refine_list(1);
    fprintf('                                                              \r');
    fprintf('cluster: %d -- ', r);
    rsz = size(Clust(r).avg);
    for i=keep_list

        if i == r
            continue;
        end

        isz = size(Clust(i).avg);

        %determine and scale the cluster average corresponding to rsz
        SclR = imresize(Clust(r).avg, isz, scale_method);
        
        [match, Mi, Mr, d] = euc_match(Clust(i).avg, SclR, distance_thresh);
        if match
            fprintf('found a match between cluster %d and %d\r', r, i);
            if display_images
                clf;
                subplot(1,2,1), imshow(Clust(r).avg), xlabel('r');
                subplot(1,2,2), imshow(Clust(i).avg), xlabel('i');
                drawnow;
                pause(2);
            end

            %now add the scaled component to the normal-sized cluster
            Clust(i) = add_and_reaverage(Clust(i), Clust(r), Mi, Mr);
            pos = find(keep_list == r);
            keep_list = [keep_list(1:pos-1), keep_list(pos+1:end)];
            chg_list = [chg_list, i];
            break;
        end
    end
    refine_list = refine_list(2:end);
end
fprintf('\n');

%refinement complete, now delete those items not on the keep_list
[Dummy, Dummy, chg_list] = intersect(unique(chg_list), keep_list);
Clust = Clust(keep_list);
