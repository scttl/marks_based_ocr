function evaluate_clusters(Clust, Comps, varargin)
% EVALUATE_CLUSTERS Get statistics on incorrectly clustered components
%
% evaluate_clusters(CLUST, COMPS, [VAR1, VAL1]...)
%
% This function uses ground truth component labels to gather statistics on how
% well the Components were clustered together.
%
% Clust is a struct like that returned from cluster_comps().  Comps is a struct 
% like that returned from get_comps().
%
% Optional arguments specified in LOCAL VARS below can be overridden

% CVS INFO %
%%%%%%%%%%%%
% $Id: evaluate_clusters.m,v 1.1 2006-10-18 15:56:12 scottl Exp $
%
% REVISION HISTORY
% $Log: evaluate_clusters.m,v $
% Revision 1.1  2006-10-18 15:56:12  scottl
% initial check-in.
%

% LOCAL VARS %
%%%%%%%%%%%%%%
ligature_val = 255;

%show we draw (and possibly) save an image of the mismatch elements?
draw_mismatches = false;
save_mismatches = false;

% CODE START %
%%%%%%%%%%%%%%
if nargin < 2
    error('must specify cluster and component structs!');
elseif nargin > 2
    process_optional_args(varargin{:});
end

if ~ Comps.found_true_labels
    error('can only gather stats if component ground truth labels exist');
end

fprintf('%d components\n', Comps.max_comp);
fprintf('%d clusters\n', Clust.num);

fprintf('%d ligature components\n', sum(Comps.truth_label == ligature_val));

misclust_count = 0;
misclust_idx = [];
misclust_comps = [];
clust_labels = zeros(Clust.num,1);
for ii=1:Clust.num
    cc = Clust.comps{ii};
    lbls = Comps.truth_label(cc);
    idx = find(lbls ~= ligature_val);
    nonlig_lbls = double(lbls(idx));
    cc = cc(idx);
    if ~ isempty(nonlig_lbls)
        clust_labels(ii) = mode(nonlig_lbls);
        idx = find(nonlig_lbls ~= clust_labels(ii));
        if ~isempty(idx)
            misclust_idx = [misclust_idx; ii];
            misclust_comps = [misclust_comps; cc(idx)];
            misclust_count = misclust_count + length(idx);
            fprintf('Clust %d misclassified components (mode: %c):\n', ...
                    ii, char(clust_labels(ii)));
            fprintf('   %d --> %c\n', [cc(idx)'; char(nonlig_lbls(idx))']);

        end
    end
end
fprintf('%d total misclustered non-lig components\n', misclust_count);

if draw_mismatches && misclust_count ~= 0
    fprintf('displaying mismatched clusters\n');
    display_cluster_elements(Clust, Comps, misclust_idx, 'comp_list', ...
         misclust_comps, 'display_avg', true, 'save_elements', save_mismatches);
    if save_mismatches
        fprintf('saved mismatched clusters to disk\n');
    end
end
