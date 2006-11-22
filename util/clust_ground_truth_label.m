function Clust = clust_ground_truth_label(Clust, Comps, varargin)
% CLUST_GROUND_TRUTH_LABEL  Label images for cluster ground truth data
%
% Clust = clust_ground_truth_label(Clust, Comps, [VAR1, VAL1]...)
%
% This function can be used to hand-label each of the Clusters passed, for
% later use as training data, or estimating performance.
%
% Clust is a struct like that returned from cluster_comps.
% Comps is a struct like that returned from get_comps.
%
% Optional arguments specified in LOCAL VARS below can be overridden
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: clust_ground_truth_label.m,v 1.1 2006-11-22 17:23:35 scottl Exp $
%
% REVISION HISTORY
% $Log: clust_ground_truth_label.m,v $
% Revision 1.1  2006-11-22 17:23:35  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%
start_clust = 1;

%how often should we write the truth labels out to disk (set to [] to disable)
save_every = 30;
global MOCR_PATH;
save_file = [MOCR_PATH, '/results/truth_clust.mat'];

% CODE START %
%%%%%%%%%%%%%%
if nargin < 2
    error('must specify the cluster and components structs!');
elseif nargin > 2
    process_optional_args(varargin{:});
end

if start_clust == 1
    Clust.truth_label = cell(Clust.num, 1);
end

for ii=start_clust:Clust.num
    subplot(1,2,1); imshow(Clust.avg{ii});
    title(['Cluster ', num2str(ii)]);
    subplot(1,2,2); display_neighbours(Comps, Clust.comps{ii}(1));
    xx = input('enter the keyboard symbol that represents this symbol\n', ...
              's');
    Clust.truth_label{ii} = xx;
    if ~isempty(save_every) && rem(ii,save_every) == 0
        save(save_file, 'Clust', 'ii');
    end
end
Clust.found_true_labels = true;
