function cluster_err_analysis(Clust, Syms, map, varargin)
% CLUSTER_ERR_ANALYSIS    Plot positional features of incorrectly predicted sym
%
%     cluster_err_analysis(CLUST, SYMS, MAP, [VAR1, VAL1]...)
%
% CLUST and SYMS should be structs as returned by cluster_comps and
% create_alphabet respectively.
% 
% NOTE: we assume CLUST has already been sorted based on component occurence
%
% MAP should be the cluster symbol index mapping of the closest match
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: cluster_err_analysis.m,v 1.1 2007-02-01 18:04:02 scottl Exp $
%
% REVISION HISTORY
% $Log: cluster_err_analysis.m,v $
% Revision 1.1  2007-02-01 18:04:02  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%should we save the plot we generate?
save_plot = true;
plot_file_prefix = ['./plot_'];
plot_file_suffix = '.eps';
plot_driver = '-depsc2'; %help print for other choices
plot_res = '-r300';  %output resolution DPI

%how many of the mismatches should we plot?
num_to_plot = 3;

%we can optionally pass in an ordering and score to determine where the 
%correct match actually occured, and how far the scores deviated
order = [];
score = [];

%by overridding the value for this parameter, you can restrict mismatches to
%particular clusters.  Useful for just considering digit or upper case character
%mismatches etc.
r_clust = [];

% CODE START %
%%%%%%%%%%%%%%

if nargin < 4
    error('incorrect number of arguments passed');
elseif nargin > 4
    process_optional_args(varargin{:});
end

if isempty(r_clust)
    r_clust = 1:Clust.num;
end
mismatches = find(~strcmp(Syms.val(map(r_clust)), Clust.truth_label(r_clust)'));
mismatches = r_clust(mismatches);

if mismatches == 0
    fprintf('no mismatches!\n');
    return;
elseif mismatches < num_to_plot
    num_to_plot = mismatches;
end

cl_cts = cell2mat(Clust.pos_count);
sm_cts = cell2mat(Syms.pos_count);

for ii=1:num_to_plot
    cl_idx = mismatches(ii);
    corr_sym = Clust.truth_label{cl_idx};
    gen_sym = Syms.val{map(cl_idx)};
    gen_sym_idx = map(cl_idx);

    %locate the index of the correct symbol
    corr_sym_idx = strmatch(corr_sym, Syms.val, 'exact');
    if isempty(corr_sym_idx)
        error('could not find symbol index');
    end
    corr_sym_idx = corr_sym_idx(1);

    %calculate some error statistics
    fprintf('Incorrect symbol %s.  Correct symbol %s\n', gen_sym, corr_sym);
    fprintf('    generated %d times. For %.4f document percent\n',...
            Clust.num_comps(cl_idx), ...
            (Clust.num_comps(cl_idx)/sum(Clust.num_comps(2:end))*100));
    if ~isempty(order)
        pos = find(order{cl_idx} == corr_sym_idx);
        fprintf('    actual symbol was the %d closest positional match\n', pos);
        if ~isempty(score)
            fprintf('    predicted distance: %.4f, Actual distance: %.4f\n', ...
            score{cl_idx}([1; pos]));
        end
        fprintf('\n');
    end

    %plot the feature vectors
    ft = [cl_cts(cl_idx,:); sm_cts(gen_sym_idx,:); sm_cts(corr_sym_idx,:)];
    lgnd_str = {['Cluster: ''', corr_sym, ''''];
                ['Symbol Generated: ''', gen_sym, ''''];
                ['True Symbol: ''', corr_sym, ''''];};
    
    pos_feature_plot(ft, 'plot_marker_string', '+:', ...
                     'plot_legend_strings', lgnd_str);
    if save_plot
        fprintf('saving plot to disk\n');
        print(gcf, plot_driver, plot_res, [plot_file_prefix, gen_sym, '_', ...
              corr_sym, plot_file_suffix]);
    end
end

