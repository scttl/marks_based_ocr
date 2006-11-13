function [trans, emit, init, d_idx] = char_HMM_init(D, clust_counts, varargin)
% CHAR_HMM_INIT Determine parameters for a character HMM model
%
%   [trans_prob, emit_prob, init_prob, chars] = CHAR_HMM_INIT(D, ....
%                clust_counts, [VAR1, VAL1]...)
%
%   D is a struct containing character frequency counts taken from a text
%   corpus (see create_dictionary)
%
%   clust_counts should be a vector listing the number of components belonging
%   to that cluster.  This should be passed in as the Clust.num_comps (see
%   cluster_comps, and the Clust struct for more info)
%
%   the transition probabilites and initial start state probabilites are
%   estimated directly from the counts in D (potentially with smoothing counts
%   added)
%
%   the emission probabilities are estimated by a gaussian distribution over
%   clusters (which are ordered by frequency).  The mean is set to
%   approximately the spot it which the underlying character appears (when
%   characters are ordered by frequency using counts from D).  The variance is
%   a tunable parameter.
%
%   NOTE: We will strip out non-printable characters like newlines
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: char_HMM_init.m,v 1.2 2006-11-13 17:57:23 scottl Exp $
%
% REVISION HISTORY
% $Log: char_HMM_init.m,v $
% Revision 1.2  2006-11-13 17:57:23  scottl
% no change.
%
% Revision 1.1  2006/10/18 16:02:23  scottl
% Initial check-in.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%smoothing counts to add
trans_pseudo_counts = 1;
init_pseudo_counts = 1;

%the standard deviation should cover what portion of the numbers
emit_std_dev_pct = .05;

strip_chars = [10];  %line feed characters are removed

% CODE START %
%%%%%%%%%%%%%%
tic;

if nargin < 2
    error('incorrect number of arguments');
elseif nargin > 2
    process_optional_args(varargin{:});
end

num_chars = length(D.char);
num_clust = length(clust_counts);

rem_idx = [];
chars = double(D.char);
for ii=1:length(strip_chars)
    rem_idx = [rem_idx, find(strip_chars(ii) == chars)];
end
d_idx = setdiff(1:num_chars, rem_idx);
num_chars = length(d_idx);

trans = D.char_bigram(d_idx,d_idx) + trans_pseudo_counts;
trans = trans ./ repmat(sum(trans,2),1,size(trans,2));
fprintf('%.2fs: transition probabilites computed\n');

init = D.first_count(d_idx) + init_pseudo_counts;
init = init ./ sum(init);
fprintf('%.2fs: initial state probabilites computed\n');


emit = zeros(num_clust, num_chars);

[ch_cnts,ch_idx] = sort(D.char_count(d_idx), 'descend');
[cl_cnts,cl_idx] = sort(clust_counts, 'descend');
sigma = emit_std_dev_pct * num_clust;
emit = emit + (1/(sigma * sqrt(2 * pi)));
for ii=1:num_chars
    %note: we fix the distribution of the most frequent character, since this
    %will correspond to the space character.  Its typically more than twice as
    %frequent as the next most frequent character.  Same for the cluster
    %frequencies
    if ch_idx(ii) == 1
        %find most frequent cluster index
        mfcl_idx = find(cl_idx == 1);
        emit(mfcl_idx,ii) = 1;
        emit([1:mfcl_idx-1,mfcl_idx+1:num_clust], ii) = 0;
    else
        %@@emit(:,ii) = ones(num_clust,1) .* (1/num_clust);
        mu = ceil(ch_idx(ii)/num_chars * num_clust);
        x = cl_idx;
        emit(:,ii) = emit(:,ii) .* exp(- ((x-mu).^2/(2*sigma^2)));
    end
end

%because we've discretized the distribution we must renormalize it to ensure
%everything adds up correctly
emit = emit ./ repmat(sum(emit), num_clust, 1);
