function [char_seq,ll] = viterbi_decode(obs_sequence, trans_probs, ...
                         emit_probs, init_probs)
% VITERBI_DECODE Run Viterbi decoding to get most likely character sequence
%
%   [char_seq,log_likehood] = VITERBI_DECODE(obs_sequence, trans_probs, ...
%                             emit_probs, init_probs)
%
%   obs_sequence should be vector of cluster id's
%   trans_probs, emit_probs, and init_probs should be probability tables like
%   those returned from char_HMM_init.
%
%   the most likely corresponding character sequence is returned in char_seq
%   the score of the sequence is also returned in log_likelihood
%   


% CVS INFO %
%%%%%%%%%%%%
% $Id: viterbi_decode.m,v 1.1 2006-10-18 16:02:23 scottl Exp $
%
% REVISION HISTORY
% $Log: viterbi_decode.m,v $
% Revision 1.1  2006-10-18 16:02:23  scottl
% nitial check-in.
%

% LOCAL VARS %
%%%%%%%%%%%%%%

% CODE START %
%%%%%%%%%%%%%%
num_obs = length(obs_sequence); 
[num_clust,num_chars] = size(emit_probs); 

if(size(init_probs,2)~=1) 
    init_probs=init_probs(:); 
end

%ensure arguments are ok
if num_obs <= 0 || length(init_probs) ~= num_chars || ...
   max(obs_sequence) > num_clust || min(obs_sequence) < 1
   error('invalid parameter');
end

% setup output, intermediate terms that can be precalculated
delta=zeros(num_chars,num_obs);  
psi=zeros(num_chars,num_obs); 
char_seq=zeros(1,num_obs);

% grab the distribution over characters of our observed cluster sequence
char_dist = emit_probs(obs_sequence,:)';

%augment parameters so we can work in the log domain
I = log(init_probs+eps);
T = log(trans_probs+eps);
char_dist = log(char_dist+eps);

delta(:,1) = I+char_dist(:,1); 
psi(:,1)=0;

for tt=2:num_obs
  [delta(:,tt),psi(:,tt)] = max((delta(:,tt-1)*ones(1,num_chars)+T)',[],2);
  delta(:,tt) = delta(:,tt)+char_dist(:,tt);
end

[ll,char_seq(num_obs)] = max(delta(:,num_obs));

for tt=(num_obs-1):-1:1
  char_seq(tt)=psi(char_seq(tt+1),tt+1);
end

