function [val, df] = poisson_mix(Theta, S, varargin)
% POISSON_MIX Determine partial derivatives and value of a 2-Poisson mixture
%
%   [val, df] = POISSON_MIX(Theta, S, [VAR1, VAL1])
%
%   This function evaluates the 2-Poisson mixture model objective function,
%   returning its value (in val), and vector of partial derivatives (in df)
%   for use in Carl Rasmussen's minimize utility.
%
%   The model is described in Huang et Al's "Cryptogram Decoding for Optical
%   Character Recognition" paper.
%
%   The objective function is defined by the likelihood of the data and is
%   as follows:
%
%     f(l1,l2,c) = \prod_i sig(s_i,c)pos(s_i,l1) + (1-sig(s_i,c)pos(s_i,l2)
%   where
%     s_i represents the ith space interval (S is the vector of s_i's passed)
%     c represents the threshold beyond which a space interval constitutes an
%       interword space (values < c constitute intercharacter spaces)
%     sig(s_i,c) is a shifted sigmoid = 1/(1+e^{c-s_i})
%     pos(s_i,l) is the Poisson distribution = (e^{-l} l^{s_i})/s_i!
%
%   we actually work with the negative log of this function for numerical ease 
%   (still valid since monotonic), and to give us a function to minimize 
%   (instead of maximizing the original).
%
%   Theta should be a column vector containing the logarithm of the 3 
%   parameter values (to constrain them to be positive), lnl1, lnl2, lnc, 
%   (in that order)
%
%   S should be a vector of space intervals (the s_i's)
%
%   val, is the value of this function given the parameters
%
%   df is the vector of partial derivatives. ie. df = [df/dlnl1; df/dlnl2; 
%   df/dlnc]
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: poisson_mix.m,v 1.2 2007-04-10 15:46:20 scottl Exp $
%
% REVISION HISTORY
% $Log: poisson_mix.m,v $
% Revision 1.2  2007-04-10 15:46:20  scottl
% working implementation of Hunag space model implemented.
%
% Revision 1.1  2007-03-06 20:12:10  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%


% CODE START %
%%%%%%%%%%%%%%

if nargin < 2
    error('incorrect number of arguments passed');
elseif nargin > 2
    process_optional_args(varargin{:});
end
if length(Theta) ~= 3
    error('Theta must contain 3 parameter values');
end
ll1 = Theta(1);
ll2 = Theta(2);
lc = Theta(3);

%cleanup the input to ensure we have a sequence of positive integers that we
%can take the factorial of
S = round(S);
S(S <= 0) = 1;

%initialize the output and precompute values required at multiple points
num = length(S);
df = zeros(3,1);

if isinf(ll1)
    error('lambda1 is 0');
end
if isinf(ll2)
    error('lambda2 is 0');
end
sumS = sum(S);
c = exp(lc);
ecs = exp(c - S);
ex = exp(S - c - exp(ll1) + S*ll1 + exp(ll2) - S*ll2);
ex1 = (ex + 1);
sumsigex = sum(ex ./ ex1);
sumsigSex = sum((S .* ex) ./ ex1);
sumsigcex = sum((c * ex) ./ ex1);

%compute the log factorial of S in a stable manner
%easy to see that log(s_i!) = sum(log([1:s_i]))
logfactS =  zeros(num,1);
for ii=1:num
    logfactS(ii) = sum(log([1:S(ii)]));
end
sumlogfactS = sum(logfactS);

%evaluate the function
val = sumlogfactS + sum(log(1 + ecs)) - sum(log(exp(S*ll1 - exp(ll1)) + ...
      exp(c - S + S*ll2 - exp(ll2))));



%evaluate the derivative with respect to l1
df(1) = exp(ll1) * sumsigex - sumsigSex;

%evaluate the derivative with respect to l2
df(2) = num * exp(ll2) - sumS - exp(ll2) * sumsigex + sumsigSex;

%evaluate the derivative with respect to c
df(3) = sum((c*ecs) ./ (1+ecs)) - sum((c*exp(c-S+S*ll2-exp(ll2)))./ ...
        (exp(S*ll1-exp(ll1)) + exp(c-S+S*ll2-exp(ll2))));
