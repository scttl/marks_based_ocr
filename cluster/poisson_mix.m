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
%   (instead of maximizing the original)
%
%   Theta should be a column vector containing the 3 parameter values: l1, l2, 
%   and c (in that order)
%
%   S should be a vector of space intervals (the s_i's)
%
%   val, is the value of this function given the parameters
%
%   df is the vector of partial derivatives. ie. df = [df/dl1; df/dl2; df/dc]
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: poisson_mix.m,v 1.1 2007-03-06 20:12:10 scottl Exp $
%
% REVISION HISTORY
% $Log: poisson_mix.m,v $
% Revision 1.1  2007-03-06 20:12:10  scottl
% initial revision.
%


% LOCAL VARS %
%%%%%%%%%%%%%%
max_intrvl_len = 170;  %larger vals will have problems when factorial is taken


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
l1 = Theta(1);
l2 = Theta(2);
c = Theta(3);

%cleanup the input to ensure we have a sequence of positive integers that we
%can take the factorial of
S = round(S);
S(S <= 0) = 1;
S(S > max_intrvl_len) = max_intrvl_len;

%initialize the output and precompute values required at multiple points
num = length(S);
df = zeros(3,1);
ll1 = log(l1);
ll2 = log(l2);
sumS = sum(S);
ecs = exp(c - S);
ex = exp(S - c - l1 + S .* ll1 + l2 - S .* ll2);
ex1 = (ex + 1);
sumsigex = sum(ex ./ ex1);
sumsigSex = sum((S .* ex) ./ ex1);

%evaluate the function
val = num * l2 + sum(log(1 + exp(c - S))) + sum(log(factorial(S))) - ...
      ll2 * sumS - sum(log(ex1));

%evaluate the derivative with respect to l1
df(1) = sumsigex - 1/l1 * sumsigSex;

%evaluate the derivative with respect to l2
df(2) = num - 1/l2 * (sumS - sumsigSex) - sumsigex;

%evaluate the derivative with repsect to c
df(3) = sum(ecs ./ (1+ecs)) + sumsigex;
