function Comps = ground_truth_label(Comps, varargin)
% GROUND_TRUTH_LABEL  Allow the user to type in the character label of images.
%
% Comps = ground_truth_label(Comps, [VAR1, VAL1]...)
%
% This function can be used to hand-label each of the Components passed, for
% later use as training data, or estimating performance.
%
% Comps is a struct like that returned from get_comps.
%
% Optional arguments specified in LOCAL VARS below can be overridden

% CVS INFO %
%%%%%%%%%%%%
% $Id: comps_ground_truth_label.m,v 1.1 2006-10-18 15:56:12 scottl Exp $
%
% REVISION HISTORY
% $Log: comps_ground_truth_label.m,v $
% Revision 1.1  2006-10-18 15:56:12  scottl
% initial check-in.
%

% LOCAL VARS %
%%%%%%%%%%%%%%
start_comp = 1;

ligature_val = 255;


% CODE START %
%%%%%%%%%%%%%%
if nargin < 1
    error('must specify the components struct!');
elseif nargin > 1
    process_optional_args(varargin{:});
end

if start_comp == 1
    Comps.truth_label = zeros(Comps.max_comp, 1, 'uint8');
end

for ii=start_comp:Comps.max_comp
    display_neighbours(Comps, ii);
    xx = input('enter the keyboard symbol that represents this character\n', ...
              's');
    if length(xx) > 1
        Comps.truth_label(ii) = ligature_val;
    else
        Comps.truth_label(ii) = uint8(xx);
    end
    if rem(ii,30) == 0
        save 'truth_comps.mat' Comps ii;
    end
end
Comps.found_true_labels = true;
