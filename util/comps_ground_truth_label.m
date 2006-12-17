function Comps = comps_ground_truth_label(Comps, varargin)
% COMPS_GROUND_TRUTH_LABEL  type in the character label of images.
%
% Comps = comps_ground_truth_label(Comps, [VAR1, VAL1]...)
%
% This function can be used to hand-label each of the Components passed, for
% later use as training data, or estimating performance.
%
% Comps is a struct like that returned from get_comps.
%
% Optional arguments specified in LOCAL VARS below can be overridden

% CVS INFO %
%%%%%%%%%%%%
% $Id: comps_ground_truth_label.m,v 1.2 2006-12-17 19:59:00 scottl Exp $
%
% REVISION HISTORY
% $Log: comps_ground_truth_label.m,v $
% Revision 1.2  2006-12-17 19:59:00  scottl
% renamed from ground_truth_label.m, changed to storing strings of
% characters (instead of singles), so that multiple character blobs can
% be correctly stored.
%
% Revision 1.1  2006/10/18 15:56:12  scottl
% initial check-in.
%

% LOCAL VARS %
%%%%%%%%%%%%%%
start_comp = 1;


% CODE START %
%%%%%%%%%%%%%%
if nargin < 1
    error('must specify the components struct!');
elseif nargin > 1
    process_optional_args(varargin{:});
end

if start_comp == 1
    Comps.truth_label = cell(Comps.max_comp, 1);
end

for ii=start_comp:Comps.max_comp
    display_neighbours(Comps, ii);
    Comps.truth_label(ii) = input(...
              'enter the characters that represent this component\n', 's');
    if rem(ii,30) == 0
        save 'truth_comps.mat' Comps ii;
    end
end
Comps.found_true_labels = true;
