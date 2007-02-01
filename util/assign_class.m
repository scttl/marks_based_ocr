function class_idcs = assign_class(dsc_offs, asc_offs, varargin)
% ASSIGN_CLASS   Categorizes the line offsets or symbols passed into classes
%
%   CLASS_IDCS = ASSIGN_CLASS(DCS_OFFS, ASC_OFFS, [VAR1, VAL1]...)
%
%   DCS_OFFS can be either a vector listing pixel baseline offsets (can be +ve
%   or -ve), or it can be a cell array listing characters to be classified
%
%   if DCS_OFFS is a vector, ASC_OFFS must also be specified as a vector of the
%   same length, listing x-line offsets (can be +ve or -ve).  If DCS_OFFS is a 
%   cell array of symbol values, this parameter's value is ignored.
%
%   VAR1 and VAL1 are optional, and can be used to override the default values
%   for the LOCAL VARS defined below.  Each VAR1 should be a string giving the
%   exact name of the variable to override.  Each VAL1 should be the updated
%   value to use for that variable.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: assign_class.m,v 1.1 2007-02-01 18:04:02 scottl Exp $
%
% REVISION HISTORY
% $Log: assign_class.m,v $
% Revision 1.1  2007-02-01 18:04:02  scottl
% initial revision.
%
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%how many different classes are there?
num_classes = 4;

%what symbols belong to which classes?
%the "ascenders"
class_1 = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZbdfhkl?![](){}@%$#/"^<>';
class_1_min_asc_thresh = 5;
%the "descenders"
class_2 = 'gjpqy';
class_2_min_dsc_thresh = 5;
%symbols that are short and near the baseline
class_3 = '.,_';
class_3_max_asc_thresh = -5;
%the rest (lie mostly between the baseline and x-height line)
class_def = ' aceimnorstuvwxz:;&-=*+~'; %all remainder that don't match above


% CODE START %
%%%%%%%%%%%%%%

%process input arguments
if nargin < 1
    error('must list either descender offsets or symbol values');
elseif nargin > 2
    process_optional_args(varargin{:});
end

if isnumeric(dsc_offs)
    %ensure ascender offsets of the same length are specified
    if nargin < 2 || length(dsc_offs) ~= length(asc_offs)
        error('incorrect ascender offset values');
    end
    class_idcs = zeros(length(dsc_offs),1);
    rem = 1:length(dsc_offs);
    cl_1_idx = asc_offs(rem) >= class_1_min_asc_thresh;
    cl_1_idx = rem(cl_1_idx);
    class_idcs(cl_1_idx) = 1;
    rem = setdiff(rem, cl_1_idx);
    cl_2_idx = dsc_offs(rem) >= class_2_min_dsc_thresh;
    cl_2_idx = rem(cl_2_idx);
    class_idcs(cl_2_idx) = 2;
    rem = setdiff(rem, cl_2_idx);
    cl_3_idx = asc_offs(rem) <= class_3_max_asc_thresh;
    cl_3_idx = rem(cl_3_idx);
    class_idcs(cl_3_idx) = 3;
    %remainder fall into default class
else
    %assume character values are passed
    if ~iscell(dsc_offs)
        error('if passing symbols, should store one per cell array entry');
    end
    class_idcs = zeros(length(dsc_offs),1);
    for ii=1:length(dsc_offs)
        if any(dsc_offs{ii} == class_1)
            class_idcs(ii) = 1;
        elseif any(dsc_offs{ii} == class_2)
            class_idcs(ii) = 2;
        elseif any(dsc_offs{ii} == class_3)
            class_idcs(ii) = 3;
        end
        %otherwise we retain our default '0' value for the class
    end
end

