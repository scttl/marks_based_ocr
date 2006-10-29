function seq = get_cluster_seq(Comps, line_nums, varargin)
% GET_CLUSTER_SEQ  Get sequence of cluster id's for each line passed.
%
%   SEQ = GET_CLUSTER_SEQ(COMPS, LINE_NUMS, ['var1', val1]...)
%
%  Comps should be a struct containing component information, and should
%  include a line field listing which line each component belongs to.  It
%  should also include a clust field listing which cluster each component 
%  belongs to.
%
%  line_nums should be a scalar or vector listing the line numbers for which
%  the sequence of clusters should be returned (in left-right "reading" order)
%
%  seq returned will be a cell array, each entry contianing row vectors listing
%  the sequence of cluster id's


% CVS INFO %
%%%%%%%%%%%%
% $Id: get_cluster_seq.m,v 1.1 2006-10-29 17:14:43 scottl Exp $
%
% REVISION HISTORY
% $Log: get_cluster_seq.m,v $
% Revision 1.1  2006-10-29 17:14:43  scottl
% initial check-in.
%



% LOCAL VARIABLES %
%%%%%%%%%%%%%%%%%%%



% CODE START %
%%%%%%%%%%%%%%
if nargin < 2
    error('must pass comps struct and listing of line nums');
elseif nargin > 2
    process_optional_args(varargin{:});
end

if ~isfield(Comps, 'found_lines') || ~ Comps.found_lines
    error('Comps.line must exist and contain lines');
end

if size(line_nums,1) > 1
    line_nums = line_nums';
end

seq = cell(length(line_nums), 1);

for ii=1:length(line_nums)
    comp_list = find(Comps.line == line_nums(ii));
    if isempty(comp_list)
        seq{ii} = [];
        continue;
    end
    [idx, idx] = min(Comps.pos(comp_list, 1));
    out_list = comp_list(idx);
    while Comps.nb(out_list(end),3) ~= 0
        out_list = [out_list, Comps.nb(out_list(end),3)];
    end
    seq{ii} = Comps.clust(out_list);
end


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
