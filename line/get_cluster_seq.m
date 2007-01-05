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
% $Id: get_cluster_seq.m,v 1.2 2007-01-05 17:15:52 scottl Exp $
%
% REVISION HISTORY
% $Log: get_cluster_seq.m,v $
% Revision 1.2  2007-01-05 17:15:52  scottl
% added ability to limit lines to particular regions.  Fix potential bug
% leading to infinite loop if a neighbour is self-referential.
%
% Revision 1.1  2006-10-29 17:14:43  scottl
% initial check-in.
%



% LOCAL VARIABLES %
%%%%%%%%%%%%%%%%%%%

%this is used to further limit the components returned to those that lie 
%entirely within the contained LTRB region on each page.  This is only really
%useful if all lines to process lie on the same page (as the same regions are
%used on each page.  Use a value of 0 to represent one of the extremities of 
%the page
keep_region = [0 0 0 0];


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
keep_lines = [];

for ii=1:length(line_nums)
    comp_list = find(Comps.line == line_nums(ii));
    if isempty(comp_list)
        seq{ii} = [];
        keep_lines = [keep_lines, ii];
        continue;
    end
    [idx, idx] = min(Comps.pos(comp_list, 1));
    idx = comp_list(idx);
    if (keep_region(2) ~= 0 && Comps.pos(idx,2) < keep_region(2)) || ...
       (keep_region(4) ~= 0 && Comps.pos(idx,4) > keep_region(4))
        %this line lies above or below the region boundary.  Skip it
        continue;
    end
    while keep_region(1) ~= 0 && idx ~= 0 && Comps.pos(idx,1) < keep_region(1)
        idx = Comps.nb(idx,3);
    end
    if idx == 0
        %entire line lies left of the left region boundary.  Skip it
        seq{ii} = [];
        keep_lines = [keep_lines, ii];
        continue;
    end
    out_list = idx;
    while Comps.nb(out_list(end),3) ~= 0 && (keep_region(3) == 0 || ...
          Comps.pos(out_list(end),3) <= keep_region(3))
        if out_list(end) == Comps.nb(out_list(end),3)
            %repeated neighbours, skip to the next neighbour listing this one
            next = find(Comps.nb(:,1) == out_list(end));
            if isempty(next)
                break;
            else
                out_list = [out_list, next(end)];
            end
        else
            out_list = [out_list, Comps.nb(out_list(end),3)];
        end
    end
    seq{ii} = Comps.clust(out_list);
    keep_lines = [keep_lines, ii];
end
seq = seq(keep_lines);


% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
