function [imgs, norms] = get_comp_imgs(Comps, indices, varargin)
% GET_COMP_IMGS   Return cell array of component images from indices passed
%
%   [IMGS, NORMS] = GET_COMP_IMGS(COMPS, INDICES, ['var1', new_val1]...)
%
%   COMPS is a struct meant to hold information on all the connected components
%   found, and should be like that returned via get_comps().
%
%   INDCS should either be a vector listing which components (from COMPS)
%   should be returned as an image.  We assume each entry in this vector
%   references a valid component in COMPS.
%
%   IMGS will then be a cell array containing binary images of each component
%
%   NORMS are the associated L2 norm squared of each image (useful for
%   Euclidean distance calculations)
%
%   Note: any of the LOCAL VARS below can be overriden by passing in the name
%   of the variable to override as a string, and its new value as the second
%   parameter.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: get_comp_imgs.m,v 1.1 2006-10-18 15:54:21 scottl Exp $
%
% REVISION HISTORY
% $Log: get_comp_imgs.m,v $
% Revision 1.1  2006-10-18 15:54:21  scottl
% initial check-in.
%


% LOCAL VARS %
%%%%%%%%%%%%%%
num_dirs = 8;  %for connected components labelling

fg_val = 1;  %pixel value for 'on' pixels

% CODE START %
%%%%%%%%%%%%%%

%process arguments
if nargin < 2
    error('must specify the Comps struct and a list of indices!');
elseif nargin > 2
    process_optional_args(varargin{:});
end

%argument sanity checking
if any(indices > Comps.max_comp | indices < 1)
    error('indices out of range of Comps struct');
end

imgs = cell(length(indices),1);
norms = NaN(length(indices),1);

[pgs,idx] = sort(Comps.pg(indices));
pos = Comps.pos(indices,:);

for pp=unique(pgs)'
    M = bwlabel(~imread(Comps.files{pp}), num_dirs);

    for ii=idx(pgs == pp)'
        l = pos(ii,1); t = pos(ii,2); r = pos(ii,3); b = pos(ii,4);
        imgs{ii} = zeros(b-t+1, r-l+1);

        %there is a slight chance more than one connected component exists
        %inside the bounding box (can happen after merging characters).
        %However we want to remove any accidental overlaps of other components,
        %so don't include any components that also exist in other bounding
        %boxes
        lbl_part = M(t:b,l:r);
        on_idx = find(lbl_part ~= 0);

        if isempty(on_idx)
            warning('MBOCR:EmptyComp', 'empty component found\n');
            norms(ii) = 0;
            continue;
        end

        on_vals = lbl_part(on_idx);
        if ~ all(on_vals == on_vals(1))
            ovr_vals = [];
            if l > 1
                ovr_vals = [ovr_vals;M(t:b,l-1)];
            end
            if t > 1
                ovr_vals = [ovr_vals;M(t-1,l:r)'];
            end
            if r < Comps.pg_size(pp,2)
                ovr_vals = [ovr_vals;M(t:b,r+1)];
            end
            if b < Comps.pg_size(pp,1)
                ovr_vals = [ovr_vals;M(b+1,l:r)'];
            end
            keep_vals = setdiff(on_vals, ovr_vals);
            keep_idx = logical(zeros(length(on_idx),1));
            for jj=1:length(keep_vals)
                keep_idx = keep_idx | on_vals == keep_vals(jj);
            end
            on_idx = on_idx(keep_idx);
        end
        imgs{ii}(on_idx) = fg_val;

        %the norm_squared value is proportional to the number of 'on' pixels 
        %since each gets an initial intensity of fg_val.
        norms(ii) = fg_val * length(on_idx);
    end
end
