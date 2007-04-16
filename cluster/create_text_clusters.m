function [Clust, Comps, Lines] = create_text_clusters(Files, varargin)
% CREATE_TEXT_CLUSTERS   Create cluster data from ASCII text files
%
%   [CLUST, COMPS, LINES] = create_text_clusters(FILES, [VAR1, VAL1]...)
%   FILES should be either a string or cell array listing the full path and 
%   name of the plain-text files to be used in clustering.
%
%   CLUST will be a struct like that returned in cluster_comps()
%   COMPS will be a struct like that returned in get_comps()
%   LINES will be a struct like that returned in get_lines()
%
%   VAR1 and VAL1 are optional, and can be used to override the default values
%   for the LOCAL VARS defined below.  Each VAR1 should be a string giving the
%   exact name of the variable to override.  Each VAL1 should be the updated
%   value to use for that variable.
%


% CVS INFO %
%%%%%%%%%%%%
% $Id: create_text_clusters.m,v 1.3 2007-04-16 03:30:22 scottl Exp $
%
% REVISION HISTORY
% $Log: create_text_clusters.m,v $
% Revision 1.3  2007-04-16 03:30:22  scottl
% convert logical array to true or false to make things clearer
%
% Revision 1.2  2007-02-05 21:33:45  scottl
% added density field.
%
% Revision 1.1  2007-02-01 18:04:48  scottl
% initial revision
%


% LOCAL VARS %
%%%%%%%%%%%%%%

%by defeault turn-off the override warning display since it will be generated
%for each line we create.
override_display = 'OFF';   %other legal value is 'ON'

%up to what word length should positional counts be taken?
max_word_len = 10;

%default font and size to use for density calculations
density_img_font = 'times-roman';
density_img_font_size = '36';


% CODE START %
%%%%%%%%%%%%%%
tic;
%process input arguments
if nargin < 1
    error('Must specify files to cluster');
elseif nargin > 1
    process_optional_args(varargin{:});
end
warning(override_display, 'MBOCR:override');

%initialize the output structs (note that not every field will be populated
Clust = init_clust();
Comps = init_comps(Files);
Lines = init_lines(Comps.files);

%this will hold the mapping from characters to cluster indices
clust_map = [];

%iterate through each page adding the necessary entries to each struct
for pp=1:length(Comps.files)

    fprintf('Processing page %d\n', pp);
    fid = fopen(Comps.files{pp});
    if fid == -1
        error('unable to open file');
    end
    C = fread(fid);
    fclose(fid);
    line_idcs = find(C == 10);  %10 == newline

    %fill in line properties
    Lines.pg(Lines.num+1:Lines.num+length(line_idcs)) = pp;
    curr_line = Lines.num+1;
    Lines.num = Lines.num + length(line_idcs);

    num_new_comps = length(C) - length(line_idcs);
    Comps.clust = [Comps.clust; zeros(num_new_comps,1,'uint16')];
    Comps.pg = [Comps.pg; pp+zeros(num_new_comps,1,'uint32')];
    Comps.nb = [Comps.nb; zeros(num_new_comps,4)];
    Comps.pos = [Comps.pos; inf+zeros(num_new_comps,4, 'uint16')];
    Comps.line = [Comps.line; zeros(num_new_comps,1,'uint64')];
    first_comp = true;
    for ii=1:length(C);
        if ~isempty(line_idcs) && ii == line_idcs(1)
            curr_line = curr_line + 1;
            first_comp = true;
            line_idcs = line_idcs(2:end);
            continue;
        end
        Comps.max_comp = Comps.max_comp + 1;
        cc = Comps.max_comp;

        %update cluster properties
        if C(ii) > length(clust_map) || clust_map(C(ii)) == 0
            %first time seeing this cluster
            Clust.num = Clust.num + 1;
            this_clust = Clust.num;
            Clust.num_comps(this_clust,:) = 1;
            Clust.comps{this_clust} = cc;
            Clust.truth_label{this_clust} = char(C(ii));
            clust_map(C(ii)) = this_clust;
        else
            %add this component to this cluster's list
            this_clust = clust_map(C(ii));
            Clust.num_comps(this_clust,:) = Clust.num_comps(this_clust,:) + 1;
            Clust.comps{this_clust} = [Clust.comps{this_clust}; cc];
        end

        %fill in component properties
        Comps.clust(cc) = this_clust;
        if ~first_comp
            Comps.nb(cc,1) = cc-1;
        else
            Comps.pos(cc,1) = 0;  %to ensure its the left most component
        end
        if ~isempty(line_idcs) && ii+1 ~= line_idcs(1)
            Comps.nb(cc,3) = cc+1;
        end
        Comps.line(cc) = curr_line;
        Comps.truth_label{cc} = char(C(ii));
        first_comp = false;
    end
end

%fill out the non-used fields to ensure sorting etc. is ok
Clust.mode_num = zeros(Clust.num,1);
Clust.avg = cell(Clust.num,1);
Clust.norm_sq = zeros(Clust.num,1);
Clust.refined = true(Clust.num,1);
Clust.changed = false(Clust.num,1);
Clust.ascender_off = zeros(Clust.num,1,'int16');
Clust.descender_off = zeros(Clust.num,1,'int16');
if any(strcmp(Clust.truth_label, ' '))
    Clust.model_spaces = true;
    Comps.model_spaces = true;
end
Clust.class = assign_class(Clust.truth_label);
Clust.density = assign_density(char(Clust.truth_label), 'img_font', ...
                density_img_font, 'img_font_sz', density_img_font_size);

%add cluster position counts
[Clust, Comps] = create_cluster_dictionary(Clust, Comps, 'max_word_len', ...
                 max_word_len);

%sort the clusters.  
[Clust, Comps] = sort_clusters(Clust, Comps);

%if warnings have been turned off, turn them back on
if strcmp(override_display, 'OFF')
    warning('ON', 'MBOCR:override');
end

fprintf('TOTAL NUMBER OF COMPONENTS NOW: %d\n', Comps.max_comp);
fprintf('TOTAL NUMBER OF CLUSTERS FOUND: %d\n', Clust.num);



% SUBFUNCTION DECLARATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize an empty cluster struct
function Clust = init_clust()
Clust.num = 0;
Clust.num_comps = [];
Clust.mode_num = [];
Clust.comps = {};
Clust.avg = {};
Clust.norm_sq = [];
Clust.refined = logical([]);
Clust.changed = logical([]);
Clust.found_offsets = false;
Clust.descender_off = int16([]);
Clust.ascender_off = int16([]);
Clust.num_trans = 0;
Clust.bigram = [];
Clust.pos_count = {};
Clust.pos_total = 0;
Clust.found_true_labels = true;
Clust.truth_label = {};
Clust.model_spaces = false;
Clust.class = uint16([]);
Clust.density = [];

%this function creates an empty component stucture
function Comps = init_comps(Files)
Comps.max_comp = 0;
if ~iscell(Files)
    Comps.files = {Files};
else
    Comps.files = Files;
end
Comps.pg_size = uint16([]);
Comps.clust = uint16([]);
Comps.pos = uint16([]);
Comps.pg = uint32([]);
Comps.nb = [];
Comps.nb_dist = uint16([]);
Comps.regions = [];
Comps.found_lines = true;
Comps.line = uint64([]);   %uint32([]);
Comps.modal_height = NaN;
Comps.scale_factor = [];
Comps.descender_off = int16([]);
Comps.ascender_off = int16([]);
Comps.found_true_labels = true;
Comps.truth_label = cell(0);
Comps.model_spaces = false;

%this function creates an empty line structure
function Lines = init_lines(files)
Lines.num = 0;
Lines.pg = uint32([]);
Lines.files = files;
Lines.pos = uint16([]);
Lines.baseline = uint16([]);  %must be same class as pos, so +,- work
Lines.xheight = uint16([]);
