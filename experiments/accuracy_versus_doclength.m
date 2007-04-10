function accuracy_versus_doclength(Files, varargin)
% ACCURACY_VERSUS_DOCLENGTH  Collect stats and create plot of accuracy
%
%    accuracy_versus_doclength(FILES, [VAR1, VAL1]...)
%
%    FILES should be either a string or a cell array listing the input files to
%    process and collect statistics from.  Each file should be plain ASCII text
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: accuracy_versus_doclength.m,v 1.1 2007-04-10 17:54:38 scottl Exp $
% REVISION HISTORY
% $Log: accuracy_versus_doclength.m,v $
% Revision 1.1  2007-04-10 17:54:38  scottl
% initial check-in.
%

global MOCR_PATH;  %used to determine where to save results

% LOCAL VARS %
%%%%%%%%%%%%%%

%up to what length words should we include stats for
max_word_len = 15;

%at what document lengths (number of characters) should we take positional
%counts, and generate output
doc_len_intervals = [100:100:2000, 3000:1000:50000, 100000:50000:1000000];

%this file should point at a file containing the list of possible symbols to use
sym_in_file = [MOCR_PATH, '/data/input_utf8_syms.txt'];

%alternately if this file isn't empty, it will specify a Syms struct already created to be used for symbol information
syms_struct_file = '';

%where should output text, result reports, and plots be kept?
save_res = true;
out_dir = [MOCR_PATH, '/results/acc_v_doclen_15_diffsyms'];
save_figure = true;
figure_driver = '-depsc2'; %help print for other choices
figure_res = '-r300';  %output resolution DPI
figure_file_prefix = [out_dir, '/plot_'];
figure_file_suffix = '.eps';

%where should we generate output and reports if not saving results
work_dir = '/tmp/acc_v_doc_tmp';

%OCRtk defines a hard maximum on the number of characters that can be processed %at any time
ocr_max_filesize = 65000;

%should we run through word lookups to improve the mappings?
word_lookup = true;


% CODE START %
%%%%%%%%%%%%%%
if nargin < 1
    error('must specify at least one file to collect stats for');
elseif nargin > 1
    process_optional_args(varargin{:});
end

if save_res
    work_dir = out_dir;
end

if ~iscell(Files)
    Files = {Files};
end

%create the working dir if it doesn't alrady exist
unix(['mkdir -p ', work_dir]);

%first generate symbol counts from the files (or load an existing struct)
if ~isempty(syms_struct_file)
    load(syms_struct_file);
else
    Syms = create_alphabet(sym_in_file, 'corpora_files', Files, ...
           'use_srilm', false, 'max_word_len', max_word_len);
    %must fixup the normalization term
    for ii=1:length(Syms.pos_count)
        val = Syms.pos_count{ii} .* Syms.pos_total;
        norms = sum(val,2);
        norms(norms == 0) = 1;  %to prevent dividing by 0
        Syms.pos_norms{ii} = repmat(norms, 1, size(val,2));
        Syms.pos_count{ii} = val ./ Syms.pos_norms{ii};
    end
end

%now load the text files into memory
C = [];
for ii=1:length(Files)
    fid = fopen(Files{ii});
    if fid == -1
        warning('MBOCR:FileOpen', 'Unable to open file: %s\n', Files{ii});
    else
        C = [C; fread(fid)];
        fclose(fid);
    end
end

prev_idx = 1;
last_iter = false;

%now loop over each chunk, clustering
for ii=1:length(doc_len_intervals)
    this_len = doc_len_intervals(ii);
    if this_len > length(C)
        this_len = length(C);
        last_iter = true;
    end

    %write out this piece so it can be clustered
    in_file = [work_dir, '/gt_', sprintf('%04d', ii), '.txt'];
    out_file = [work_dir, '/clust_', sprintf('%04d', ii), '.txt'];
    fid = fopen(in_file, 'w');
    if fid == -1
        error('problem opening file: %s for writing', in_file);
    end
    count = fwrite(fid, C(prev_idx:this_len));
    if count ~= this_len - prev_idx + 1;
        error('problem writing to file: %s', in_file);
    end

    %cluster the piece
    [Clust, Comps, Lines] = create_text_clusters(in_file, ...
                            'max_word_len', max_word_len);
    %now renormalize positional counts
    for jj=1:length(Clust.pos_count);
        val = Clust.pos_count{jj} ./ Clust.pos_total;
        norms = sum(val,2);
        norms(norms == 0) = 1;  %to prevent dividing by 0
        Clust.pos_count{jj} = val ./ repmat(norms, 1, size(val,2));
    end

    %determine positional mappings
    [order, score] = positional_learn_mappings(Clust, Syms, ...
                     'dist_metric', 'euc', 'weight_per_symbol', false, ...
                     'weight_proportion', 0);

    %set the final mapping
    if word_lookup
        map = word_lookup_map(Clust, Comps, Syms, 'order', order, ...
              'restrict_order_to_class', false, 'calc_valid_acc', false);
    else
        map = cell2mat(order);
        map = map(:,1);
    end

    %print the resultant text to file
    print_ocr_text(1:Lines.num, Comps, Syms, map, ...
                  'display_text', false, 'save_results', true, ...
                  'save_file', out_file);

    %generate character and word accuracy reports by processing 
    %pieces of the chunk (if its too large), and combining the reports
    fid = fopen(out_file);
    if fid == -1
        error('unable to open clustered text output');
    end
    D = fread(fid);
    fclose(fid);
    count = 1;
    %accuracy can't handle long files so we must break them into temp files
    while count < length(D) 
        end_count = min(count+ocr_max_filesize, length(D));
        cl_file = '/tmp/tmp.cl_txt';
        gt_file = '/tmp/tmp.gt_txt';
        fidcl = fopen(cl_file, 'w');
        fidgt = fopen(gt_file, 'w');
        if fidcl ==1 || fidgt == 1
            error('problems opening tmp file');
        end
        fwrite(fidcl, D(count:end_count));
        fwrite(fidgt, C(count:end_count));
        fclose(fidcl); fclose(fidgt);
        cmd = ['accuracy ', gt_file, ' ', cl_file, ' /tmp/tmp_', ...
               num2str(count), '.char_rprt'];
        s = unix(cmd);
        if s ~= 0
            error('prob running accuracy. cmd: %s', cmd);
        end
        cmd = ['wordacc ', gt_file, ' ', cl_file, ' /tmp/tmp_', ...
               num2str(count), '.word_rprt'];
        s = unix(cmd);
        if s ~= 0
            error('prob running word accuracy. cmd: %s', cmd);
        end
        count = count + ocr_max_filesize + 1;;
    end
    %merge these temp reports into a final pair of reports
    char_rprts = dir('/tmp/tmp_*.char_rprt');
    word_rprts = dir('/tmp/tmp_*.word_rprt');
    char_rprt_list = '';
    word_rprt_list = '';
    if length(char_rprts) > 1
        for jj=1:length(char_rprts)
            char_rprt_list = [char_rprt_list,'/tmp/', char_rprts(jj).name, ' '];
            word_rprt_list = [word_rprt_list,'/tmp/', word_rprts(jj).name, ' '];
        end
        cmd = ['accsum ', char_rprt_list, ' > ', work_dir, '/', ...
               sprintf('%04d', ii), '.chartot_rprt'];
        s = unix(cmd);
        if s ~= 0
            error('prob running accsum. cmd: %s', cmd);
        end
        cmd = ['wordaccsum ', word_rprt_list, ' > ', work_dir, '/', ...
               sprintf('%04d', ii), '.wordtot_rprt'];
        s = unix(cmd);
        if s ~= 0
            error('prob running word accsum. cmd: %s', cmd);
        end
    else
        %just copy the files
        cmd = ['cp /tmp/', char_rprts(1).name, ' ', work_dir, ...
               '/', sprintf('%04d', ii), '.chartot_rprt'];
        s = unix(cmd);
        if s ~= 0
            error('prob running cp. cmd: %s', cmd);
        end
        cmd = ['cp /tmp/', word_rprts(1).name, ' ', work_dir, ...
               '/', sprintf('%04d',ii),  '.wordtot_rprt'];
        s = unix(cmd);
        if s ~= 0
            error('prob running cp. cmd: %s', cmd);
        end
    end

    %cleanup temp files
    unix('rm -f /tmp/*_rprt /tmp/tmp*_txt');
end

%create a plot showing the resultant accuracy for each document length, broken
%up by the type of symbol.
if save_figure
    unlv_char_ocr_analysis(work_dir, 'renorm_by_total', true);
    print(gcf, figure_driver, figure_res, ...
          [figure_file_prefix, 'char', figure_file_suffix]);
    figure;
    unlv_word_ocr_analysis(work_dir);
    print(gcf, figure_driver, figure_res, ...
          [figure_file_prefix, 'word', figure_file_suffix]);
else
    unlv_char_ocr_analysis(work_dir, 'renorm_by_total', true);
    figure;
    unlv_word_ocr_analysis(work_dir);
end

%if we aren't saving results, cleanup work dir
if ~save_res
    unix(['rm -rf ', work_dir]);
end
