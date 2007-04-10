%this script will use tesseract to attempt to recognize the UNLV ISRI
%OCR dataset (the long articles in the 'L' dataset in particular)

global MOCR_PATH;  %used to determine where to save results

%set the following line to true to record process
create_diary = true;
if create_diary
    diary_file = [MOCR_PATH, ...
    '/results/tess_unlv_L_deskew_all.diary'];
    if exist(diary_file)
        delete(diary_file);
    end
    diary(diary_file);
    diary on;
    fprintf('EXPERIMENT STARTED: %s\n', datestr(now));
end

run_comps=true;
run_lines=true;
run_ocr=true;
run_ocr_analysis=true;

%this file should point at a file containing the list of pages to run
pg_file = [MOCR_PATH, '/data/unlv_ocr/L/LONG_PAGES'];

%these determine where to find the pages listed above
pg_prefix = [MOCR_PATH, '/data/unlv_ocr/L/L_B/'];
pg_suffix = '.FF';  %use fine-mode fax to compare with Nagy paper

%these determine where to find the ground truth text files for OCR analysis
gt_prefix = [MOCR_PATH, '/data/unlv_ocr/L/L_GT/'];

%this should give the path to the base part of where results will be kept
res_base = [MOCR_PATH, '/results/tess_L_long_deskew'];
if ~exist(res_base, 'dir')
    [s,w] = unix(['mkdir -p ', res_base]);
    if s~=0
        error('problem creating dir: %s', res_base);
    end
end

%open and read the list of pages from the pg_file
imgs = textread(pg_file, '%s');
docs = unique(regexprep(imgs, '(\w*)\-.*', '$1'));
num_docs = length(docs);
%convert any  '-' in the listed files to '_'
imgs = regexprep(imgs, '-', '_');

%give the path to the tesseract executable
tess_exe = '${HOME}/src/tesseract_trunk/tesseract';

tic;
for ii=1:num_docs
    fprintf('%.2f: Processing document: %s\n', toc, docs{ii});
    res_dir = [res_base, '/', docs{ii}];
    if ~exist(res_dir, 'dir')
        fprintf('%.2f: Creating new dir\n', toc);
        [s,w] = unix(['mkdir -p ', res_dir]);
        if s~=0
            warning('MBOCR:NoDir', 'problem creating dir: %s\n', res_dir);
            continue;
        end
    end

    res_datafile = [res_dir, '/data.mat'];
    idx = strmatch(docs{ii}, imgs);
    Files = regexprep(imgs(idx), '(.*)', [pg_prefix, '$1', pg_suffix]);

    %get the components
    if run_comps
        Comps = get_comps(Files, 'min_elem_width',2, 'min_elem_height',2, ...
                'max_elem_width',100, 'max_elem_height',100, 'deskew_pages',...
                true);
        fprintf('components complete: %f\n', toc);
        save(res_datafile, 'Comps');
    end

    load(res_datafile);

    %determine line boundaries
    if run_lines
        [Lines, Comps] = get_lines(Comps, 'base_thresh', .2, ...
                         'xheight_thresh', .3);
        fprintf('lines complete: %f\n', toc);
        save(res_datafile, 'Comps', 'Lines');
    end

    load(res_datafile);

    %now recognize the regions
    if run_ocr
        pgs = unique(Comps.pg);
        for jj=1:length(pgs)
            rgns = find(Comps.regions(:,1) == pgs(jj));
            this_pg = imgs{idx(jj)};
            pg_img = imread(Comps.files{pgs(jj)});
            for kk=1:length(rgns)
                tr = Comps.regions(rgns(kk),3:6);
                in_file = '/tmp/tess_region.tiff';
                imwrite(pg_img(tr(2):tr(4), tr(1):tr(3)), in_file, 'TIFF');
                cmd = [tess_exe, ' ', in_file, ' /tmp/tess_tmp batch'];
                s = unix(cmd);
                if s ~= 0
                    error('prob running tesseract. cmd: %s', cmd);
                end
                txt_file = [this_pg, '.Z', sprintf('%02d', ...
                            Comps.regions(rgns(kk),2))];
                res_txtfile = [res_dir, '/', txt_file];
                cmd = ['mv /tmp/tess_tmp.txt ', res_txtfile];
                s = unix(cmd);
                if s ~= 0
                    error('prob moving outputfile. cmd: %s', cmd);
                end
                delete('/tmp/tess*');
            end
        end
        fprintf('Tesseract OCR complete: %f\n', toc);
    end

    load(res_datafile);

    %print out mapped ground truth to a text file, and determine accuracy stats
    if run_ocr_analysis
        pgs = unique(Comps.pg);
        for jj=1:length(pgs)
            lines = find(Lines.pg == pgs(jj));
            rgns = find(Comps.regions(:,1) == pgs(jj));
            this_pg = imgs{idx(jj)};
            for kk=1:length(rgns)
                txt_file = [this_pg, '.Z', sprintf('%02d', ...
                            Comps.regions(rgns(kk),2))];
                res_txtfile = [res_dir, '/', txt_file];
                res_char_rprtfile = [res_txtfile, '.char_rprt'];
                res_word_rprtfile = [res_txtfile, '.word_rprt'];
                cmd = ['accuracy ', gt_prefix, txt_file, ' ', ...
                       res_txtfile, ' ', res_char_rprtfile];
                s = unix(cmd);
                if s ~= 0
                    error('prob running accuracy. cmd: %s', cmd);
                end
                cmd = ['wordacc ', gt_prefix, txt_file, ' ', ...
                       res_txtfile, ' ', res_word_rprtfile];
                s = unix(cmd);
                if s ~= 0
                    error('prob running word accuracy. cmd: %s', cmd);
                end
            end
        end
        %combine all the report files in this directory into a single
        %cumulative report
        char_rprts = dir([res_dir, '/*.char_rprt']);
        word_rprts = dir([res_dir, '/*.word_rprt']);
        char_rprt_list = '';
        word_rprt_list = '';
        if length(char_rprts) > 1
            for jj=1:length(char_rprts)
                char_rprt_list = [char_rprt_list, res_dir, '/', ...
                                  char_rprts(jj).name, ' '];
                word_rprt_list = [word_rprt_list, res_dir, '/', ...
                                  word_rprts(jj).name, ' '];
            end
            cmd = ['accsum ', char_rprt_list, ' > ', res_dir, '/', ...
                   docs{ii}, '.chartot_rprt'];
            s = unix(cmd);
            if s ~= 0
                error('prob running accsum. cmd: %s', cmd);
            end
            cmd = ['wordaccsum ', word_rprt_list, ' > ', res_dir, '/', ...
                   docs{ii}, '.wordtot_rprt'];
            s = unix(cmd);
            if s ~= 0
                error('prob running wordaccsum. cmd: %s', cmd);
            end
        else
            %just copy the single file for the total count
            cmd = ['cp ', res_dir, '/', char_rprts(1).name, ' ', res_dir, ...
                   '/', docs{ii}, '.chartot_rprt'];
            s = unix(cmd);
            if s ~= 0
                error('prob running cp. cmd: %s', cmd);
            end
            cmd = ['cp ', res_dir, '/', word_rprts(1).name, ' ', res_dir, ...
                   '/', docs{ii}, '.wordtot_rprt'];
            s = unix(cmd);
            if s ~= 0
                error('prob running cp. cmd: %s', cmd);
            end
        end

        fprintf('ocr analysis complete: %f\n', toc);
    end

end

if create_diary
    fprintf('EXPERIMENT ENDED: %s\n', datestr(now));
    diary off;
end
