%this script will load about 5 articles (each with ~ 8 pages of images) from
%the NIPS 2001 training set, then cluster the blobs of ink found on each page,
%determine the lines of the page, and create training data out of each of them
%(both the lines, and the cluster counts for the "language model")

global MOCR_PATH;  %used to determine where to save results
%diary([MOCR_PATH, '/results/nips_5_articles.diary']);
%diary on;

run_comps=true;
run_lines=true;
run_cluster=true;
run_training_data=true;
run_dictionary=true;
run_ocr=true;

Files = ls('~scottl/research/TRAINING_DATA/nips_2001/AA0[1-7]*.tif');


%must convert Files from a single long newline delimited string into a cell 
%array
idx = strfind(Files, char(10));
for ii=1:length(idx)
    if ii == 1
        len(ii) = idx(ii);
    else
        len(ii) = idx(ii) - idx(ii-1);
    end
end
Files = strtrim(mat2cell(Files, 1, len)');

%get the components
tic;
if run_comps
    Comps = get_comps(Files);
    fprintf('components complete: %f\n', toc);
    save([MOCR_PATH, '/results/nips_5_articles.mat'], 'Comps');
end

load([MOCR_PATH, '/results/nips_5_articles.mat']);

%determine line boundaries
if run_lines
    [Lines, Comps] = get_lines(Comps);
    fprintf('lines complete: %f\n', toc);
    save([MOCR_PATH, '/results/nips_5_articles.mat'], 'Comps', 'Lines');
end

load([MOCR_PATH, '/results/nips_5_articles.mat']);

%now cluster the components
if run_cluster
    if run_lines
        [Clust, Comps] = cluster_comps(Comps, 'Lines', Lines);
        save([MOCR_PATH, '/results/nips_5_articles.mat'], 'Clust', 'Comps', ...
             'Lines');
    else
        [Clust, Comps] = cluster_comps(Comps);
        save([MOCR_PATH, '/results/nips_5_articles.mat'], 'Clust', 'Comps');
    end
    fprintf('clustering complete: %f\n', toc);
end

load([MOCR_PATH, '/results/nips_5_articles.mat']);

%now build up a set of training cases (50)
if run_training_data
    imgs = create_cluster_training_data(Comps, 50);
    fprintf('training data images complete: %f\n', toc);
    save([MOCR_PATH, '/results/nips_5_articles.mat'], 'Clust', 'Comps', 'imgs');
end

load([MOCR_PATH, '/results/nips_5_articles.mat']);

%now create a language model from the same dataset
if run_dictionary
    [Clust, Comps] = create_cluster_dictionary(Clust, Comps);

    fprintf('dictionary complete: %f\n', toc);
    save([MOCR_PATH, '/results/nips_5_articles.mat'], 'Clust', 'Comps', 'imgs');
end

load([MOCR_PATH, '/results/nips_5_articles.mat']);

if run_ocr
    %try solving the training image lines on the top 50 clusters (characters)
    [Clust, Comps] = sort_clusters(Clust, Comps);
    [vals, segs, scores] = do_ocr(imgs, Clust.avg(1:50), ...
                           double(Clust.offset(1:50)), ...
                           Clust.bigram(1:50,1:50));
    save([MOCR_PATH, '/results/nips_5_articles.mat'], 'Clust', 'Comps', ...
         'imgs', 'vals', 'segs', 'scores');
end

%diary off;
