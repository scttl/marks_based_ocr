%this script will load about 5 articles (each with ~ 8 pages of images) from
%the NIPS 2001 training set, then cluster the blobs of ink found on each page,
%determine the lines of the page, and create training data out of each of them
%(both the lines, and the cluster counts for the "language model")

%diary('../results/nips5_articles.diary');
%diary on;

run_cluster=true;
run_training_data=true;
run_dictionary=true;
run_ocr=true;

Files = ls('/p/learning/scottl/research/TRAINING_DATA/infinity_book/*.tiff');

%must convert Files from a single long newline delimited string into a cell 
%array
len = [];
idx = strfind(Files, char(10));
for ii=1:length(idx)
    if ii == 1
        len(ii) = idx(ii);
    else
        len(ii) = idx(ii) - idx(ii-1);
    end
end
Files = strtrim(mat2cell(Files, 1, len)');

%now cluster the components
tic;
if run_cluster
    [Clust, Comps] = cluster_comps(Files);
    fprintf('clustering complete: %f\n', toc);
    save('../results/infinity_book.mat', 'Clust', 'Comps');
end

load('../results/infinity_book.mat');

%now build up a set of training cases (50)
if run_training_data
    imgs = create_cluster_training_data(Comps, 50);
    fprintf('training data images complete: %f\n', toc);
    save('../results/infinity_book.mat', 'Clust', 'Comps', 'imgs');
end

load('../results/infinity_book.mat');

%now create a language model from the same dataset
if run_dictionary
    [Clust, Comps] = create_cluster_dictionary(Clust, Comps);

    fprintf('dictionary complete: %f\n', toc);
    save('../results/infinity_book.mat', 'Clust', 'Comps', 'imgs');
end

load('../results/infinity_book.mat');
if run_ocr
    %try solving the training image lines on the top 50 clusters (characters)
    [Clust, Comps] = sort_clusters(Clust, Comps);
    [vals, segs] = do_ocr(imgs, Clust.avg(1:50), double(Clust.offset(1:50)), ...
                   Clust.bigram(1:50,1:50));
    fprintf('testing complete: %f\n', toc);
    save('../results/infinity_book.mat', 'Clust', 'Comps', 'imgs', 'vals', ...
         'segs');
end

diary off;
