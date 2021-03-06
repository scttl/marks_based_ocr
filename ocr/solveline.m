function [bestpath,bestseg,bestscore]=solveline(data,models,bigram,delprob,...
                                      insprob,Wmin,Wmax,sl)
% SOLVELINE  Determine best sequence of character models to explain a line image
%
%  [bestpath,bestseg,bestscore] = SOLVELINE(data, models, bigram, delprob, 
%                                 insprob, Wmin, Wmax, [short_list])
%
%  data should be an image array representing a line of text to be analyzed.
%
%  models should be a cell array of character/symbol images, normalized and
%  background padded so that they each are the same height as the image in
%  data.
%
%  bigram should be a square matrix whose entries represent the transition
%  probabilities between characters, based on an underlying language model (like
%  English).  The number of characters must match the length of the models cell
%  array.
%
%  delprob and insprob, represent the score (or probability) of a mismatch.
%  Delprob refers to the score of having to remove a foreground pixel from the
%  model to match a corresponding background pixel in the data.  Likewise,
%  insporb refers to the score assigned to adding a foreground pixel to the
%  model to match a corresponding foreground pixel in the data.  Foreground
%  pixels in the model that match foreground pixels in the data are given a
%  score of 1-delprob.  Background pixels in the model that match background
%  pixels in the data are given a score of 1-insprob.  Since the log of these 
%  values will be taken they should lie in the range (0...1) to ensure scores 
%  remain well-defined.
%
%  Wmin and Wmax should be integers representing the minimum and additional 
%  maximum (beyond the width of the character) distance (in pixels) between 
%  the placement of successive model characters.
%
%  short_list is optional and if specified should be a cell array containing
%  column vectors (one such vector for every column of data passed) which lists
%  the indices of the models which are valid for that particular column.  If
%  not passed, all models will be considered for all columns.  Note that
%  specifying a small subset of valid models for each column vastly improves
%  the speed of this algorithm and is recommended.  Also note that any empty
%  vectors will indicate that no models match this column (for instance if we
%  are in a column between words).  Single valued vectors (scalars) will ensure
%  that only that particular model matches, thus constraining the solver.
%
%  bestpath should be a vector of character model indices, and bestseg should
%  be a vector of the same length, giving the columns in data at which to start
%  placing the corresponding model character from bestpath.  bestscore
%  similarly will be score of placing each character in bestpath at the
%  bestpath columns


% CVS INFO %
%%%%%%%%%%%%
% $Id: solveline.m,v 1.6 2006-09-16 22:31:46 scottl Exp $
%
% REVISION HISTORY
% $Log: solveline.m,v $
% Revision 1.6  2006-09-16 22:31:46  scottl
% implemented return of scores corresponding to the selected characters (based on
% their mismatch with the underlying data).  Useful for generating heatmap result
% display.
%
% Revision 1.5  2006/09/05 15:53:56  scottl
% implemented short-lists to constrain line solver.
%
% Revision 1.4  2006/08/30 17:36:33  scottl
% implemented ability to pass the post-line transition character as a parameter.
%
% Revision 1.3  2006/08/14 01:22:59  scottl
% suspended use of weighted penalties for now.
%
% Revision 1.2  2006/08/05 17:33:08  scottl
% added comment header.  Weighted penalities by their distance to closest match
% in the data.
%


% LOCAL VARS %
%%%%%%%%%%%%%%
end_val = 1;  %default model to transition to after completing the data.

% CODE START %
%%%%%%%%%%%%%%
if nargin < 7 || nargin > 8
    error('incorrect number of arguments specified');
elseif nargin == 8
    short_list = sl;
end

[hh,N]=size(data);
K=length(models);
if nargin ~= 8
    %construct a short list containing all models over each column of data
    short_list = mat2cell(repmat((1:K)', 1, N), K, ones(1,N));
end
if N ~= length(short_list)
    error('short_list length doesnt match number of columns of data!');
end
%determine the log insertion and deletion probabilities for each model as well
%as the maximum character width.
maxcharwidth=0;
mincharwidth=inf;
for kk=1:K
  maxcharwidth=max(maxcharwidth,size(models{kk},2));
  mincharwidth=min(mincharwidth,size(models{kk},2));
  [logmp{kk},logmpbar{kk}]=logmprob(models{kk},delprob,insprob);
end

if Wmin > mincharwidth + Wmax
    error('Wmin is too large (or Wmax is too small) for the model sizes!');
end

%initialize the costs for the image, ensuring that we end in the last column.
%This is accomplished by including a perfect score for a transition to the
%end_val model in the next column beyond the edge of the image.
costs=Inf(K,N+maxcharwidth+Wmax); costs(end_val,N+1)=0;
bestpaths=zeros(K,N);
bestdeltas=zeros(K,N);
bestscores=zeros(K,N);

data=[data,zeros(hh,maxcharwidth+Wmax)];
%the lines below assign increasing distance to pixels further away from the
%closest 'on' pixel (or 'off' pixel for insdist).  To prevent this behaviour
%uncomment the two lines below
insdist=data;
deldist=(1-data);
%insdist=bwdist(1-data);
%deldist=bwdist(data);

bgcost=-log(bigram);


%calculate the costs for each column, starting from the last.  
thispos=N;
while(thispos>=1)
  fprintf('Computing costs for column %d          \r',thispos);
  newcost=Inf(K,1);
  for cc=short_list{thispos}'
    thisww=size(logmp{cc},2);
    %determine cost of placing character cc at this column
    sc=scorecols(data(:,thispos:(thispos+thisww+Wmax-1)),...
       insdist(:,thispos:(thispos+thisww+Wmax-1)),...
       deldist(:,thispos:(thispos+thisww+Wmax-1)),...
       logmp{cc},logmpbar{cc},insprob);

    %add to this score the previous costs of completing the rest of the line,
    %plus the costs according to the bigram model.
    tmpcost=costs(:,(thispos+Wmin):(thispos+thisww+Wmax))+... %prev cost
            bgcost(cc,:)'*ones(1,thisww+Wmax-Wmin+1)+...      %bigram cost
            ones(K,1)*sc(Wmin:end);                           %model score

    %calculate the minimal cost, whose location gives rise to the best 
    %character to transition to as well as the offset (column) at which to 
    %place this transitional character.
    [newcost(cc),minidx]=min(tmpcost(:));
    [bestkk,bestrr]=ind2sub([K,thisww+Wmax-Wmin+1],minidx(1)); 
                                                  % ARBITRARILY TAKE FIRST!
    bestpaths(cc,thispos)=bestkk;
    bestdeltas(cc,thispos)=bestrr+Wmin-1;
    bestscores(cc,thispos)=sc(bestrr+Wmin-1);
  end
  costs(:,thispos)=newcost;
  thispos=thispos-1;
end

%Use the calculated costs to determine the best path and segmentation points in
%the data, starting by placing the first character in the first column.
[bestpath,bestpath]=min(costs(:,1));
bestseg=1;
bestscore=[];
while(bestseg(end)<N)
  bestscore=[bestscore,bestscores(bestpath(end),bestseg(end))];
  nextdelta=bestdeltas(bestpath(end),bestseg(end));
  nextcolumn=bestseg(end)+nextdelta;
  bestseg=[bestseg,nextcolumn];
  [tmpcost,choice]=min(costs(:,nextcolumn));
  bestpath=[bestpath,choice];
end

%since we anchor the last model at the first column beyond the line, we can
%safely prune this from our path and segmentation sequences
bestpath=bestpath(1:end-1);
bestseg=bestseg(1:end-1);
