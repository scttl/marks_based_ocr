function [bestpath,bestseg]=solveline(data,models,bigram,delprob,insprob,Wmin,Wmax)
%[bestpath,bestseg]=solveline(data,models,bigram,delprob,insprob,Wmin,Wmax)


more off;

[hh,N]=size(data);
K=length(models);
maxcharwidth=0;
for kk=1:K
  maxcharwidth=max(maxcharwidth,size(models{kk},2));
  [logmp{kk},logmpbar{kk}]=logmprob(models{kk},delprob,insprob);
end

costs=NaN(K,N+maxcharwidth+Wmax); costs(:,(N+2):end)=Inf; costs(K,N+1)=0;
bestpaths=zeros(K,N);
bestdeltas=zeros(K,N);

data=[data,zeros(hh,maxcharwidth+Wmax)];

bgcost=-log(bigram);


thispos=N;
while(thispos>=1)
  fprintf('Computing costs for column %d          \r',thispos);
  newcost=zeros(K,1);
  for cc=1:K
    thisww=size(logmp{cc},2);
    sc=scorecols(data(:,thispos:(thispos+thisww+Wmax-1)),...
		logmp{cc},logmpbar{cc},insprob);
    tmpcost=costs(:,(thispos+Wmin):(thispos+thisww+Wmax))+...
	    bgcost(cc,:)'*ones(1,thisww+Wmax-Wmin+1)+ones(K,1)*sc(Wmin:end);
    [newcost(cc),minidx]=min(tmpcost(:));
    [bestkk,bestrr]=ind2sub([K,thisww+Wmax-Wmin+1],minidx(1)); 
                                                  % ARBITRARILY TAKE FIRST!
    bestpaths(cc,thispos)=bestkk;
    bestdeltas(cc,thispos)=bestrr+Wmin-1;
  end
  costs(:,thispos)=newcost;
  thispos=thispos-1;
end

[tmpcost,bestpath]=min(costs(:,1));
bestseg=1;
while(bestseg(end)<N)
  nextdelta=bestdeltas(bestpath(end),bestseg(end));
  nextcolumn=bestseg(end)+nextdelta;
  bestseg=[bestseg,nextcolumn];
  [tmpcost,choice]=min(costs(:,nextcolumn));
  bestpath=[bestpath,choice];
end

