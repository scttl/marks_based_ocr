function [logmp,logmpbar] = logmprob(model,delprob,insprob)
% [logmp,logmpbar] = logmprob(model,delprob,insprob)

mprob=double(model);
mprob(find(mprob(:)==0))=insprob;
mprob(find(mprob(:)==1))=1-delprob;

logmp=log(mprob);
logmpbar=log(1-mprob);
