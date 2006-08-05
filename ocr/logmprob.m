function [logmp,logmpbar] = logmprob(model,delprob,insprob)
% LOGMPROB    Determine log insertion and deletion probs. for a character model
%
%  [logmp,logmpbar] = logmprob(model,delprob,insprob)
%
%  model should be a binary image representing a single character or symbol.
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
%  logmp is the resulting matrix giving the log score for each pixel based on
%  the insprob and delprob value passed; essentially computing a false 
%  positive score (having additional foreground pixels that aren't present in
%  the data). logmpbar uses the values to compute a false negative score
%  (missing foreground pixels that are present in the data).


% CVS INFO %
%%%%%%%%%%%%
% $Id: logmprob.m,v 1.2 2006-08-05 17:31:01 scottl Exp $
%
% REVISION HISTORY
% $Log: logmprob.m,v $
% Revision 1.2  2006-08-05 17:31:01  scottl
% added comment header.
%

% LOCAL VARS %
%%%%%%%%%%%%%%


% CODE START %
%%%%%%%%%%%%%%
mprob=double(model);
mprob(find(mprob(:)==0))=insprob;
mprob(find(mprob(:)==1))=1-delprob;

logmp=log(mprob);
logmpbar=log(1-mprob);
