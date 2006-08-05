function cost=scorecols(data,insdist,deldist,logmp,logmpbar,insprob)
% SCORECOLS    Determine the score of placing a character over the image passed.
%
%  cost = scorecols(data, insdist, deldist, logmp, logmpbar, insprob)
%
%  data should be an image array representing one portion of a line of text
%  that we will score the model against.
%
%  insdist and deldist should be the same dimensions as the data, and provide a
%  weighting to penalize false positives and negatives based on their distance
%  from the closest match.
%
%  logmp and logmpbar are log scores of insertion and deletion probabilities
%  (essentially false positive and false negative scores) for placing a
%  particular character symbol down.  These should be calculated via a call to
%  logmprob.m
%
%  insprob is the score assigned to a mismatch whereby a foreground pixel must
%  be added to the model to match a corresponding foreground pixel in the data.
%  See logmprob.m for further details.  It should lie in the range (0...1).
%
%  the cost returned is a vector denoting the cumulative sum of laying down
%  each column of the model overtop of the data.


% CVS INFO %
%%%%%%%%%%%%
% $Id: scorecols.m,v 1.2 2006-08-05 17:34:22 scottl Exp $
%
% REVISION HISTORY
% $Log: scorecols.m,v $
% Revision 1.2  2006-08-05 17:34:22  scottl
% added comment header.  Weighted penalties, and took the cost of placing
% the rest of the character down into account.
%


% LOCAL VARS %
%%%%%%%%%%%%%%

% CODE START %
%%%%%%%%%%%%%%
if nargin ~= 6
    error('incorrect number of arguments specified');
end

[wdata,wdata]=size(data);
[wmodel,wmodel]=size(logmp);

if(wdata>wmodel)
  falsepos_cost = (1-data(:,1:wmodel)).*deldist(:,1:wmodel).*logmpbar;
  cost=-[sum(data(:,1:wmodel).*insdist(:,1:wmodel).*logmp+falsepos_cost,1),...
         sum(data(:,(wmodel+1):end).*insdist(:,(wmodel+1):end)*log(insprob)+...
         (1-data(:,(wmodel+1):end)).*deldist(:,(wmodel+1):end)*...
         log(1-insprob),1)];
elseif(wdata<wmodel)
    %@@@ NOTE: this case is not possible based on how scorecols is currently
    %called.  If we do reach here, this code doesn't currently take into account
    %additional 'on' pixels in the model that aren't present in the data due to
    %the difference in length.
  cost=-sum(data.*logmp(:,1:wdata)+(1-data).*logmpbar(:,1:wdata),1);
else
  falsepos_cost = (1-data).*deldist.*logmpbar;
  cost=-sum(data.*insdist.*logmp + falsepos_cost,1);
end

cost=cumsum(cost);

%also need to penalize 'on' pixels in the model that don't match the data in
%those columns beyond those being considered.
model_completion_cost=zeros(1,wdata);
model_completion_cost(1:wmodel) = -sum(falsepos_cost,1);
model_completion_cost=cumsum(model_completion_cost(wdata:-1:1));

cost = cost + model_completion_cost;
