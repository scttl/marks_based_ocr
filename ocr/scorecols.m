function cost=scorecols(data,logmp,logmpbar,insprob)
% cost=scorecols(data,logmp,logmpbar,insprob)

[hh,wdata]=size(data);
[hh2,wmodel]=size(logmp);

if(wdata>wmodel)
  cost=-[sum(data(:,1:wmodel).*logmp+(1-data(:,1:wmodel)).*logmpbar,1),...
	 sum(data(:,(wmodel+1):end)*log(insprob),1)];
elseif(wdata<wmodel)
  cost=-sum(data.*logmp(:,1:wdata)+(1-data).*logmpbar(:,1:wdata),1);
else
  cost=-sum(data.*logmp + (1-data).*logmpbar,1);
end

cost=cumsum(cost);
