function [thresh,vv,smvels,meanvels,totmeanvel] = getThreshold(tracklist,nsm,varargin)
% determine threshold level based on a combination of universal threshold
% to remove noise + one standard deviation of smoothed velocities
% nsm is the wavelet half-span to use for smoothing
% vv is the mean smoothed velocity in each dimension
% smvels are the smoothed velocities used to determine threshold

% how to weigh the standard deviations from flow variation vs noise
flowweight = 1;
noiseweight=1;
for vc = 1:2:length(varargin)
    switch(varargin{vc})
        case('flowweight')
            flowweight = varargin{vc+1};
        case('noiseweight')
            noiseweight = varargin{vc+1};
    end
end

%%
tracklens = cellfun(@(x) size(x,1),tracklist);

nvals = nsm;
BN = BNEWobj(nvals);

% get wavelet coefficients
BN = BN.getCoefficients();

% run wavelet analysis on
BN = BN.analyzeTracks(tracklist);

BN20 = BN;

%% smoothed velocities
smvels = BN.SmVel(:,end);
allvel = []; meanvel = []; allnvel = [];
for tc = 1:length(smvels)
    allvel = [allvel; smvels{tc}];    
end
vv=mean(allvel);
flowthresh = sqrt(sum(var(allvel)));

%%
if (nargout>3)
    meanvels = zeros(1,length(smvels));
    totmeanvel = 0; ct=0;
    for tc = 1:length(smvels)
        tmp = bsxfun(@minus,smvels{tc},vv);
        nvel = tmp(:,1).^2+tmp(:,2).^2;
        meanvels(tc) = mean(nvel);
        totmeanvel = totmeanvel + sum(smvels{tc}(:,1).^2 + smvels{tc}(:,2).^2);        
        ct = ct + length(nvel);
    end
    totmeanvel = sqrt(totmeanvel/ct);
end
%% noise estimate
nvals = 1;
BN1 = BNEWobj(nvals,'wavetype','haar');

% get wavelet coefficients
BN1 = BN1.getCoefficients();

% run wavelet analysis on
BN1 = BN1.analyzeTracks(tracklist);

%%
smvels1 = BN1.SmVel(:,end);
allvel = []; allnvel = [];
for tc = 1:length(smvels1)
    nvel = sqrt(smvels1{tc}(:,1).^2+smvels1{tc}(:,2).^2);
    allnvel = [allnvel;nvel];
    allvel = [allvel; smvels1{tc}];    
end
tmp = bsxfun(@minus,allvel, vv);
noiseest = median(sqrt(tmp(:,1).^2+tmp(:,2).^2))/0.6745;

% final threshold
thresh = sqrt(flowweight*flowthresh.^2+noiseweight*(sqrt(2*log(tracklens-2*nsm))*noiseest/sqrt(nsm)).^2);

end