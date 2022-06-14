%% generate a tabulation of f(alpha), g(alpha) for a range of alpha values
% to use in extracting diffusion coefficient and localization error in the
% case of fractional Brownian motion

alphalist = linspace(0.2,1,50);

nvals = 3:13;
optionsHaar = struct('kmax',1,'kmaxsave',1,'wavetype','haar');
optionsSVG = struct('kmax',0.74,'kmaxsave',0.74,'wavetype','svg');

D=0.25;
locE=0;

cmat = lines(3);
clear allDfitSVG allDfitHaar
for ac = 1:length(alphalist)
    alpha = alphalist(ac);
            
    %get rescaled MSD, Haar wavelets
    [MSDscl,tscl,allMSD,alltscl] = fBMmsdRescale(nvals,alpha,D,locE,optionsHaar);
        
    % fit power law   
    fitfunc = @(D,t) D(1)*(t).^D(3) + D(2)^2;    
    [Dfit,R,J,CovFit,MSE,ErrorModelInfo] = nlinfit(alltscl,allMSD,fitfunc,[1,0.1,alpha]);
    allDfitHaar(ac,:) = Dfit;
    
    % SVG wavelets     
    [MSDscl,tscl,allMSD,alltscl] = fBMmsdRescale(nvals,alpha,D,locE,optionsSVG);
    
    % fit power law        
    fitfunc = @(D,t) D(1)*(t).^D(3) + D(2)^2;    
    [Dfit,R,J,CovFit,MSE,ErrorModelInfo] = nlinfit(alltscl,allMSD,fitfunc,[1,0.1,alpha]);
    allDfitSVG(ac,:) = Dfit
end
%%
alphavals = alphalist;
fvalsSVG = allDfitSVG(:,1);
fvalsHaar = allDfitHaar(:,1);
gvalsSVG = allDfitSVG(:,2);
gvalsHaar = allDfitHaar(:,2);

%%
save('fgfunc.mat','alphavals','gvalsHaar','gvalsSVG','fvalsSVG','fvalsHaar','nvals','optionsHaar','optionsSVG')

%% plot f(alpha) and g(alpha)
plot(alphavals,fvalsSVG,'b--')
hold all
plot(alphavals,fvalsHaar,'b-')
plot(alphavals,gvalsSVG,'r--')
plot(alphavals,gvalsHaar,'r-')
hold off
axis square