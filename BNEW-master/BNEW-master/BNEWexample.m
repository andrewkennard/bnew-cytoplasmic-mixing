% example run of BNEW analysis on a set of tracks

%% load pre-generated tracks
% diffusive motion
% simulated with alpha=1, D=5, epsilon = 1, gamma=1, tau=100
load('exampletracks.mat')
%
%%
% subdiffusive motion
% simulated with alpha=0.5, D=5, epsilon=1, gamma=1, tau=100
load('exampletracksFBM.mat')
tracklist = tracklistFBM;
load('fgfunc.mat')
%% alternately, simulate some tracks with diffusion, localization error, and persistent random drift
ntrack = 200; % number of tracks
tracklen = 200; %track length
tau = 100; % correlation time for drift
gam =1; % magnitude for drift
D = 5; %diffusion coefficient
locE = 1; % localization error
tracklist = simTracks(ntrack,tracklen,tau,gam,D,locE);

%% Run BNEW analysis for tracks with diffusive motion

% set of wavelet spans to use
nvals=4:4:20;

% make a BNEW analysis object
% using Haar wavelet type
% possible wavelet shapes are: 
% 'haar', 'svg' (Savitzky-Golay), and 'mean' (sliding-mean)
% the degree is only used for the Savitzky-Golay wavelet
BN = BNEWobj(nvals,'wavetype','svg','wavedeg',3);

% get wavelet coefficients 
BN = BN.getCoefficients();

% run wavelet analysis on 
BN = BN.analyzeTracks(tracklist);

% plot adjusted MSD
BN.plotMSD()

%%
% rescale to collapse everything to universal line, using A^n_k and B^n_k functions
BN = BN.rescaleData();

% fit coefficients (using 1<= k <= n)
[BN,parfit] = BN.fitDcoeff(struct('kmax',0.74));
% parameters are: D, epsilon, alpha
alphafit = parfit(3);

if (alphafit<1)
    % for subdiffusive motion get corrected diffusion and localization
    % error estimates
    fa = interp1(alphavals,fvalsSVG,alphafit); % assumes SVG wavelets
    ga = interp1(alphavals,gvalsSVG,alphafit);    
else
    fa = 1;
    ga = 0;
end

Dfit = parfit(1)/fa;
locEfit = sqrt(parfit(2)^2-ga^2);

% plot rescaled MSD and fitted power law
[tlist,fitvals] = BN.plotRescaled('kmax',0.74);

% write out fit parameters
tmin = min(tlist);
mmax = max(fitvals);
mmin = min(fitvals);

hold all
text(tmin,mmax-(mmax-mmin)*0.1,sprintf('D=%0.2f\n eps=%0.2f \n alpha=%0.2f',...
    Dfit,locEfit,alphafit))
hold off;

disp(sprintf('alpha fit: %f', alphafit))
disp(sprintf('D fit: %f', Dfit))
disp(sprintf('eps fit: %f', locEfit))