function [tracklist,tracksmooth] = simTracks(ntrack,tracklen,tau,gam,D,locE)
% simulate many two-dimensional tracks consisting of diffusion, localization error, and
% persistent random walk drift
% ---------
% inputs:
% --------
% ntrack = number of tracks
% tracklen = length of each track
% tau = decorrelation time of drift velocities
% gam = magnitude of drift velocities
% D = diffusion coefficient
% locE = localization error (epsilon), in each dimension
% -----------
% outputs:
% ----------
% tracklist: cell array of tracks of size tracklen x 2
% tracksmooth: cell array of corresponding tracks with drift only (without
% diffusion and localization error)
% -----------


tracklist = {}; 
tracksmooth = {};
% persistence length for wormlike chain (persistent random walk)
% divide by 2 so that velocity correlations scale as exp(-t/tau)
lt = 0.5*tau;
dels = 1;
npt = tracklen;


for tc = 1:ntrack
    
    % WLC drift
    beads = sampleWLC2D(lt,gam,npt,dels);    
    driftvel = diff(beads');
   
    % diffusive velocities
    diffvel = randn(tracklen-1,2)*sqrt(2*D*dels);
    
    % localization errors
    locerr = randn(tracklen,2)*locE;
    
    tracklist{tc} = zeros(tracklen,8);
    tracklist{tc}(:,1:2) = [0 0;cumsum(driftvel+diffvel)] + locerr;   
    tracksmooth{tc}(:,1:2) = [0 0;cumsum(driftvel)];
end