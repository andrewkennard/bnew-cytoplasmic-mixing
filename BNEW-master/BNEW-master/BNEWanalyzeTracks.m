function [MSD,MSDcnt,MSDerr,smvelsave] = BNEWanalyzeTracks(tracklist,nvals,ws, varargin)
% perform BNEW analysis on many tracks
% calculate MSD of adjusted tracks (with smoothed velocity subtracted out)
% this function is self-sufficient for analyzing a set of tracks,
% or can be called on a BNEWobj object
% ----------
% inputs:
% ----------
% tracklist: cell array of tracks. Each track is an Nx2 array, with N the
% number of time-points (assumed to be equispaced in time)
% nvals: list of spans n to use for the wavelet smoothing
% ws: cell array of coefficients defining the wavelet for each span n;
% (2n+1)x1 array of coefficients w_(-n), ..., w_0, ... w_n for each span
% ----------
% outputs
% ----------
% MSD = time and population avg of mean sq disp of corrected tracks
% MSDcnt = data counts associated with each point
% MSDerr = standard error associated with each point (not really accurate
% since assumes no correlation in time avg)
% smvel = cell array of smoothed velocity approximations for timepoints
% nmax+1 to (tracklength - nmax)

% ----------
% Optional parameters
% ----------
% how to do the time average
% 0: no time averaging at all (NOT SET UP)
% 1: fully overlapping intervals for time average (all intervals of size k)
% 2: intervals separated by k
% 3: intervals separated by k+2n
timeavg = 1;


for vc = 1:2:length(varargin)    
    switch (varargin{vc})                    
        case('timeavg')
            timeavg = varargin{vc+1};
    end 
end

%%

% set up all wavelet coeefficients in one large array
nmax = max(nvals);
allW = zeros(length(nvals),2*nmax+1);
for nc = 1:length(nvals)
    nn = nvals(nc);
    allW(nc,-nn+nmax+1:nn+nmax+1) = ws{nc}';
end

tracklens = cellfun(@(x) size(x,1),tracklist);
longind = find(tracklens>2*nmax+1);

noisetracklist = cell(length(tracklist),length(nvals));
smtracklist = noisetracklist;
smvelsave = noisetracklist;

for tc = longind
    % get noise tracks and smoothed tracks
    track = tracklist{tc};
    [trackcorr,tracksm,smvel] = applyWavelets(track,nvals,allW);
    for nc = 1:length(nvals)
        noisetracklist{tc,nc} = trackcorr{nc};
        smtracklist{tc,nc} = tracksm{nc};
        smvelsave{tc,nc} = smvel{nc};
    end
end    
    
MSD = cell(length(nvals),1);
MSDcnt = MSD;
MSDerr = MSD;

for nc = 1:length(nvals)
    nn = nvals(nc);
    
    %% MSD of noise tracks
    if (timeavg == 0)
        overlap = @(k) -1; % no time averaging
    elseif (timeavg == 1)
        overlap = @(k) 1;  % overlapping intervals
    elseif (timeavg(1)==2)
        overlap = @(k) k; % intervals separated by k
    elseif (timeavg(1)==3)
        overlap = @(k) 2*nn+k;
    else
        error('MSD timeaverage of this type not yet set up: %d', timeavg(1));
    end
    [MSD{nc},MSDcnt{nc},MSDerr{nc}] = MSDensemble(noisetracklist(:,nc),'overlap',overlap);
end

end