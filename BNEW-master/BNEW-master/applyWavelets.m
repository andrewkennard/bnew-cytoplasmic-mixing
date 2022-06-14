function [trackcorr,tracksm,smvelsave] = applyWavelets(track,nvals,allW)
% for a given track (Nx2 array), apply the wavelets defined by (position) weights ws
% returns the smoothed track and the leftover noise track
% nn is span of wavelet

% set up the wavelet coefficients in one large array
nmax = max(nvals);

%% apply wavelets to whole track 
nt = size(track,1);
smvelsave = cell(length(nvals),1);
for dim = 1:2
    smvels = zeros(length(nvals),nt);
    for j = -nmax:nmax
        smvels(:,nmax+1:nt-nmax) = smvels(:,nmax+1:nt-nmax)  + allW(:,j+nmax+1)*track(nmax+1+j:nt-nmax+j,dim)';
    end
    
    for nc = 1:length(nvals)
        smvelsave{nc}(:,dim) = smvels(nc,nmax+1:nt-nmax)';
    end
    
end
%%
tracksm = cell(length(nvals),1); 
trackcorr = tracksm;

for nc = 1:length(nvals)
    % cumulative track
    tracksm{nc} = cumsum([track(nmax+1,1:2);smvelsave{nc}(1:end-1,1:2)]);
    
    %% corrected track
    trackcorr{nc} = track(nmax+1:nt-nmax,1:2) - tracksm{nc};
end
