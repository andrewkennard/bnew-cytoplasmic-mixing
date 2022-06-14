function [track,velconv1,velconv2] = genTrackFBMdrift(nsamp,FBMpow,FBMweight,driftstd,per,locE,tstep,driftmean,driftmean2)
% generate a simulated track consisting of FBM overlayed on periodic drift
% (in 2D)
% nsamp = number of points to sample
% FBMpow, FBMweight = exponent and prefactor for FBM such that
% MSD~FBMweight*t^FBMpow
% driftmag = standard dev of drift magnitudes (for constant and oscillating
% component)
% per = period (in steps) for the drift oscillation
% locE = localization error
% tstep = time difference associated with each step

tlist = (1:nsamp)*tstep;

H = FBMpow/2;

FBM1 = genfbm1d(H,nsamp);
FBM1 = FBM1'*sqrt(FBMweight)*(tstep*nsamp)^H;
FBM2 = genfbm1d(H,nsamp);
FBM2 = FBM2'*sqrt(FBMweight)*(tstep*nsamp)^H;

velFBM1 = diff(FBM1)/tstep;
velFBM2 = diff(FBM2)/tstep;

if (nargin<8)
    driftmean = [0 0 0];
end
if (nargin<9)
    driftmean2 = driftmean;
end

mags1 = randn(1,3).*driftstd+driftmean;
%mags1 = 0.5*driftmag;
velconv1 = mags1(1)+mags1(2)*cos(2*pi*tlist/per);
%velconv1 = mags1(1)+mags1(2)*tlist +mags1(3)*tlist.^2;
mags2 = randn(1,3).*driftstd+driftmean2;
%mags2 = 0.3*driftmag;
velconv2 = mags2(1)+mags2(2)*cos(2*pi*tlist/per);
%velconv2 = mags2(1)+mags2(2)*tlist+mags2(3)*tlist.^2;

%[mags1;mags2]

veltot1 = velFBM1+velconv1;
veltot2 = velFBM2 + velconv2;

error1 = randn(1,nsamp+1)*locE;
error2 = randn(1,nsamp+1)*locE;

possim1 = cumsum([0 veltot1])+error1;
possim2 = cumsum([0 veltot2])+error2;
possim1conv = cumsum(velconv1);
possim2conv = cumsum(velconv2);

%plot(tlist(1:per*2),possim1(1:per*2),'.-')

track = [possim1(1:end-1); possim2(1:end-1)]';
end