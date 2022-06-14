function [bias,errvals,toterrvals,trscl,MSDrscl] = getFitParamErrors(nvals,gam,D,locE,ntrack,datamats1,datamats2,fu1u2mat,varargin)
% get the errors in the fitted parameters from BNEW analysis
% assuming drift composed of 2 independent persistent random walks
% (with different exponential velocity decorrelation times)
% ----------
% inputs:
% ----------
% nvals is list of n values used for fitting
% gam: magnitude of each drift (length 2 vector)
% D: diffusion coefficient
% locE: localization error in one dimension (epsilon in the paper)
% ntrack: number of independent tracks that go into the calculation
% datamats1: structure of data matrices calculated with fortran code
% obtain from output txt file using covartxt2mat.m
% these matrices are calculated for a particular correlation time of the drift velocity
% and for a particular length of each track
% datamats2: same thing for the second drift velocity
% fu1u2mat: coupling matrix between the drift velocities
% assumes the coupling matrix is calculated with the same k values as
% datamats1
% ----
% optional arguments (keyword, value pairs)
% kscl: for each span n, use time separations 1<= k <= floor(kscl*n) for
% the fitting; default kscl=1
% del: time step; default del=1; if changing del then the data matrices
% correspond to a different correlation time, since they were calculated
% using a preset value of tau/del
% fitD0: starting point to use for the fits
% ---------
% outputs:
% ---------
% bias: bias in fitted parameter (Dfit, locEfit,alphafit)
% errvals: sampling error in fitted parameters
% toterrvals: overall root mean square error in fitted params, including
% both bias and sampling error component
% trscl: rescaled time values used for fitting
% MSDrscl: rescaled MSD values used for fitting
% -------------------

% limit of k for each span n is floor(kscl*n)
kscl = 1;
% time step (for each drift velocity)
del = 1;
% starting point for fitting
fitD0 = [D,min(locE,0.1),1.5];

for vc=1:2:length(varargin)
    switch (varargin{vc})
        case('del')
            del = varargin{vc+1};
        case('fitD0')
            fitD0 = varargin{vc+1};
        case('kscl')
            kscl = varargin{vc+1};
    end
end


% indices of the data matrices to use
nkinds1 = []; nkinds2 = [];
for nc = 1:length(nvals)
    nn = nvals(nc);
    klim = floor(kscl*nn);
    inds1 = datamats1.nistart(nn):datamats1.nistart(nn)-1+klim;
    nkinds1 = [nkinds1, inds1];
    inds2 = datamats2.nistart(nn):datamats2.nistart(nn)-1+klim;
    nkinds2 = [nkinds2, inds2];
end


%%
fvals1 = datamats1.fvals(nkinds1);
avals = datamats1.avals(nkinds1);
bvals = datamats1.bvals(nkinds1);

Hu1 = datamats1.Hu(nkinds1,nkinds1);
Hv = datamats1.Hv(nkinds1,nkinds1);
Hxi = datamats1.Hxi(nkinds1,nkinds1);
FuFv1 = datamats1.FuFv(nkinds1,nkinds1);
FuFxi1 = datamats1.FuFxi(nkinds1,nkinds1);
FvFxi = datamats1.FvFxi(nkinds1,nkinds1);

Hu2 = datamats2.Hu(nkinds2,nkinds2);
FuFv2 = datamats2.FuFv(nkinds2,nkinds2);
FuFxi2 = datamats2.FuFxi(nkinds2,nkinds2);
fvals2 = datamats2.fvals(nkinds2);

% full covariance matrix containing all components for drift, diffusion,
% localization error
Mcov1 = (gam(1)*del)^4*Hu1 + (2*D*del(1))^2*Hv + locE^4*Hxi + 4*D*del*(gam(1)*del)^2*FuFv1 ...
    + 2*(gam(1)*del)^2*locE^2*FuFxi1 + 4*D*del*locE^2*FvFxi;
Mcov2 = (gam(2)*del)^4*Hu2  + 4*D*del*(gam(2)*del)^2*FuFv2 + 2*(gam(2)*del)^2*locE^2*FuFxi2 ;

Mcov = Mcov1 + Mcov2 + 2*(gam(1)*del)^2*(gam(2)*del)^2*fu1u2mat(nkinds1,nkinds1);

MSD = del^2*(gam(1)^2*fvals1+gam(2)^2*fvals2)+4*D*del*avals+4*locE^2*bvals;

%% rescaled values
MSDrscl = MSD./(bvals);
trscl = avals./(bvals);
Mcov = Mcov./(bvals'*bvals);

fitfunc = @(D,t) (4*D(1)*(del*t).^D(3)+4*D(2).^2);

Dfit= nlinfit(trscl,MSDrscl,fitfunc,fitD0);
Dfit(2) = abs(Dfit(2));

% matrix of derivatives of fit function wrt parameters
Z = zeros(length(trscl),3);
Z(:,1) = 4*(del*trscl).^(Dfit(3));
Z(:,2) = 8*Dfit(2);
Z(:,3) = 4*Dfit(1)*del*log(del*trscl).*(del*trscl).^(Dfit(3));
invZZ = inv(Z'*Z);

errmat = invZZ*Z'*Mcov*Z*invZZ/ntrack;
toterr = sqrt(diag(errmat)'+(Dfit-[D,locE,1]).^2);
errvals = sqrt(diag(errmat)');
bias = Dfit-[D,locE,1];
toterrvals = toterr;
end