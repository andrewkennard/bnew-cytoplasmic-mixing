function [bias,errvals,toterrvals,trscl,MSDrscl] = getFitParamErrors(nvals,gam,D,locE,ntrack,datamats,varargin)
% get the errors in the fitted parameters from BNEW analysis
% assuming single persistent random walk drift velocity
% ----------
% inputs:
% ----------
% nvals is list of n values used for fitting
% gam: magnitude of drift
% D: diffusion coefficient
% locE: localization error in one dimension (epsilon in the paper)
% ntrack: number of independent tracks that go into the calculation
% datamats: structure of data matrices calculated with fortran code
% obtain from output txt file using covartxt2mat.m
% these matrices are calculated for a particular correlation time of the drift velocity
% and for a particular length of each track
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
% time step
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
            klimscl = varargin{vc+1};
    end
end


% indices of the data matrices to use
nkinds = [];
for nc = 1:length(nvals)
    nn = nvals(nc);
    klim = floor(kscl*nn);
    inds = datamats.nistart(nn):datamats.nistart(nn)-1+klim;
    nkinds = [nkinds, inds];
end


%%
fvals = datamats.fvals(nkinds);
avals = datamats.avals(nkinds);
bvals = datamats.bvals(nkinds);

Hu = datamats.Hu(nkinds,nkinds);
Hv = datamats.Hv(nkinds,nkinds);
Hxi = datamats.Hxi(nkinds,nkinds);
FuFv = datamats.FuFv(nkinds,nkinds);
FuFxi = datamats.FuFxi(nkinds,nkinds);
FvFxi = datamats.FvFxi(nkinds,nkinds);

% full covariance matrix containing all components for drift, diffusion,
% localization error
Mcov = (gam*del)^4*Hu + (2*D*del)^2*Hv + locE^4*Hxi + 4*D*del*(gam*del)^2*FuFv ...
    + 2*(gam*del)^2*locE^2*FuFxi + 4*D*del*locE^2*FvFxi;

MSD = (gam*del)^2*fvals+4*D*del*avals+4*locE^2*bvals;

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