function [MSDscl,tscl,allMSD,alltscl] = fBMmsdRescale(nvals,alpha,D,locE,options)
% do analytical calculations for rescaled MSD after BNEW analysis
% of trajectories involving constant drift, localization error, and 
% fractional brownian motion
% ----------
% inputs
% ------------
% nvals: list of wavelet spans to use for the smoothing
% alpha: scaling coefficient for the fractional Brownian motion
% D: effective diffusion coefficient, in units of length^2 / timestep
% locE: localization error (units of length)
%
% options: structure with possible fields....
%
% wavetype: wavelet shape to use (options are svg [Savitzky-Golay], haar, mean)
% wavedeg: degree of wavelet for svg wavelets
% default is 3rd order svg wavelet
% gam: magnitude of drift velocity (length / timestep)
% tau: correlation time for drift velocity (in timesteps)
% kmax: save k values 1<=k<= floor(kmax*n) in the combined array allMSD for fitting
% kmaxsave: use k values 1 <= k <= floor(kmaxsave*n) for saving individual
% MSD arrays (for plotting)
% array allMSD
%
% --------------
% outputs
% --------------
% MSDscl: rescaled mean square displacement (cell array giving the rescaled
% MSD for each wavelet span)
% tscl: rescaled time (cell array)
% allMSD: all the rescaled MSD from all spans compiled together
% alltscl: all the rescaled times from all spans compiled together


% default options
opt = struct();

% scale k as a function of n, for saving in allMSD, to be used in fitting
opt.kmax = 1;
% scale k for saving in individual MSDrscl arrays
opt.kmaxsave = 1;

% wavelet type and degree
opt.wavetype = 'svg';
opt.wavedeg = 3;

opt.gam = 0;
opt.tau = 100;

% input options
if (nargin>1)
    opt = copyStruct(options,opt);
end
            
            
% inter-velocity correlation function
vvcorr = @(k) 2*D*((k+1)^alpha - 2*k^alpha+(k-1)^alpha);
alltscl=[]; allMSD=[];
for nc=1:length(nvals)
    n = nvals(nc);
    
    switch(opt.wavetype)
        case('svg')
            [ws, what] = svgwavelet(nvals(nc),opt.wavedeg);
        case('haar')
            [ws,what] = haarwavelet(nvals(nc));
        case('mean')
            [ws,what] = meanwavelet(nvals(nc));
        otherwise
            error('Wavetype %s is invalid. Must be svg, haar, or mean')
    end
           
    klist = 1:floor(n*opt.kmaxsave);
    A = zeros(1,length(klist));
    B = zeros(1,length(klist));
    MSDfbm = A; 
    
    for kc = 1:length(klist)
        k = klist(kc);
        [cts,ctshat] = getCvals(what,n,k);
        A(kc) = sum(cts.^2);
        B(kc) = 0.5*sum(ctshat.^2);
        
        MSDfbm(kc) = 4*D*A(kc)+4*locE^2*B(kc);
        
        % get the fbm correlation terms
        % factor of 2 b/c two dimensional trajectories
        for ic = -n:k+n-2
            for jc = ic+1:k+n-2
                MSDfbm(kc) = MSDfbm(kc) + 2*cts(ic+n+1)*cts(jc+n+1)*vvcorr(jc-ic);
            end
        end
        
        if (opt.gam~=0)
            % get the drift correlation terms            
            for ic = -n:k+n-2
                for jc = -n:k+n-2
                    MSDfbm(kc) = MSDfbm(kc) + 2*cts(ic+n+1)*cts(jc+n+1)*opt.gam^2*exp(-abs(jc-ic)/opt.tau);
                end
            end
        end
    end
    
    % rescale
    MSDscl{nc} = MSDfbm./B;
    tscl{nc} = A./B;
    
    timax = floor(n*opt.kmax);
    
    alltscl = [alltscl, tscl{nc}(1:timax)];
    allMSD = [allMSD,MSDscl{nc}(1:timax)];    
end

end