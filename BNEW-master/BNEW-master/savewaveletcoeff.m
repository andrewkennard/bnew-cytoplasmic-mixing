% calculate relevant coefficients for different wavelet shapes and spans
% and save to file

%% save wavelet information for svg3 wavelets
klist = 1:3000;
nvals = 2:1:100;

% type of wavelet
wavetype='svg';
% polynomial degree (for velocity smoothing)
wavedeg=3;

Afunc = zeros(length(nvals),length(klist));
Bfunc = Afunc;
ws = {}; what = {};

for nc = 1:length(nvals)     
    nn = nvals(nc)    
    [ws{nn}, what{nn}] = svgwavelet(nn,wavedeg);
    for kc = 1:length(klist)
        kk = klist(kc);
        [cts,ctshat] = getCvals(what{nn},nn,kk);
        Afunc(nn,kk) = sum(cts.^2);
        Bfunc(nn,kk) = 0.5*sum(ctshat.^2);
    end
end

save('svg3waveletcoeff.mat','wavetype','wavedeg','Afunc','Bfunc','ws','what')
