function [that,MSDhat,Ank,Bnk] = rescaleMSD(tvals,MSDtot,what,nn)
% given times and corrected mean square displacements, rescale values in
% accordance with wavelet defined by what (step coefficients)
% so that it should collaps to universal line (assuming diffusive Brownian
% motion and constant drift)
% nn is wavelet span
% outputs:
% that = rescaled time
% MSDhat= rescaled MSD
% Ank, Bnk = rescaling functions

Ank = zeros(1,length(tvals));
Bnk = Ank;
for kc = 1:length(tvals)
    k = tvals(kc);
    [cts,ctshat] = getCvals(what,nn,k);
    Ank(kc) = sum(cts.^2);
    Bnk(kc) = 0.5*sum(ctshat.^2);
end

that = Ank./Bnk;
MSDhat = MSDtot./Bnk;

end