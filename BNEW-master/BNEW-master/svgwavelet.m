function [ws, what] = svgwavelet(n,deg)
% get wavelet coefficients for savitsky-golay filter of span n, degree deg
% wavelet coefficients are for extracting a first derivative (velocity)

jvals = (-n:n);%-0.5;

A = ones(2*n+1,deg+1);
for i = 1:deg
    A(:,i+1) = jvals.^i; 
end

G = A'*A;
g = zeros(deg+1,1);
g(2) = 1;

lam = G\g;

ws = A*lam;


% coefficients for the steps
what = zeros(1,2*n);
for j = -n:n-1
    what(j+n+1) = -sum(ws(1:j+n+1));
end

ws = ws'; 

end