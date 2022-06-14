function [ws, what] = meanwavelet(n)
% get wavelet coefficients for sliding mean wavelet of span n

jvals = (-n:n);
ws = zeros(1,2*n+1);
ws(1)=-1/(2*n);
ws(end) = 1/(2*n);

% coefficients for the steps
what = zeros(1,2*n);
for j = -n:n-1
    what(j+n+1) = -sum(ws(1:j+n+1));
end


end