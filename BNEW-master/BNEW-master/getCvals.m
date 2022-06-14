function [cts,ctshat] = getCvals(what,n,k)
% get coefficients (cts) for each step s_i (-n to n+k-2)
% based on wavelet vector what, with span n, separation k
% coefficients in final expression for adjusted x_k - x_0
% ctshat gives the coefficients for each point p

cts = zeros(1,2*n+k);

for j = -n:n+k-1
    ind = max(j-k+1,-n):min(j,n-1);
    ind = ind+n+1;
    cts(j+n+1) = sum(what(ind));
end

% plain velocities minus corrected ones
ind = (0:k-1)+n+1;
cts(ind) = cts(ind)-1;

cts = -cts;

ctshat = zeros(1,2*n+k+1);
ctshat(1) = -cts(1);
ctshat(2:end-1) = cts(1:end-1)-cts(2:end);
ctshat(end) = cts(end);
end