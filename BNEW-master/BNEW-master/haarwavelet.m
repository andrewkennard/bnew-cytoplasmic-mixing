function [ws, what] = haarwavelet(nn)
% get the haar wavelet for half-span nn
% ws are coefficients for points, what are coefficients for velocities

ws = zeros(1,2*nn+1);
ws(1:nn) = -1;
ws(nn+2:end) = 1;
ws = ws/(nn*(nn+1));


% coefficients for the steps
what = zeros(1,2*nn);
for j = -nn:nn-1
    what(j+nn+1) = -sum(ws(1:j+nn+1));
end

end