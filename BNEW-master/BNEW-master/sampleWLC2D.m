function beads = sampleWLC2D(lp,gam,npt,del)
% sample wormlike chain (or persistent random walk) in two dimensions
% lp2: persistence length (time)
% gam: instantaneous velocity magnitude 
% npt: number of points to sample
% del: length (time) interval
% energy function for step is exp(-theta^2*del/2/lp)

%% sample turning angles
dtvals = randn(1,npt-2)*sqrt(del/lp);

thetavals = zeros(1,npt-1);

% uniformly sample original orientation
thetavals(1) = rand()*2*pi;

thetavals(2:end) = thetavals(1)+cumsum(dtvals);

% actual bead positions
beads = zeros(2,npt);

for pc = 1:npt-1
    beads(:,pc+1) = beads(:,pc)+gam*del*[cos(thetavals(pc));sin(thetavals(pc))];
end
end