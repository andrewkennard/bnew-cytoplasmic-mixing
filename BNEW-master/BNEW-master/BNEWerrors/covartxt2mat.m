function datamats = covartxt2mat(name, directory)
% convert a txt file containing covariance matrix information into a .mat
% file for easier load-in
% --------
% inputs:
% ---------
% name: name of .txt data file (without the .txt extension, without the
% path)
% directory: path to look for file; default is current path; must end in /
% ---------
% outputs: datamats structure containing the following fields
% ----------
% nmax: maximal wavelet span
% nkmax: maximal k value for each n
% nistart: linearized index for each n, k=1
% Hu: 4-th order drift covariance matrix
% Hv: 4-th order diffusion covariance matrix
% Hxi: 4-th order localization error covariance
% FuFv: coupling of drift and diffusion
% FuFxi: coupling of drift and localization error
% FvFxi: coupling of diffusion and localization error
% avals, bvals: rescaling functions for the wavelet
% fvals: bias component from drift
% ------------
% the datamats structure is saved in directory/name.mat


if (nargin<2)
    directory = './';
end

data = dlmread([directory, sprintf('%s.txt',name)]);

datamats.nmax = data(1,1);
datamats.nkmax = data(2,1:datamats.nmax);
datamats.nistart = cumsum([1,datamats.nkmax(1:end-1)]);

amax = sum(datamats.nkmax);
data = data(3:end,:);
datamats.Hu = data(1:amax,:);
datamats.Hv = data(amax+1:2*amax,:);
datamats.Hxi = data(2*amax+1:3*amax,:);
datamats.FuFv = data(3*amax+1:4*amax,:);
datamats.FuFxi = data(4*amax+1:5*amax,:);
datamats.FvFxi = data(5*amax+1:6*amax,:);
datamats.avals = data(6*amax+1,:);
datamats.bvals = data(6*amax+2,:);
datamats.fvals = data(6*amax+3,:);

save([directory sprintf('%s.mat',name)],'datamats')
