% set up paths
restoredefaultpath
addpath('./cellobj')
addpath(genpath('./tools'))
% path to the BNEW code
addpath('../BNEW/git/BNEW')


%% play around with segmentation gui
tuneSegmentationParams('mnum',15)