% load in all cell objects along with their lysotracker data

% ------------------------
% Uncomment these lines to load in preprocessed cell objects and
% trajectories
% --------------------------
% will load files with names lyso*
% The second list includes cells that were judged bad for one reason or another and will be excluded

%celldir = '../Cells_HL60_intpix/';
%Cells = loadCellObj(celldir,{'lyso*'},{'lysozplane10_z0_2','lysozplane7_z0_2','lysozplane9_z0_2'});

% -----------------
% Uncomment these lines to load in results of example processing
% ------------------
celldir = './example/results/';
Cells = loadCellObj(celldir,{'example15_*'},{});


AllCells= Cells;


%% keep only cells with at least 50 tracks of minimal length
ntracks = zeros(length(AllCells),1);
alldt = ntracks;
for ccell = 1:length(AllCells)
    CL = AllCells(ccell);
    tracklist = CL.getTracklist('splitskip',1,'mintracklen',42);
    ntracks(ccell) = length(tracklist);
    alldt(ccell) = mean(diff([CL.Frames.time]));
end
goodind = find(ntracks>=50 & alldt<0.051);
Cells = AllCells(goodind);
alldt = alldt(goodind);

% ----------------------------------
%% Particle trajectories from example cell 

% which cell are we interested in?
ccell = 3
CL = Cells(ccell);

% interpolate the cell frame of reference to frames between the contours
CL.interpCellFrameRef();
%
cmat = jet(CL.NFrame);
% get tracks relative to the cell frame of reference
% can set flag to 0 to get lab frame of reference
tracklist = CL.getTracklist('relCellPos',0);
for tc = 1:length(tracklist)
    track = tracklist{tc};
    plot(track(:,1),track(:,2),'Color',cmat(track(1,6),:))
    hold all
end
title(sprintf('Cell %d: %s', ccell,CL.Name))
set(gca,'Visible','off')
hold off
axis equal

%% BEWARE: 
% scale bar 
um_per_px = 0.06435006; % um per pixel conversion

hold all
% draw a 10-um scale bar
plot([-350,-350+10/um_per_px],[200,200],'k','LineWidth',2)
%plot([600,600+20/(0.086*1.5)],[800,800],'k','LineWidth',2)
hold off

