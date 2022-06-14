function [cellmv,jumps,alignopt] = cellAlignDriver(cellmv,options)
% driver routine for dealing with stage jumps by aligning images
% using approximate segmented contours to pick regions for
% cross-correlation alignment

% default parameters
opt = struct();
% maximal pixel displacement for COM to identify potential stage jump
% regions
opt.maxdisp=40;

% superset of all frames to check for stage jumps; by default, do all
opt.frames = 1:length(cellmv);

% clear all previous cell positions
opt.clearpos = 1;
% input options
opt = copyStruct(options,opt,'addnew',1);

%% make sure segmentation contours have been set
valframes = find(arrayfun(@(x) ~isempty(x.segcont),cellmv));
if (isempty(valframes))
    error('cannot align movie without segmentation contours')
end

%% link up individual cells based on COM
% tracks will break at stage jumps
if (opt.clearpos)
for cc = 1:length(cellmv)
    cellmv(cc).pos = [0 0];
end
end

display('Linking up cell COMs')
% link up individual cell tracks within the movie
[cellind, COMlist, cellmv] = linkCells(cellmv,'memory',1,'maxdisp',opt.maxdisp,'goodenough',2);

%%
% frames that will not be checked for stage jump (no stagejump preceeding this frame)
checkforjump = zeros(1,length(cellmv));
checkforjump(opt.frames) = 1;

checkforjump(1) = 0;
for tc = 1:length(COMlist)
    inds = COMlist{tc}(1,6)+1:COMlist{tc}(end,6);
    checkforjump(inds) = 0;
end

% also get rid of all frames before first track and after last
minframe = min(cellfun(@(x) x(1,6),COMlist));
maxframe = max(cellfun(@(x) x(end,6),COMlist));
if (minframe>1)
    checkforjump(1:minframe-1) = 0;
end
if (maxframe<length(cellmv)-1)
    checkforjump(maxframe+1:end) = 0;
end
disp('Will check following frames for stage jumps:')
find(checkforjump)

opt.checkforjump = checkforjump;

%% check for stage jumps by aligning images
if (nnz(checkforjump)>0)
    [jumps,cellmv,alignopt] = checkStageJumps(cellmv,opt);
else
    disp('No frames to be checked for stage jump')
    alignopt = opt;
    jumps = [];
end

%% copy over additional options used in this driver file
alignopt = copyStruct(opt,alignopt,'addnew',1);
end