function [cellmv,segopt] = runSegmentation(cellmv,savefile,options)
% driver for segmenting cell movie

%% default options
opt=struct();
% cutoff for brightness jump indicating phase
% if negative, don't worry about brightness levels at all
opt.brjumpcutoff = 10;
% width of phase stretch (in frames)
opt.phasefr = 6;
% clear all segcont fields before starting
opt.clearsegcont=0;
% frames to use for segmentation (phase frames will be excluded)
% default = every 10 frames + first and last
opt.segframelist = 10:10:length(cellmv);
if (opt.segframelist(1)>1)
    opt.segframelist=[1 opt.segframelist];
end
if (opt.segframelist(end)<length(cellmv))
    opt.segframelist=[opt.segframelist length(cellmv)];
end

% input options
% only copy over the ones appropriate for this function
inputopt = fieldnames(options);
for c = 1:length(inputopt)
    s = inputopt(c); s=s{1};
    if (isfield(opt,s))
        opt.(s) = options.(s);
    end
end

%% pick out phase frames
if (opt.brjumpcutoff>0)
diffs = diff([cellmv.meanbr]);
brjumps = find(diffs>opt.brjumpcutoff);
%plot(1:length(cellmv),[cellmv.meanbr],brjumps,130*ones(size(brjumps)),'r.')
phaseind = [];
for bc = 1:length(brjumps)
    phaseind = [phaseind brjumps(bc):brjumps(bc)+opt.phasefr]
end

% get rid of frames that have phase on
opt.segframelist = opt.segframelist(~ismember(opt.segframelist,phaseind));
end
%%
if (opt.clearsegcont && isfield(cellmv,'segcont'))
    disp('Clearing old segmentation contours')
    cellmv = rmfield(cellmv,'segcont')
    cellmv(1).segcont=[];
end

%% perform segmentation

[cellmv, segopt] =segmentLysotracker(cellmv,opt.segframelist,savefile,options)

% add extra fields to option structure that are relevant for this driver
% only
inputopt = fieldnames(opt);
for c = 1:length(inputopt)
    s = inputopt(c); s=s{1};
    if (~isfield(segopt,s))
        segopt.(s) = opt.(s);
    end
end


end