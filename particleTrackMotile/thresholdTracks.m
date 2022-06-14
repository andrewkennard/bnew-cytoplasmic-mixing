function [belowthresh,abovethresh] = thresholdTracks(tracklist,smvels,indshift,thresh,varargin)
% threshold tracks using smoothed velocities
% velocity indices are shifted from track indices by indshift (=n for BNEW
% smoothed velocities)
% thresh is the thresholding level for smoothed velocity magnitudes
% returns the portion
% thresh can be single value or a list of values for each track

belowthresh = tracklist;
abovethresh = tracklist;

meanshift = [0 0];
% minimal length above threshold necessary to count (in order to break up
% into below-threshold regions)
minabovelen = 1;

for vc = 1:2:length(varargin)
    switch(varargin{vc})
        case('meanshift')
            meanshift = varargin{vc+1};
        case('minabovelen')
            minabovelen = varargin{vc+1};
    end
end

if (length(thresh)==1)
    thresh = thresh*ones(1,length(tracklist));
end

abovethresh = {};
for tc = 1:length(tracklist)
    %%
    track = tracklist{tc};
    vel = smvels{tc};
    vel = bsxfun(@minus, vel,meanshift);
    % corresponding indices of the track points
    trackind = (indshift+1:indshift+size(vel,1))';
    nvel = sqrt(vel(:,1).^2+vel(:,2).^2);
    % indices in velocity array that are above threshold
    aboveind = find(nvel>thresh(tc));
    
    isbelow = ones(size(track,1),1);
    %%
    if (isempty(aboveind))
        abovetracklist = {};
        keepind = [];
    else
        %% split into disjoint regions above the threshold
        abovetrack = track(trackind(aboveind),:);
        abovetrack(:,3) = trackind(aboveind)';
        
        abovetracklist = splitTrackSkip({abovetrack});
        % keep only those regions longer than a cut off length
        abovelen = cellfun(@(x) size(x,1), abovetracklist);
        keepind = find(abovelen>=minabovelen);
        
        isbelow = ones(size(track,1),1);
        %isbelow(1:indshift) = 0; isbelow(end-indshift+1:end) = 0;
        for tca = 1:length(keepind)
            isbelow(abovetracklist{keepind(tca)}(:,3)) = 0;
        end
    end
    %%    
    belowthresh{tc} = track(find(isbelow),:);
    
    abovethresh = [abovethresh abovetracklist(keepind)];    
%     
%     
%     diffind = diff(aboveind);
%     
%     
%    belowind = find(nvel<=thresh(tc));
%    belowind = [1:indshift,belowind',size(track,1)-indshift+1:size(track,1)];
%    belowthresh{tc} = tracklist{tc}(trackind(belowind),:);
%     
%     abovethresh{tc} = tracklist{tc}(trackind(aboveind),:);    
%     belowthresh{tc} = tracklist{tc}(trackind(belowind),:);
end

%% break up tracks and skipped frames
abovethresh = splitTrackSkip(abovethresh);
belowthresh = splitTrackSkip(belowthresh);

end