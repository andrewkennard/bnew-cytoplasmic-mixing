function newtracklist = breakupTracks(tracklist,badind,varargin)
% break up tracks, excluding designated set of frames
% ie: all tracks are broken at those frames and data for those frames is
% not included

% This is currently written to run quite slow (does not take consecutive
% stretches of frames into account)

%badind = [1:94, 312:717];
%badind = [4:5,10:15,41:57];
%badind = 100:110
%%
% put the numbering label for the continuous good region in column 9
labelregion = 1;
if (nargin>2)
    for vc = 1:2:length(varargin)
        switch varargin{vc}
            case('labelregion')
                labelregion = varargin{vc+1};
        end
    end
end

maxfr = max(cellfun(@(x) x(end,6),tracklist));

% number continuous good stretches
if (isempty(badind))
    goodstretchnum = ones(1,maxfr);
else
    diffs = diff(badind);
    ct = 0;
    goodstretchnum = zeros(1,maxfr);
    if (badind(1)>1)
        ct = ct+1;
        goodstretchnum(1:badind(1)-1) = ct;
    end
    for c = 1:length(diffs)
        if (diffs(c)>1) % continuous good stretch
            ct = ct+1;
            inds = (badind(c)+1):(badind(c+1)-1);
            goodstretchnum(inds) = ct;
        end
    end
    if (badind(end)<maxfr)
        ct = ct+1;
        goodstretchnum(badind(end)+1:maxfr) = ct;
    end
end

newtracklist = {};
newct = 0;

for tc = 1:length(tracklist)
    %%
    track=tracklist{tc};
    goodind = find(~ismember(track(:,6),badind));
    % track does not overlap with good region
    if (isempty(goodind)); continue; end 
    
    % break up discontinuous stretches
    dind = diff(goodind);
    
    brind = find(dind>1);
    
    if (isempty(brind))
        % single chunk of track falls in good region
        newtracklist{newct+1} = track(goodind,:);
    else 
        newtracklist{newct+1} = track(goodind(1:brind(1)),:);
        newct = newct+1;
        for bc = 2:length(brind)        
            newtracklist{newct+1} = track(goodind(brind(bc-1)+1:brind(bc)),:);
            newct = newct+1;
        end
        newtracklist{newct+1} = track(goodind(brind(end)+1:end),:);
    end            
    newct = newct+1;  
end

% stick the numbered good stretch into column 9
if (labelregion)
    for tc = 1:newct
        newtracklist{tc}(:,9) = goodstretchnum(newtracklist{tc}(:,6));
    end
end
