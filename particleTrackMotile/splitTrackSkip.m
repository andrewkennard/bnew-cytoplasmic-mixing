function newtracklist = splitTrackSkip(tracklist)
% split up tracks whenever there is a skipped frame
newtracklist = {};
newcnt=0;
for tc = 1:length(tracklist)
    track = tracklist{tc};
    if (isempty(track)); continue; end
    df = diff(track(:,6));
    
    skipind = find(df>1);
    skipind = [skipind; size(track,1)];
    
    newcnt = newcnt +1;
    newtracklist{newcnt} = track(1:skipind(1),:);
    for ic = 2:length(skipind)
        newcnt = newcnt + 1;
        newtracklist{newcnt} = track(skipind(ic-1)+1:skipind(ic),:);
    end
    
end

end