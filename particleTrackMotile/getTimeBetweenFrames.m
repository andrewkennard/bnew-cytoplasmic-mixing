function [time, eltime] = getTimeBetweenFrames(imfile)
%%% Get the time between frames in a movie (in seconds)
%%% Input: 
%%%   - imfile: path to an image file that can be read by Bioformats
%%% Output:
%%%   - time: time elapsed since the beginning of the experiment (in
%%%   seconds) with time(1) == 0
%%%   - eltime: time elapsed (in ms) since the beginning of the experiment total
%%%   (including any delay at the beginning of the experiment)

reader = bfGetReader(imfile);
omeMeta = reader.getMetadataStore();

nPlanes = omeMeta.getPlaneCount(0);
sizeT = omeMeta.getPixelsSizeT(0).getValue();
time = zeros(sizeT, 1);
eltime = zeros(sizeT, 1);

for k = 0:nPlanes-1
    iC = omeMeta.getPlaneTheC(0, k).getValue();
    iZ = omeMeta.getPlaneTheZ(0, k).getValue();
    iT = omeMeta.getPlaneTheT(0, k).getValue();
    if (iC == 0) && (iZ == 0)
        disp(k);
        time(iT+1) = omeMeta.getPlaneDeltaT(0, k).value(ome.units.UNITS.SECOND);
        %set first timepoint to time 0
        time(iT+1) = time(iT + 1) - time(1);
        eltime(iT+1) = omeMeta.getPlaneDeltaT(0,k).value(ome.units.UNITS.MILLISECOND);
    end
end
eltime(2:end) = time(2:end) - time(1:end-1);

reader.close()


