imFile = fullfile('example', ...
       '20141223-Zyla_FN_LysoTracker_100nM_Diff_Z_Planes', ...
       'M1_100x_1.0x_FN_LysoTrackerDeepRed_Ctrl_Cell15_50ms_1 (Z+0) - Viscous_ Good', ...
       'M1_100x_1.0x_FN_LysoTrackerDeepRed_Ctrl_Cell15_50ms_1_MMStack.ome.tif');
imFile2 = fullfile('F:\DATA 01\AndrewK\exampleMovies\mitored_example\SD-mitored-3min.ome.tif');
imFile3 = fullfile('C:\Users\pseudopod\Downloads\SD-mitored-3min.nd2');
DT = getDeltaT(imFile3);

%%
nPlanes = omeMeta.getPlaneCount(0)
timeElapsed = zeros(nPlanes,1);
dt = zeros(nPlanes,1);
for k = 1:nPlanes
    iC = omeMeta.getPlaneTheC(0,k-1).getValue();
    iT = omeMeta.getPlaneTheT(0,k-1).getValue();
    iZ = omeMeta.getPlaneTheZ(0,k-1).getValue();
    if (iC==0) && (iZ == 0)
        timeElapsed(k) = omeMeta.getPlaneDeltaT(0,k-1).value(ome.units.UNITS.SECOND);
    end
end
dt(2:end) = timeElapsed(2:end) - timeElapsed(1:end-1);