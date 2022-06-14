function CL = convert2cellObj(cellmv,name, dirname)
% convert a cellmv structure to a single cell object, no particles
%%

CL = CellObj(name);    

img = imread(cellmv(1).fname,cellmv(1).fpage);
imsize  = size(img);

CL.ImgSize = imsize;
CL.DirName = dirname;

%% data on each frame
CL.NFrame = length(cellmv);
CL.Frames = cellmv;
CL.clearParticles;

for fc = 1:CL.NFrame
    CL.Frames(fc).COM = imsize/2;    
    CL.Frames(fc).segcont = [1 1; 1 imsize(1); imsize(2) imsize(1); imsize(2) 1; 1 1];
    CL.Frames(fc).rect = [1 1 imsize(2) imsize(1)];
end
end