# bnew-cytoplasmic-mixing
Update of code for BNEW analysis of cytoplasmic mixing (edits to facilitate different input file formats)

See also https://github.com/lenafabr/BNEW

See also Koslover LF _et al._ (2017). Cytoplasmic Flow and Mixing Due to Deformation of Motile Cells. _Biophys. J._ 113(9): 2077-2087. https://www.sciencedirect.com/science/article/pii/S0006349517310202 


# To Run
To get you started, got into this script and go through it section-by-section:
particleTrackMotile/example/process_example.m
This will load an example movie, let you manually adjust some segmentation parameters, segment the cell contour, extract particle trajectories, save everything in a separate .mat file for each individual cell. A couple of example movies are included to practice on.

Separately, if you want to start with trajectories and plot them or play around with them, go to
particleTrackMotile/viewLysotrackerTrajectories.m

At the top you can choose to look at the trajectories from the example movies you processed, or at the full set of trajectories from all the lysotracker movie data I have. Those are stored in individual .mat files (one per cell) in the directory Cells_HL60_intpix
