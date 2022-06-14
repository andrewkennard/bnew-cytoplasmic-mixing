classdef FeatureObj
    % feature object, corresponding to a single particle on a single frame
    
    properties 
        % coordinates
        XY = []
        % cell index
        CI = 0; 
        % index within Frames array of a cell
        FI = 0;
        % index within list of points (tracks) of a given frame
        FIC = 0;
        % index within array of particle tracks (Particles)
        PI = 0;
        % index within a particular particle track
        PIC = 0;                 
    end
    
    methods
        function FL = FeatureObj(XY)
            % create a feature object with the given coordinates
            FL.XY = XY;
        end
    end
    
end