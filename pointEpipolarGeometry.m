function [epi_planeG,epi_lineG] = pointEpipolarGeometry(epi_geo,pt2DI,cam)

% This function will take a two dimensional camera point in homogenous
% image co-ordinates i.e. (x,y,1)', and then compute the epipolar line,
% and plane normal in world co-ordinates
% ---------------------------INPUTS---------------------------------------
% epi_geo   =  computed from the mayacams and
%               "epipolarGeometryfromMayaCam.m"
% pt2DI     =   3x1 homogenous co-ordinate of the point in the image
%               co-ordinates.
% cam       =   either scalar 1 or 2. Gives the camera number of the image 
%               point, to calculate the geometry for the other camera view 
%               (i.e. put 1 if the image point is in 1, to calculate geometry in 2).
% --------------------------OUTPUTS---------------------------------------
% epi_planeG =  normal of the epipolar plane containing the 2D image point
% epi_lineG  =  the epipolar line in world co-ordinates of the 2D point.
%               Pre-multiply by the camera projection matrix (P) to convert
%               to the relevant 2D image.
% 
% ------------------------------------------------------------------------
% All computations are derived from Multiple View Geometry, 2nd Edition.
% Hartley and Zisserman.
% L. Welte Sept/2019


if cam == 1
    epi_lineI = epi_geo.F * pt2DI; % the epipolar line in image 2
    epi_planeG = epi_geo.P2'*epi_lineI;
%     epi_planeG = epi_planeG/epi_planeG(4);
    epi_planeG = unit(epi_planeG(1:3,1));
elseif cam == 2
    
    epi_lineI = epi_geo.F' * pt2DI; % the epipolar line in image 1
    epi_planeG = epi_geo.P1' * epi_lineI; % the plane of line solutions
%     epi_planeG = epi_planeG/epi_planeG(4);
    epi_planeG = unit(epi_planeG(1:3,1));
else
    
    error('Incorrect value for cam variable in pointEpipolarGeometry. Select either 1 or 2.')
end

EPP = cross(epi_geo.pvec1G,epi_planeG(1:3)); % cross of principal vector and epipolar plane normal 
epi_lineG  = EPP/norm(EPP); % the epipolar (line) vector - it originates at E1n and is in the image plane
