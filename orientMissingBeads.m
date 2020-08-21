function [ptsB,refPtB,T_GtB,L_CT] = orientMissingBeads(ptsCT,bead3dG,refPtG,plotFlag)

% Use this function when two beads are visible in both cameras, but the
% third is only visible in one. It will take the beads and re-orient them
% into a bead co-ordinate system in preparation to be input to
% "beadLocationProbability.m". It only does ONE frame at a time.


beadnan = isnan(bead3dG(1,:));

bI = find(beadnan == 0);
bN = find(beadnan == 1);

pt1G = bead3dG(:,bI(1));
pt2G = bead3dG(:,bI(2)); 

pt1CT = ptsCT(:,bI(1));
pt2CT = ptsCT(:,bI(2));
pt3CT = ptsCT(:,bN(1)); % the missing bead

refpt1G = refPtG(:,bI(1)); % the beads in the previous frame
refpt2G = refPtG(:,bI(2));
refpt3G = refPtG(:,bN(1)); 


% ideal lengths of the sides
L_12 = norm(pt1CT - pt2CT);
L_13 = norm(pt1CT - pt3CT);
L_23 = norm(pt3CT - pt2CT);

% create the bead co-ordinate system

xt = pt2G-pt1G;
yt = [0 1 0]';
zt = cross(xt,yt);
yt2 = cross(zt,xt);


xu = unit(xt);
yu = unit(yt2);
zu = unit(zt);


orig = (pt2G+pt1G)/2 - L_12/2 * xu; % this combines the tracking errors associated with beads 1 and 2 and doesn't bias it toward one or the other
% transform from bead to global
T_BtG = eye(4);
T_BtG(1:3,:) = [xu yu zu orig];
T_GtB = invTranspose(T_BtG);        % global to bead

pt1B = transformPoints(T_GtB,pt1G);
pt2B = transformPoints(T_GtB,pt2G);
% 
refpt1B = transformPoints(T_GtB,refpt1G); % the beads in the previous frame
refpt2B = transformPoints(T_GtB,refpt2G);
refpt3B = transformPoints(T_GtB,refpt3G); 

refPtB = [refpt1B,refpt2B,refpt3B];

% we know l12, l23 and l13 from the CT scan. 
% determine the angles at the triangle corners at pt1 and pt2

gamma = acosd((L_12^2+L_13^2-L_23^2)/(2*L_12*L_13));
beta = acosd((L_23^2+L_12^2-L_13^2)/(2*L_23*L_12));

if gamma > 90 % pt3 lies in the -x direction
    pt3B(1,1) = -L_13*cosd(180-gamma);
    pt3B(2,1) = L_23*sind(beta);

else % pt 3 x is positive on the x axis
    pt3B(1,1) = L_13*cosd(gamma);
    pt3B(2,1) = L_13*sind(gamma);
  
end

pt3B(3,1) = 0;

ptsB = [pt1B,pt2B,pt3B];

L_CT = [L_12;L_13;L_23];


if plotFlag == 1
    
    figure;
    hold on;
    plot3quick(pt1B,'g','o')
    plot3quick(pt2B,'b','o')
    plot3quick(pt3B,'k','o')
    
    plot3quick(refpt1B,'g','diamond')
    plot3quick(refpt2B,'b','diamond')
    plot3quick(refpt3B,'k','diamond')
    
    plotvector3(refpt1B,vB(:,1),'g')
    plotvector3(refpt2B,vB(:,1),'b')
    plotvector3(refpt3B,vB(:,1),'k')
    
%     plot3quick(refPtB,'k','diamond')
    axis equal
    title('Bead co-ordinate system')
    xlabel('X');ylabel('Y');zlabel('Z');
    grid on
end
end

