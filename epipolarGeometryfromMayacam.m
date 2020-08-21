function epi_geo = epipolarGeometryfromMayacam(camFolder,plotflag)

% determines the epipolar geometry according to Richard Hartley and Andrew
% Zimmerman - Multiple View Geometry in Computer Vision, 2nd Edition

% L. Welte, Sept/2019

% load the mayacam files
cam.c1 = importMayaCam2( [camFolder, ls( [camFolder '*C1S*' ])] );
cam.c2 = importMayaCam2( [camFolder, ls( [camFolder '*C2S*' ])] );


T_c1_WtC = (cam.c1.T_WtC);
T_c1_CtW = invTranspose(T_c1_WtC);
K1 = cam.c1.M_int;
P1 =  cam.c1.M_int * T_c1_WtC(1:3,:);
P1p = P1'*inv(P1*P1'); % the P+
M1 = P1(1:3,1:3);
C1_cent = [[-inv(M1) * P1(:,4)];1];


T_c2_WtC = (cam.c2.T_WtC);
T_c2_CtW = invTranspose(cam.c2.T_WtC);
K2 = cam.c2.M_int;
P2 =  cam.c2.M_int * T_c2_WtC(1:3,:);
P2p = P2'*inv(P2*P2'); % the P+
M2 = P2(1:3,1:3);
C2_cent =[[-inv(M2) * P2(:,4)];1];


T_2to1 = T_c2_WtC*T_c1_CtW ;
R = T_2to1(1:3,1:3);
t = T_2to1(1:3,4);

F = inv(K2)'*makeSkewMatrix(t)*R/(K1);

e1 =(K1*R'*t); % the image coordinates of the epipoles
e2 = (K2*t);

e1 = e1(1:3)/e1(3); % turn to Euclidean co-ords from homogenous
e2 = e2(1:3)/e2(3);

E1n = P1p * e1; % transform to 3D using camera matrix
E2n = P2p * e2;

E1n = E1n(1:3)/E1n(4);
E2n = E2n(1:3)/E2n(4);

prin_vec1 = M1(3,:)'; % principal camera axes
prin_vec2 = M2(3,:)';

PVn1 = prin_vec1;%/(prin_vec1(4));
PVn2 = prin_vec2;%/(prin_vec2(4));

prin_pt1 =  M1 * M1(3,:)';
prin_pt2 = M2 * M2(3,:)';

prin_pt1n = prin_pt1/prin_pt1(3);
prin_pt2n = prin_pt2/prin_pt2(3);
prin_pt1n = P1p * prin_pt1n; % convert to 3D
prin_pt2n = P2p * prin_pt2n ;
prin_pt1n = prin_pt1n/prin_pt1n(4);
prin_pt2n = prin_pt2n/prin_pt2n(4);


epi_geo.e1I = e1; % image co-ordinates of epipoles
epi_geo.e2I = e2;
epi_geo.E1G = E1n; % world co-ordinates of epipoles
epi_geo.E2G = E2n;
epi_geo.P1 = P1; % 3D to 2D image co-ordinate transform
epi_geo.P2 = P2;
epi_geo.P1p = P1p; % pseudo inverse of P1 - 2Dimage  to 3D coords transform
epi_geo.P2p = P2p;
epi_geo.K1 = K1; % intrinsic camera matrix
epi_geo.K2 = K2;
epi_geo.C1G = C1_cent; % camera centre in world coords
epi_geo.C2G = C2_cent;
epi_geo.M1 = M1; % 3x3 matrix within P
epi_geo.M2 = M2;
epi_geo.pvec1G = prin_vec1; % principal camera vector in 3D world coords
epi_geo.pvec2G = prin_vec2;
epi_geo.pPt1G = prin_pt1n; % principal camera point in 3D world coords
epi_geo.pPt2G = prin_pt2n;
epi_geo.F = F; % fundamental matrix

if plotflag == 1
    figure; hold on;
    plot3quick( E1n,'r','+') % epipoles
    plot3quick( E2n,'r','+')
    text(E1n(1),E1n(2),E1n(3),'e1','color','r') % label epipoles
    text(E2n(1),E2n(2),E2n(3),'e2','color','r')
    plotvector3(C1_cent(1:3),PVn1(1:3)*800,'k') % plot principal vectors
    plotvector3(C2_cent(1:3),PVn2(1:3)*800,'k')
    plot3quick(prin_pt1n,'g','x') % plot the principal points
    plot3quick(prin_pt2n,'b','x')
    plot3(C1_cent(1,1),C1_cent(2,1),C1_cent(3,1),'go') % plot camera centres
    plot3(C2_cent(1,1),C2_cent(2,1),C2_cent(3,1),'bo')
    text(C1_cent(1,1),C1_cent(2,1),C1_cent(3,1)-100,'C1','color','g') % add text to camera centres
    text(C2_cent(1,1),C2_cent(2,1),C2_cent(3,1)-100,'C2','color','b')
    plot3([C1_cent(1,1);C2_cent(1,1)],[C1_cent(2,1);C2_cent(2,1)],[C1_cent(3,1);C2_cent(3,1)],'r') % plot the baseline
    
    % plot the image planes - cam 1
    imPl = @(P,X,Z,pt) (P(1) * (X - pt(1)) - P(2) *  pt(2) + P(3)*(Z - pt(3)) ) / (- P(2));
    
    x = [-1000 1000];
    z = [-1000 1000];
    
    [X,Z] = meshgrid(x,z);
    Y1 = imPl(P1(3,:),X,Z,prin_pt1n) ;
    reOrder = [1 2 4 3];
    patch(X(reOrder),Y1(reOrder),Z(reOrder),'g');
    alpha(0.2);
    
    % plot the image planes - cam 2
    x = [-200 200];
    z = [-1000 500];
    
    [X,Z] = meshgrid(x,z);
    Y2 = imPl(P2(3,:),X,Z,prin_pt2n);
    patch(X(reOrder),Y2(reOrder),Z(reOrder),'b');
    alpha(0.2);
    
    axis equal
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
end


