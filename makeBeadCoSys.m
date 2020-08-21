
function T_BtG = makeBeadCoSys(pt1,pt2,pt3)

xt = pt2-pt1;
yt = pt3-pt1;
zt = cross(xt,yt);
yt2 = cross(zt,xt);


xu = unit(xt);
yu = unit(yt2);
zu = unit(zt);


% orig = (pt2+pt1)/2 - L_12/2 * xu; % this combines the tracking errors associated with beads 1 and 2 and doesn't bias it toward one or the other
% transform from bead to global
T_BtG = eye(4);
T_BtG(1:3,:) = [xu yu zu pt1];
% T_GtB = invTranspose(T_BtG); 
end
