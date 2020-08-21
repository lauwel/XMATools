function [beadB,rotVal] = beadLocationProbability(limits , N,  epi_geoB, ptsB,refPtB, sigma,epi_error,analysisType,rotValInit,plotFlag)
%       beadB = beadLocationProbability(limits , N,  epi_geo, posB,refB, sigma)
%
% This function determines the probable location of a tantalum bead based
% on the known co-ordinates in CT space and knowledge of the bead in one
% camera.
% The beads must be put into the bead co-ordinate system which is defined
% by the two known 3D positions of the other beads. These beads form the
% x-axis in bead space, and the third, unknown, bead makes a triangle off
% that axis. Use "orientMissingBeads.m" to do this.
%
% --------------------------INPUTS-----------------------------------------
%       limits      =  [xlims;ylims;zlims] 3x2 array containing the
%                       minimum and maximum sampling limits in each
%                       dimension
%       N           =   [nx; ny; nz] 3x1 array containing the number of
%                       sample points ineach dimension
%       epi_geoB     =  structure containing the epipolar geometry in bead space.
%                       epi_geoB.epipole is a 3x1 point of the epipole in
%                       bead space
%                       epi_geoB.plane_norm is a 3x1 vector of the epipolar
%                       plane in bead space
%       ptsB        =   3x1 array with the x and y co-ordinate of the bead
%                       of interest relative to the origin in the bead
%                       co-ordinate system. (z coord is defined as 0)
%       refPtB      =   3x3 array with the x, y and z co-ordinates in
%                       rows and the beads in the previous frame in columns 
%                       1 through 3 (i.e. dim x bead).
%       sigma       =   3 x 1 array containing the standard deviations
%                       for each gaussian function - [sx;sr;sepi]. sx
%                       is along the cylinder, sr is the radial direction
%                       of the cylinder, sepi is standard deviation of the
%                       distance from the epipolar plane.
%     epi_error     =   error off of the epipolar plane, due to errors in
%                       triangulation. Give in mm. 
%   analysisType    =   either 'general' or 'refine'. If 'general' is
%                       selected, then ALL solutions will be found and the
%                       optimal solution selected using the position of the
%                       beads in the previous frame. If 'refine' is
%                       selected, it's assumed that only one solution is
%                       contained within the solution space.
%     plotFlag  =       0 or 1; plot the solution after computation
% --------------------------OUTPUTS----------------------------------------
%       beadB       =   The probable bead location in bead space.
%
% ------------------------------------------------------------------------

% if epi_error < 0.25
%     epi_error = 0;
% end
if isempty(rotValInit)
    rotValInit = 0;
end
[nx,ny,nz] = deal(N(1),N(2),N(3));
[sx,sr,sepi] = deal(sigma(1),sigma(2),sigma(3));
[x_lims,y_lims,z_lims] = deal(limits(1,:),limits(2,:),limits(3,:));

gauss = @(x,mx,sigma) 1/sqrt(2*pi*sigma^2) * exp(-(x-mx)^2/(2*sigma^2)); % gaussian function : x is value, mx is mean, sigma is standard deviation

xs = zeros(nx,ny,nz); ys = zeros(nx,ny,nz); zs = zeros(nx,ny,nz);
ge = zeros(nx,ny,nz); gr_save = zeros(nx,ny,nz); gx = zeros(nx,ny,nz); 

pt3B = ptsB(:,3);

T_BtG_ref = makeBeadCoSys(refPtB(:,1),refPtB(:,2),refPtB(:,3));

ly = 0; lz = 0;                                     % centre of the donut

indy = 0;
for y = linspace(y_lims(1),y_lims(2),ny)
    indy = indy+1;
    indz = 0;
    
    for z = linspace(z_lims(1),z_lims(2),nz)
        
        indz = indz+1;
        r = norm([y,z]-[ly,lz]);
        
        
        gr = gauss(r,pt3B(2),sr);                % radial gaussian with mean at the y co-ordinate
        
        indx = 0;
        for x = linspace(x_lims(1),x_lims(2),nx)
            indx = indx+1;
            
            gr_save(indx,indy,indz) = gr;
            gx(indx,indy,indz) = gauss(x,pt3B(1),sx);            % x dimension gaussian with mean at the x coordinate of the missing bead
            
            xs(indx,indy,indz) = x;
            ys(indx,indy,indz) = y;
            zs(indx,indy,indz) = z;
            
            % look at the distance between the epipolar plane and the point in question
            pt_epi = closestPointonPlanealongVector([x,y,z]',epi_geoB.plane_norm,epi_geoB.epipole,epi_geoB.plane_norm);
           
   
            dist_epi = norm(pt_epi-[x,y,z]');        % distance from epipolar plane and point in question
            ge(indx,indy,indz) = gauss(dist_epi,epi_error,sepi);           % determine the value of the gaussian function at this point
            
            
            
        end
        
    end
end



gt = gr_save .* gx .* ge;
% epi_error

% if strcmp(analysisType,'refine') % we are refining a SINGLE solution
%     [~,I] = max(gt,[],'all','linear');
%     [ix_max,iy_max,iz_max] = ind2sub(N',I);
%     % % get the most probable bead location
%     beadB = [xs(ix_max,iy_max,iz_max);...
%         ys(ix_max,iy_max,iz_max);...
%         zs(ix_max,iy_max,iz_max)];
%     
%     T_fr = makeBeadCoSys( ptsB(:,1), ptsB(:,2),beadB);
%     T_test = invTranspose(T_BtG_ref) * T_fr;
%     rotVals = convertRotation(T_test,'4x4xn','helical');
%     rotVals(1)
%     plot3quick_scatter(beadB);
% elseif strcmp(analysisType,'general') % we need to find all  solutions
    
%     if epi_error < 0.1 % we need to find 2 solutions
%         numGroups = 2;
%     else
%         numGroups = 4; % we need to find 4 solutions (i.e. the bead is not exactly on the epipolar plane)
%     end
if strcmp(analysisType,'refine')
    perc = 0.99;
    if epi_error < 0.05
        numGroups = 1;
    else
        numGroups = 2;
    end
elseif strcmp(analysisType,'general')
    perc = 0.99;
    numGroups = 2;
end

X = [];
    while (length(X) < 120)
        
           
        Id = gt > max(gt,[],'all')*perc;
    
        X = xs(Id);
        
        if length(X) < 120
        perc = perc-0.05; 
        end
    end
    
    Y = ys(Id);
    Z = zs(Id);
    
    gtI = gt(Id);
    [Idx] = kmeans([X,Y,Z],numGroups,'Distance','sqeuclidean',...
    'Replicates',5);
% 
%     figure;  hold on;
    rotArray = [];
    
    for i = 1:numGroups
        
        xA{i} = X(Idx==i);
        yA{i} = Y(Idx==i);
        zA{i} = Z(Idx==i);
        [~,ImaxSol{i}] = max(gtI(Idx==i));
        
        sol{i} = [xA{i}(ImaxSol{i});yA{i}(ImaxSol{i});zA{i}(ImaxSol{i})];
        
        T_fr{i} = makeBeadCoSys( ptsB(:,1), ptsB(:,2),sol{i});
        T_test = invTranspose(T_BtG_ref) * T_fr{i};
        rotVals{i} =      convertRotation(T_test,'4x4xn','helical');
        rotArray = [rotArray, rotVals{i}(1)];
        
%         
%     plot3(X(Idx==i),Y(Idx==i),Z(Idx==i),'.')
%     plot3quick_scatter(sol{i});
    end

    rotArray
%     view([-39 59])
%     drawnow
    
    [rotVal,solInd] = min(abs(rotArray-rotValInit));

    beadB = sol{solInd};

     


if plotFlag == 1
    
    figure
    %     indg = find(gt_col >60);
    % patch(xs(indg),ys(indg),zs(indg),gt_col(indg))
    hold on
    cmap = colormap('jet');
    
    %     ivstring = createInventorHeader();
    gt_col = round(gt/max(max(max(gt))) * 62) + 1;
    for i = 1:1:length(gt(:))
        if gt_col(i)>50
            h = plot3(xs(i),ys(i),zs(i),'.');
            set(h,'color',cmap(gt_col(i),:))
            
            %                     ivstring = createInventorSphere([xs(ix,iy,iz),ys(ix,iy,iz),zs(ix,iy,iz)],0.1,cmap(gt_col(ix,iy,iz),:),0.5);
        end
    end
%     plot3(xs(I),ys(I),zs(I),'ko')
    
    plot3(x1(ImaxSol1),y1(ImaxSol1),z1(ImaxSol1),'ko')
    plot3(x2(ImaxSol2),y2(ImaxSol2),z2(ImaxSol2),'ko')
    %
    %     for ix = 1:4:indx
    %         for iy = 1:4:indy
    %             for iz = 1:4:indz
    %                 if gt_col(ix,iy,iz)>50
    %                     h = plot3(xs(I),ys(I),zs(I),'.');
    %                     set(h,'color',cmap(gt_col(I),:))
    %
    % %                     ivstring = createInventorSphere([xs(ix,iy,iz),ys(ix,iy,iz),zs(ix,iy,iz)],0.1,cmap(gt_col(ix,iy,iz),:),0.5);
    %                 end
    %             end
    %         end
    %     end
    plot3quick(beadB,'k','o')
    %     plotvector3(epi_geoB.epipole,epi_geoB.plane_norm*50,'m')
    plot3quick(ptsB(:,1),'k','o')
    plot3quick(ptsB(:,2),'k','o')
    plot3quick(ptsB(:,3),'r','o')
    
    plot3quick(refPtB,'k','diamond')
    imPl = @(P,X,Z,pt) (P(1) * (X - pt(1)) - P(2) *  pt(2) + P(3)*(Z - pt(3)) ) / (- P(2));
    
    x = [20 40];
    z = [0 10];
    
    [X,Z] = meshgrid(x,z);
    Y1 = imPl(epi_geoB.plane_norm,X,Z,epi_geoB.epipole) ;
    reOrder = [1 2 4 3];
    patch(X(reOrder),Y1(reOrder),Z(reOrder),'m');
    alpha(0.2);
    axis equal
    
    
    
end
