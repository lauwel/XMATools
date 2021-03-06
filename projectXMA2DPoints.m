function [pos3D] = projectXMA2DPoints(inFile,camFolder,subj_dir,trial_name,filter,options,bones,anim_dir)
% projectXMA2DPoints(inFile,camFolder,subj_dir,trial_name,fc,options,bones,anim_dir)
% 
% This function will project undistorted 2D points from an XMA lab file and
% will determine their positions in 3D. 
% 
% It requires matchXMAtoCT.m which will  match the 3D bead positions to the 
% beads in CT space, and will compute the transforms for each frame from CT 
% to x-ray space. The transforms can be output as a variable 'T', but are 
% also saved in a .mat file. The autoscoper .tra files for each bone, as 
% well as a combined "allBone" trial is saved. The bones are also animated 
% in Wrist Visualizer.
% 
% 
% if you don't want an argument, just put a [].
% -------------  INPUT ARGUMENTS ----------------------------------
% inFile    
% 
% options       =   either 'none' or 'optimize'. Optimize option will use
%                   the bead optimizer for missing frames
% ^  change to an optimize structure to have a Bones, 'On' etc structure
% filter        =  either [] OR filter.Type = 'none' if no filtering 
%                   is desired, otherwise, should have fields of 
%                   - Type (either 'none','2D','3D','2D-3D') which filters
%                   the 2D beads, or the 3D bead positions, or both
%                   - Fs  sample frequency in Hz (numeric)
%                   - Fc  cut-off frequency (two values for the adaptive
%                   butterworth filter)
% 
% 
addpath(genpath('P:/Code/MATLAB_SOL/'))



if ~exist('options','var') 
    options = 'none';
end




% ivDir = 'E:\SOL001_VISIT2\Models\IV\Beads\'; % where the bead IV
datBeadLocs = csvread('E:\SOL001_VISIT2\Models\bead_positions.txt',1,1);
dat = readtable('E:\SOL001_VISIT2\Models\bead_positions.txt');


[pos2D,nBones,nBeads] = loadXMA2dPoints(inFile);
nfr = size(pos2D.tib(1).cam1,1);

%% Load all the bone positions


beadStructTemp = table2struct(dat);
nbeads = length(beadStructTemp);
for b = 1:nbeads
    
    fieldname = beadStructTemp(b).bead;
    beadnum = str2double(fieldname(4));
    beadCT.(fieldname(1:3))(:,beadnum) = datBeadLocs(b,:);
    
end
bonesCell = fields(beadCT);
bone_list = char(bonesCell);


if ~exist('bones','var')
    
    bone_ind = 1:nBones;
elseif ischar(bones)
    bone_ind = find(strcmp(bonesCell,bones));
elseif isnumeric(bones) && ~isempty(bones)
    bone_ind = bones;
elseif isempty(bones)
    
    bone_ind = 1:nBones;
else
    error('The variable ''bones'' is not of supported type.');
end
%%  Determine the 3D positions of the beads. If one is missing from one
% image, then calculate the epipolar geometry and relevant plane

epi_geo = epipolarGeometryfromMayacam(camFolder,0);
res = (2048/mean([epi_geo.K1([1,5]) epi_geo.K2([1,5])])); % camera resolution

ind = 1;

for bn = 1:nBones
%     rvI(1:nfr) = 0.25; % initial rotation value;
    ind = 1;
    optimizeFrames = [];
    for fr = 1:nfr
        
        
        ptsCT = beadCT.(bone_list(bn,:));
        
        bead2d = [];
        
        nBeadsbone = size(pos2D.(bone_list(bn,:)),2);
        
        for bd = 1:nBeadsbone
            if fr == 1 && ~isempty(fc)
                pos2D.(bone_list(bn,:))(bd).cam1 = adaptiveLowPassButterworth(pos2D.(bone_list(bn,:))(bd).cam1',[20 50],250,0)';
                pos2D.(bone_list(bn,:))(bd).cam2 = adaptiveLowPassButterworth(pos2D.(bone_list(bn,:))(bd).cam2',[20 50],250,0)';
            end
            bead2d(bd,:,1) = [pos2D.(bone_list(bn,:))(bd).cam1(fr,:) 1];
            bead2d(bd,:,2) = [pos2D.(bone_list(bn,:))(bd).cam2(fr,:) 1];
        end
        
        % loop through and either triangulate if both are tracked, or
        % assign NaN values to the 3D point
        bead3dG = nan(3,nBeadsbone);
        
        ind_tracked = ~isnan(bead2d(:,1,1)) & ~isnan(bead2d(:,1,2)); % beads tracked in both camera views
        
        for bd =  find(ind_tracked)'             % loop through the tracked beads
            [bead3dG(1:3,bd),err(bd)] = triangulate(bead2d(bd,1:2,1),bead2d(bd,1:2,2),epi_geo.P1',epi_geo.P2');
            
        end
        
        if (sum(ind_tracked) >= 2 ) && ~isempty(find(bn==bone_ind)) % if this frame can be optimized
            for bd = find(~ind_tracked)' % find the missing bead(s)
                
                camflag = find(isnan(bead2d(bd,1,1:2))); % find the camera that's missing a view
                if length(camflag) == 1
                    trackedflag = find(~isnan(bead2d(bd,1,1:2))); %
                    epi_error = res *mean(err(ind_tracked));
                    optimizeFrames(ind,:) = [fr bd camflag epi_error bead2d(bd,:,trackedflag)]; %[ frame, missing bead, missing camera, error in epipolar plane, 2d coords]
                    ind = ind+1;
                end
            end
        end
        
        pos3D.(bone_list(bn,:))(:,:,fr)= bead3dG;
        
        
    end
    
    if (strcmp(options,'optimize'))
        
        for i = 1:size(optimizeFrames,1)
            
            fr = optimizeFrames(i,1); % for all the frames that can be optimized
            bd = optimizeFrames(i,2);
            camflag = optimizeFrames(i,3);
            epi_error = optimizeFrames(i,4);
            bead2d = optimizeFrames(i,5:7);
            bead3dG = pos3D.(bone_list(bn,:))(:,:,fr);
            %         if (sum(ind_tracked) >= 2 )  && (strcmp(options,'optimize')) && ~isempty(find(bn==bone_ind)) % there are AT LEAST two fully tracked beads, we can optimize
            
            
            
            %             for bd = find(~ind_tracked)'
            
            %                 camflag = find(isnan(bead2d(bd,1,1:2)));
            %                 epi_error = res *mean(err(ind_tracked));
            
            
            %                 if length(camflag) == 1         % only  one view is missing the bead
            
            %                     beadflag = bd;
            
            if camflag == 1             % MISSING bead in CAM1
                
                cam = 2;
                [epi_planeG,epi_lineG] = pointEpipolarGeometry(epi_geo,bead2d',cam);
                epi_poleG = epi_geo.E1G;
                
                
            elseif camflag == 2         % MISSING bead in CAM2
                cam = 1;
                [epi_planeG,epi_lineG] = pointEpipolarGeometry(epi_geo,bead2d',cam);
                epi_poleG = epi_geo.E2G;
                
            end
            %                 else                            % go to the next bead if both camera views are missing
            %                     continue
            %                 end
            frsTracked = find(~isnan(squeeze(pos3D.(bone_list(bn,:))(1,bd,:))));
            frsDiff = abs(frsTracked-fr);
            [M,I] = min(frsDiff);
            frRef = frsTracked(I);
            rvI = M * 0.25;
%             if ~isnan( pos3D.(bone_list(bn,:))(1,bd,fr-1))% verify the previous frame
                refPtG =  pos3D.(bone_list(bn,:))(:,:,frRef);
%             elseif ~isnan( pos3D.(bone_list(bn,:))(1,bd,fr+1))% verify the previous frame
%                 refPtG =  pos3D.(bone_list(bn,:))(:,:,fr+1);
%             end
            
            [ptsB,refPtB,T_GtB,L_CT] = orientMissingBeads(ptsCT,bead3dG,refPtG,0); % convert the points to bead co-ordinate system
            
            epi_geoB.epipole = transformPoints(T_GtB,epi_poleG);          % convert the epipole to bead coords
            epi_geoB.plane_norm = transformVectors(T_GtB,epi_planeG);        % convert the epipolar plane normal to bead coordinates
            
            x_lims = ptsB(1,3) * [0.90 1.10];                               % where to extend the probability cylinder
            y_lims = ptsB(2,3) * [-1.05 1.05];                              % needs to extend through the entire negative and positive space
            z_lims = y_lims;
            
            limits = [x_lims;y_lims;z_lims];
            sigma = [2,2,0.2];
            N = [80,80,80]';
            [beadB,~] = beadLocationProbability(limits , N,  epi_geoB, ptsB,refPtB, sigma,epi_error/2,'general',rvI,0);
            
            % run and then refine around the predicted bead space
            %                     rngs = [0.95,1.05];
            N =  [100;100;100] ;
            sigma = [1,1,0.2];
            limits = [beadB-1.5,beadB+1.5] ;
            [beadB,~] = beadLocationProbability(limits ,N, epi_geoB, ptsB,refPtB, sigma,epi_error/2,'refine',rvI,0);
            bead3dG(:,bd) = transformPoints(T_GtB,beadB',-1);
            
            
            % evaluate the found point
            pt = closestPointonPlanealongVector(beadB, epi_geoB.plane_norm,epi_geoB.epipole, epi_geoB.plane_norm);
            epiplaneDist = norm(pt-beadB);
            fprintf('Frame %i\nBone: %s.\n',fr,bone_list(bn,:))
            fprintf('Distance from epipolar plane is %0.3f mm \n',epiplaneDist)
            fprintf('Distance difference from CT beads is L13 %0.3f mm,L23 %0.3f mm . \n\n',L_CT(2)- norm(beadB-ptsB(:,1)),L_CT(3) - norm(beadB-ptsB(:,2)))
            %                     savedat.(bone_list(bn,:))(bd).results(ind,:) = [fr , epiplaneDist,abs(L_CT(2)- norm(beadB-ptsB(:,1))),abs(L_CT(3)- norm(beadB-ptsB(:,2))) ];
            %                     ind = ind+1;
            pos3D.(bone_list(bn,:))(:,bd,fr)= bead3dG(:,bd);
        end
    end
    
    
% pos3D.(bone_list(bn,:))(:,:,fr)= bead3dG;
end



% end % looping frames
% 
% 
% 
% end % looping bones

matchXMAtoCT(pos3D,beadCT,subj_dir,trial_name)
if ~isempty(fc) && isempty(anim_dir)
    
    if isnumeric(fc) && length(fc) == 2
        matchXMAtoCT(pos3D,beadCT,subj_dir,trial_name,fc)
    end
elseif ~isempty(fc) && ~isempty(anim_dir)
     if isnumeric(fc) && length(fc) == 2
        matchXMAtoCT(pos3D,beadCT,subj_dir,trial_name,fc,anim_dir)
    end
end

% write a new file so that it doesn't have to be recomputed every time 
outFile = [inFile(1:end-4) '_optimized.csv'];
writeXMA2Dfrom3DPoints(inFile,outFile,epi_geo,pos3D)

