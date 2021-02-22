function [pos3D] = processBeadData(in2DFile,beadFile,camDir,subjDir,trialName,filterOpts,optimizeOpts,saveOpts)
% This function will process the beaded data. It will take the 2D file from
% XMA lab, convert it into 3D positions, filter the data where you want it
% to, and write transforms/POS files/autoscoper files.
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
% in2DFile      = the undistorted 2D points file you want to process (string)
%                   if [], will prompt the user to select the file
% beadFile      = the CT positions of the beads in a text file (string)
% camDir        = the path to the relevant Mayacam folder (string
% subjDir       = the subject's base path
% filterOpts    = Options: either [] OR filter.Type = 'none' if no filtering 
%                   is desired, otherwise, is a structure with fields of 
%                   - Type (either 'none','2D','3D','2D3D') which filters
%                   the 2D beads, or the 3D bead positions, or both
%                   - Fs  sample frequency in Hz (numeric)
%                   - Fc  cut-off frequency (two numeric values for the adaptive
%                   butterworth filter)
% 
% optimizeOpts  = Options: either [] or optimizeOpts.Type = 'none' if none of
%                   the beads need to be optimized for missing views
%                   Otherwise is a structure with fields:
%                   - .Type: string: 'optimize','none'
%                   - .Bones : cell array with bones to optimize, with three
%                   letter code; e.g.  {'tal','tib'}
%                 
% saveOpts      = Options for saving to a specific location (such as the
%                   server). Either needs to be [] or be a structure that
%                   contains one or many of the following fields:
%                   - animDir = file path to location to save POS file
%                       -> must be specified with ivDir to link the rigid
%                       iv files
%                   - ivDir = location of the iv files to link to.
%               It's set up this way to accommodate future modifications.


% inFile    
% 
% options       =   either 'none' or 'optimize'. Optimize option will use
%                   the bead optimizer for missing frames
% ^  change to an optimize structure to have a Bones, 'On' etc structure
% filter        = 
% 
addpath(genpath('P:/Code/MATLAB_SOL/'))
% subjDir = 'E:SOL001_VISIT2\';
% beadFile = 'E:\SOL001_VISIT2\Models\bead_positions.txt';
% in2DFile = 'E:\SOL001_VISIT2\T0005_SOL001_nwalk_pref_barefoot\XMA_csv\T0005_SOL001_nwalk_pref_barefootUNDISTORTED.csv';

%% check the input variables

% get the subject directory
if isempty(subjDir)
    uiwait(msgbox('Please select the subject''s data directory.')) 
    [subjDir] = uigetdir('*');
end

if isempty(camDir)
    uiwait(msgbox('Please select the subject''s calibration directory.')) 
    [camDir] = uigetdir([subjDir '*']);
end

% get the file to process
if isempty(in2DFile)
   uiwait(msgbox('Please select the undistorted 2D XMA points file to process.')) 
    [in2DFile,in2DDir] = uigetdir([subjDir '*.csv']);;
end



if isempty(beadFile)   
    uiwait(msgbox('Please select the 3D CT bead position file to process.')) 
    [beadFile,beadDir] = uigetfile([subjDir '*']);
end
% validate the filter options
if isempty(filterOpts)
    filterOpts.Type = 'none';
end
filterTypes = {'none','2D','3D','2D3D'};
if ~ismember(filterOpts.Type,filterTypes)
    error('FilterOpts.Type is not a valid Type. Please choose from ''none'',''2D'',''3D'',''2D3D''.')
end

if ismember(filterOpts.Type,filterTypes(2:4)) % if it selects a filter, check all the inputs
    if isfield(filterOpts,'Fs')
        if ~isscalar(filterOpts.Fs)
            error('Please provide a valid sampling frequency (filterOpts.Fs).');
        end
    else
        error('Please include the sampling frequency in the filter options (filterOpts.Fs)')
    end
    if isfield(filterOpts,'Fc')
        if any(size(filterOpts.Fc) ~= [1 2]) % if it's not the right size
            error('Please provide a valid [w1 w2] cut-off frequency (filterOpts.Fc).');
        end
    else
        error('Please include the sampling frequency in the filter options (filterOpts.Fc)')
    end
end

% validate the optimizer options

if isempty(optimizeOpts)
    optimizeOpts.Type = 'none';
end


optimTypes = {'none','optimize'};
if ~ismember(optimizeOpts.Type,optimTypes)
    error('OptimizeOpts.Type is not a valid Type. Please choose from ''none'' or ''optimize''.')
end

if ismember(optimizeOpts.Type,optimTypes(2)) % if it selects an optimizer, check all the inputs
    if isfield(optimizeOpts,'Bones')
        if ~iscell(optimizeOpts.Bones)
            error('Please provide a cell array of bones in the optimizeOpts variable.');
        end
    else
        error('Please include a cell array of bones in the optimizeOpts variable.')
    end
    if ismember(filterOpts.Type,filterTypes(2:4))
        if isfield(filterOpts,'Fc')
            if any(size(filterOpts.Fc) ~= [1 2]) % if it's not the right size
                error('Please provide a valid [w1 w2] cut-off frequency (filterOpts.Fc).');
            end
        else
            error('Please include the sampling frequency in the filter options (filterOpts.Fc)')
        end
    end
end

% validate the save directories
if isempty(saveOpts)
    saveOpts.animDir = '';
end


%% Load the 2D points and 3D positions of the beads in CT

%load the 2D points and filter them if necessary
[pos2D,nBones,nBeads] = loadXMA2dPoints(in2DFile);
bonesCell = fields(pos2D);
nfr = size(pos2D.(bonesCell{1})(1).cam1,1);
if strcmp(filterOpts.Type(1:2),'2D')  % if the 2D points are being filtered
    for bn = 1:nBones    
        nBeadsbone = size(pos2D.(bonesCell{bn}),2);
        for bd = 1:nBeadsbone
            bead_init = [pos2D.(bonesCell{bn})(bd).cam1';pos2D.(bonesCell{bn})(bd).cam2'];
            
            bd1i = find(~isnan(bead_init(1,:)));
            bd2i = find(~isnan(bead_init(3,:)));
%             fprintf('Data from bone %s and bead %i. ', bonesCell{bn},bd)
            if (length(bd1i) <= 5) || (length(bd2i)<=5)
                warning('Missing data from bone %s and bead %i. ', bonesCell{bn},bd)
                continue
            end
                pos2D.(bonesCell{bn})(bd).cam1 = adaptiveLowPassButterworth(pos2D.(bonesCell{bn})(bd).cam1',filterOpts.Fc,filterOpts.Fs,0)';
                pos2D.(bonesCell{bn})(bd).cam2 = adaptiveLowPassButterworth(pos2D.(bonesCell{bn})(bd).cam2',filterOpts.Fc,filterOpts.Fs,0)';
        end
    end
end

% load the static bead positions and put in a beadCT structure
datBeadTable = readtable(beadFile);
beadPos = csvread(beadFile,1,1);
beadStructTemp = table2struct(datBeadTable);
for b = 1:nBeads
    fieldname = beadStructTemp(b).bead;
    beadnum = str2double(fieldname(4));
    beadCT.(fieldname(1:3))(:,beadnum) = beadPos(b,:);
end
bonesCTCell = fields(beadCT);
nBonesCT = length(bonesCTCell);


% make sure that there are the same beads in both the 2D points file as the
% CT bead positions

bonesRef = intersect(bonesCell,bonesCTCell); % list of bones that are found in both files

if isempty(bonesRef)
    error('No bead names match the beads in the CT position file.');
else % remove all bones that don't have correspondances  
    for bn = 1:nBones
        if ~ismember(bonesCell{bn},bonesRef)
            pos2D = rmfield(pos2D,bonesCell{bn});
            warning('The %s bone does not appear in the 2D bead position file. Skipping...',bonesCell{bn})
        end
    end
    
    for bn = 1:nBonesCT
         if ~ismember(bonesCTCell{bn},bonesRef)
            beadCT = rmfield(beadCT,bonesCTCell{bn});
            warning('The %s bone does not appear in the CT bead position file. Skipping...',bonesCTCell{bn})
         end
    end
         
end % remove beads that are not found in both files

if strcmp(optimizeOpts.Type,'optimize')
    bonesOptimize = intersect(optimizeOpts.Bones,bonesRef);
else
    bonesOptimize = {};
end

nBones = length(bonesRef);


%%  Determine the 3D positions of the beads. If one is missing from one
% image, then calculate the epipolar geometry and relevant plane

epi_geo = epipolarGeometryfromMayacam(camDir,0);
res = (2048/mean([epi_geo.K1([1,5]) epi_geo.K2([1,5])])); % camera resolution


for bn = 1:nBones
    ind = 1; 
    optimizeFrames = [];
    for fr = 1:nfr
        
        bead2d = [];
        ptsCT = beadCT.(bonesRef{bn});
        nBeadsbone = size(pos2D.(bonesRef{bn}),2); % number of beads in the bone
        
        for bd = 1:nBeadsbone
            bead2d(bd,:,1) = [pos2D.(bonesRef{bn})(bd).cam1(fr,:) 1];
            bead2d(bd,:,2) = [pos2D.(bonesRef{bn})(bd).cam2(fr,:) 1];
        end
        
        % loop through and either triangulate if both are tracked, or
        % assign NaN values to the 3D point
        bead3dG = nan(3,nBeadsbone);
        
        ind_tracked = ~isnan(bead2d(:,1,1)) & ~isnan(bead2d(:,1,2)); % beads tracked in both camera views
        
        for bd =  find(ind_tracked)'             % loop through the tracked beads
            [bead3dG(1:3,bd),err(bd)] = triangulate(bead2d(bd,1:2,1),bead2d(bd,1:2,2),epi_geo.P1',epi_geo.P2');
        end
        
        if (sum(ind_tracked) >= 2 ) && ismember(bonesRef{bn},bonesOptimize) % if this frame can be optimized
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
        
        pos3D.(bonesRef{bn})(:,:,fr)= bead3dG;
        
        
    end
    
    if (strcmp(optimizeOpts.Type,'optimize')) && ~isempty(intersect(bonesRef{bn},bonesOptimize)) % we can optimize, and it's in the selected group of bones to optimize
        
        for i = 1:size(optimizeFrames,1)
            
            fr = optimizeFrames(i,1); % for all the frames that can be optimized
            bd = optimizeFrames(i,2);
            camflag = optimizeFrames(i,3);
            epi_error = optimizeFrames(i,4);
            bead2d = optimizeFrames(i,5:7);
            bead3dG = pos3D.(bonesRef{bn})(:,:,fr);
            
            if camflag == 1             % MISSING bead in CAM1
                
                cam = 2;
                [epi_planeG,epi_lineG] = pointEpipolarGeometry(epi_geo,bead2d',cam);
                epi_poleG = epi_geo.E1G;
                
                
            elseif camflag == 2         % MISSING bead in CAM2
                cam = 1;
                [epi_planeG,epi_lineG] = pointEpipolarGeometry(epi_geo,bead2d',cam);
                epi_poleG = epi_geo.E2G;
                
            end
            
            frsTracked = find(~isnan(squeeze(pos3D.(bonesRef{bn})(1,bd,:))));
            frsDiff = abs(frsTracked-fr);
            [M,I] = min(frsDiff);
            frRef = frsTracked(I);
            rvI = M * 0.25;
            refPtG =  pos3D.(bonesRef{bn})(:,:,frRef);
            
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
            
            N =  [100;100;100] ;
            sigma = [1,1,0.2];
            limits = [beadB-1.5,beadB+1.5] ;
            [beadB,rotval_save(bd,fr,:)] = beadLocationProbability(limits ,N, epi_geoB, ptsB,refPtB, sigma,epi_error/2,'refine',rvI,0);
            bead3dG(:,bd) = transformPoints(T_GtB,beadB',-1);
            
            
            % evaluate the found point
            pt = closestPointonPlanealongVector(beadB, epi_geoB.plane_norm,epi_geoB.epipole, epi_geoB.plane_norm);
            epiplaneDist = norm(pt-beadB);
            fprintf('Frame %i\nBone: %s.\n',fr,bonesRef{bn})
            fprintf('Distance from epipolar plane is %0.3f mm \n',epiplaneDist)
            fprintf('Distance difference from CT beads is L13 %0.3f mm,L23 %0.3f mm . \n\n',L_CT(2)- norm(beadB-ptsB(:,1)),L_CT(3) - norm(beadB-ptsB(:,2)))
            %                     savedat.(bone_list(bn,:))(bd).results(ind,:) = [fr , epiplaneDist,abs(L_CT(2)- norm(beadB-ptsB(:,1))),abs(L_CT(3)- norm(beadB-ptsB(:,2))) ];
            %                     ind = ind+1;
            pos3D.(bonesRef{bn})(:,bd,fr)= bead3dG(:,bd);
            
            if strcmp(filterOpts.Type(1:2),'2D') % if the filter type is 2D, reproject, filter and then triangulate
                
                im_pt1 = epi_geo.P1 * [pos3D.(bonesRef{bn})(:,bd,fr);1];
                im_pt1 = im_pt1/im_pt1(3);
                im_pt2 = epi_geo.P2 * [pos3D.(bonesRef{bn})(:,bd,fr);1];
                im_pt2 = im_pt2/im_pt2(3);
                pos2D.(bonesRef{bn})(bd).cam1(fr,:) = im_pt1(1:2)';
                pos2D.(bonesRef{bn})(bd).cam2(fr,:) = im_pt2(1:2)';
%             end
                
        
        
        
%         if strcmp(filterOpts.Type(1:2),'2D') && ~isempty(optimizeFrames) % if the filter type is 2D, reproject, filter and then triangulate
            % refilter the 2D points
            pos2D.(bonesRef{bn})(bd).cam1 = adaptiveLowPassButterworth(pos2D.(bonesRef{bn})(bd).cam1',filterOpts.Fc,filterOpts.Fs,1)';
            pos2D.(bonesRef{bn})(bd).cam2 = adaptiveLowPassButterworth(pos2D.(bonesRef{bn})(bd).cam2',filterOpts.Fc,filterOpts.Fs,1)';
            
            for fr = 1:nfr
                
                bead2d = [];
                    bead2d(:,1) = [pos2D.(bonesRef{bn})(bd).cam1(fr,:) 1];
                    bead2d(:,2) = [pos2D.(bonesRef{bn})(bd).cam2(fr,:) 1];
              
                
                % loop through and triangulate
                bead3dG = nan(3,1);
                if ~any(isnan(bead2d),'all') % if there are no nan values - 
                    
                    [bead3dG(1:3),err] = triangulate(bead2d(1:2,1)',bead2d(1:2,2)',epi_geo.P1',epi_geo.P2');
                end
                
                pos3D.(bonesRef{bn})(:,bd,fr)= bead3dG;
            end
            end
            
        end
                
                
                
    end
    
    
end

if ~isfield(saveOpts,'animDir') ||  ~isfield(saveOpts,'animDir') 
   
        matchXMAtoCT(pos3D,beadCT,subjDir,trialName,filterOpts);
 
else
        matchXMAtoCT(pos3D,beadCT,subjDir,trialName,filterOpts,saveOpts.animDir,saveOpts.ivDir);
end

outFile = [];
% write a new file so that it doesn't have to be recomputed every time 
if  strcmp(filterOpts.Type,'none') && strcmp(optimizeOpts.Type,'optimize') 
    outFile = [in2DFile(1:end-4) '_optimized.csv'];
elseif ~strcmp(filterOpts.Type,'none') && strcmp(optimizeOpts.Type,'optimize') 
    outFile = [in2DFile(1:end-4) '_optimized_filtered' filterOpts.Type '.csv'];
elseif ~strcmp(filterOpts.Type,'none')
    outFile = [in2DFile(1:end-4) '_filtered' filterOpts.Type '.csv'];    
end
if ~isempty(outFile)
    writeXMA2Dfrom3DPoints(in2DFile,outFile,epi_geo,pos3D);
end

