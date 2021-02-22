
close all
clear
clc

% Compute all the transforms for the trials



% FILE INFORMATION

subjDir = 'E:\SOL001B\'; % subject directory
trialName = 'T0088_SOL001_nrun_rfs_minimal'; % the trial name
Fs = 250; % hz
fc = [10,20]; 
beadFile = 'E:\SOL001B\Models\bead_positions.txt'; % the bead locations from the CT scan, you get this after running sphereFitToBeads.m
camDir = fullfile(subjDir,'Calibration\Set 3\MayaCam2\'); % where your Mayacam files are located
 
%% IF YOU NEED TO INTERPOLATE: EXPORT A 2D DISTORTED POINTS FILE FROM XMALAB
in2DFile = fullfile(subjDir,trialName,'XMA_csv',[trialName 'DISTORTED.csv']); % this is the undistorted 2D points file that you exported from XMA lab

interpolateXMA2DPoints(in2DFile,[in2DFile(1:end-4) '_interp.csv'],'fib')
% you'll have to load it back into XMALab and switch the tracked marker
% points

%% FIRST PASS : OPTIMIZE TO FIND MISSING BEADS & PRODUCE UNFILTERED
% TRANSFORMS AND ANIMATIONS

in2DFile = fullfile(subjDir,trialName,'XMA_csv',[trialName 'UNDISTORTED.csv']); % this is the undistorted 2D points file that you exported from XMA lab

filterOpts.Type = 'none'; % this tells the processer to not filter 
optimizeOpts.Type ='none'; %'optimize'; % this will find positions of missing beads if there is enough other information 
optimizeOpts.Bones = {'tib','tal','cal','nav','cmm','cub','mt1'};
saveOpts.ivDir = fullfile(subjDir,'Models','IV','3Aligned Reduced',filesep);
saveOpts.animDir = fullfile(subjDir,trialName,'POS',filesep);
processBeadData(in2DFile,beadFile,camDir,subjDir,trialName,filterOpts,optimizeOpts,saveOpts);


%% SECOND PASS : FILTER:
in2DFile = fullfile(subjDir,trialName,'XMA_csv',[trialName 'UNDISTORTED.csv']); % either the same file as the one above, or if you optimized in the previous step, then use that 2D points file
filterOpts.Type = '2D';
filterOpts.Fs = Fs;
    % fc = [w1 w2] This is the filtering frequency. 
% It will filter your data somewhere between w1 and w1+w2 based on how much 
% change there is in your data.
filterOpts.Fc = fc;
% optimizeOpts.Type = 'none';
% saveOpts = [];

processBeadData(in2DFile,beadFile,camDir,subjDir,trialName,filterOpts,optimizeOpts,saveOpts);

%% If you need to add autoscoped frames to the bone transforms file

boneName = 'cal';
transformFile = 'E:\SOL001B\T0088_SOL001_nrun_rfs_minimal\BoneTransforms\T0088_SOL001_nrun_rfs_minimal_transforms2DFILT.mat';
autoscoperFile = 'E:\SOL001B\T0088_SOL001_nrun_rfs_minimal\Autoscoper\Autoscoped\T0088_SOL001_nrun_rfs_minimal_cal.tra'; % every frame has to be saved

load(transformFile);
auto_dat = dlmread(autoscoperFile);

ifr = ~isnan(auto_dat(:,1));
T4x4 = convertRotation(auto_dat,'autoscoper','4x4xn');

T.(boneName)(:,:,ifr) = T4x4(:,:,ifr);

save(transformFile,'T')

