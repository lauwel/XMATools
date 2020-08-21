
close all
clear
clc

% Compute all the transforms for the trials



% FILE INFORMATION

subjDir = 'E:\SOL001B\'; % subject directory
trialName = 'T0026_SOL001_nrun_ffs_barefoot'; % the trial name


%% IF YOU NEED TO INTERPOLATE: EXPORT A 2D DISTORTED POINTS FILE FROM XMALAB
in2DFile = fullfile(subjDir,trialName,'XMA_csv',[trialName 'DISTORTED.csv']); % this is the undistorted 2D points file that you exported from XMA lab

interpolateXMA2DPoints(in2DFile,[in2DFile(1:end-4) '_interp.csv'],'nav')
% you'll have to load it back into XMALab and switch the tracked marker
% points

%% FIRST PASS : OPTIMIZE TO FIND MISSING BEADS & PRODUCE UNFILTERED
% TRANSFORMS AND ANIMATIONS

in2DFile = fullfile(subjDir,trialName,'XMA_csv',[trialName 'UNDISTORTED.csv']); % this is the undistorted 2D points file that you exported from XMA lab
beadFile = 'E:\SOL001B\Models\bead_positions.txt'; % the bead locations from the CT scan, you get this after running sphereFitToBeads.m
camDir = fullfile(subjDir,'Calibration\Set 1\MayaCam2\'); % where your Mayacam files are located

filterOpts.Type = 'none'; % this tells the processer to not filter 

optimizeOpts.Type ='optimize'; %'optimize'; % this will find positions of missing beads if there is enough other information 
optimizeOpts.Bones = {'tib','tal','cal','nav','cmm','cub','mt1','mt5'};
saveOpts.ivDir = fullfile(subjDir,'Models','IV','3Aligned Reduced',filesep);
saveOpts.animDir = fullfile(subjDir,trialName,'POS',filesep);
processBeadData(in2DFile,beadFile,camDir,subjDir,trialName,filterOpts,optimizeOpts,saveOpts);


%% SECOND PASS : FILTER:
in2DFile = fullfile(subjDir,trialName,'XMA_csv',[trialName 'UNDISTORTED.csv']); % either the same file as the one above, or if you optimized in the previous step, then use that 2D points file
filterOpts.Type = '2D';
filterOpts.Fs = 250;
    fc = [10,30]; % fc = [w1 w2] This is the filtering frequency. 
% It will filter your data somewhere between w1 and w1+w2 based on how much 
% change there is in your data.
filterOpts.Fc = fc;
% optimizeOpts.Type = 'none';
% saveOpts = [];

processBeadData(in2DFile,beadFile,camDir,subjDir,trialName,filterOpts,optimizeOpts,saveOpts);