function interpolateXMA2DPoints(inFile,outFile,bone,bead)

% This function takes 2D points from XMA to interpolate, to then be taken
% back into XMA. If you select specific beads or bones, only those
% beads/bones will be written back to the file.
% Export the points as : 
% - combined file
% - header row
% - distorted co-ordinates
% - start co-ordinates at 0
% - y axis down
% - no leading columns

% --------------INPUTS-----------------------------------------------------
% inFile  =     path + file of export from XMA
% outfile =     optional argument to specify the save file location; 
%               default is the same filename with "_interp" at the end
% bone =        [optional] 3 letter string of one bone bead to interpolate (i.e. 'tib')
% bead =        [optional] (double) number of the bead to look at (i.e. 3)
% -------------OUTPUTS-----------------------------------------------------
% none currently
% -------------------------------------------------------------------------

% L. Welte 6/2019

% -------------------------------------------------------------------------


if exist(inFile,'file') ~= 2 % locate the file
    error('File not found.')
end

if exist('outFile','var') == 0% if outFile is not specified
    C = strsplit(inFile,'.');
    fext = ['.' C{end}];
    outFile = [C{1} '_interp' fext];
elseif isempty(outFile)
    C = strsplit(inFile,'.');
    fext = ['.' C{end}];
    outFile = [C{1} '_interp' fext];
end

[pos2D,nBones,nBeads] = loadXMA2dPoints(inFile);
bone_list = fields(pos2D);
nfr = size(pos2D.(bone_list{1})(1).cam1,1);

if ~exist('bone','var')
    bone_ind = 1:nBones;
else
    bone_ind = find(strcmp(bone_list,bone));
end



numdata = readmatrix(inFile,'Delimiter',',','Range','A2');
[nfr,n] = size(numdata);
headers = readcell(inFile,'Delimiter',',','Range',[1 1 1 n]);

new_headers = [];
ind = 1;
for bn = bone_ind
    nBeadsBone = size(pos2D.(bone_list{bn}),2);
    if ~exist('bead','var')
        bead_ind = 1:nBeadsBone;
    else
        bead_ind = bead;
    end
    
    for bd = bead_ind
        
        for c = 1:2 % each camera
            xycoords = pos2D.(bone_list{bn})(bd).(sprintf('cam%i',c));
            xy_interp = normaliseNaN(xycoords,1,nfr);
%             pos2D.(bone_list{bn})(bd).(sprintf('cam%i',c)) = xy_interp;
            
            new_headers{ind} = [bone_list{bn} sprintf('%i_cam%i_X',bd,c)];
            
            new_headers{ind+1} = [bone_list{bn} sprintf('%i_cam%i_Y',bd,c)];
%             data_rep = find(strcmp(headers, new_headers{ind}));
            new_data(:,ind:ind+1) = xy_interp;
%             numdata(:,data_rep:data_rep+1) = xy_interp;
            ind = ind+2;
        end
    end
end
% 
% 
% 
% for bd = 1:n
%     if sum(isnan(numdata(:,bd))) >= nfr-5
%         numdata_interp(:,bd) = numdata(:,bd);
%     else
%         numdata_interp(:,bd) = normaliseNaN(numdata(:,bd),1,nfr);
%     end
% end
numdata_cell = num2cell(new_data);
writecell([new_headers;numdata_cell],outFile)
fprintf('Interpolated file written to %s.\n',outFile)


