function varargout = loadXMA2dPoints(inFile)
% [pos2D] = loadXMA2dPoints(inFile);
% [pos2D,nBones,nBeads] = loadXMA2dPoints(inFile);
% This function requires R2019a. 

if verLessThan('matlab','9.6')
    error('loadXMA2dPoints.m requires MATLAB 9.6 or later (R2019a).')
end

% load the 2d data from XMA lab
numdata = readmatrix(inFile,'Delimiter',',','Range','A2');
[~,n] = size(numdata);

% rng = sprintf('1:%i',n);
headers = readcell(inFile,'Delimiter',',','Range',[1 1 1 n]);
headers_char = char(headers);
bead_list = char(unique(cellstr(headers_char(:,1:4)))); % gets the bead list in char form
nBeads = length(bead_list);
bone_list = char(unique(cellstr(headers_char(:,1:3)))); % gets the bone list in char form
nBones = length(bone_list);

% convert the XMA file

for i = 1:n % for each of the headers
   % go through each header, attribute all the data from the beads, then fill
   % in
   strs = strsplit(headers{i},'_');
    boneName = strs{1}(1:3);
    beadNum = str2double(strs{1}(4));
    camNum = str2double(strs{2}(4));
    
    switch strs{3}
        case 'X'
            pos_ind = 1;
        case 'Y'
            pos_ind = 2;
    end
  
        pos2D.(boneName)(beadNum).(sprintf('cam%i',camNum))(:,pos_ind) = numdata(:,i);

end

varargout{1} = pos2D;
varargout{2} = nBones;
varargout{3} = nBeads;

