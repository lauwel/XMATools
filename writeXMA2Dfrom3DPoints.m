function writeXMA2Dfrom3DPoints(inFile,outFile,epi_geo,pos3D)
% inFile is an example file from XMA Lab




if verLessThan('matlab','9.6')
    error('loadXMA2dPoints.m requires MATLAB 9.6 or later (R2019a).')
end



numdata = readmatrix(inFile,'Delimiter',',','Range','A2');
[nFr,n] = size(numdata);

% rng = sprintf('1:%i',n);
headers = readcell(inFile,'Delimiter',',','Range',[1 1 1 n]);
% headers_char = char(headers);
% bead_list = char(unique(cellstr(headers_char(:,1:4)))); % gets the bead list in char form
% nBeads = length(bead_list);
% bone_list = char(unique(cellstr(headers_char(:,1:3)))); % gets the bone list in char form
% nBones = length(bone_list);

% convert the XMA file

for i = 1:4:n % for each of the headers by bone/bead X/Y
for fr = 1:nFr 
   % go through each header, attribute all the data from the beads, then fill
   % in
   strs = strsplit(headers{i},'_');
    boneName = strs{1}(1:3);
    beadNum = str2double(strs{1}(4));
%     camNum = str2double(strs{2}(4));
%     
%     switch strs{3}
%         case 'X'
%             pos_ind = 1;
%         case 'Y'
%             pos_ind = 2;
%     end
    
   im_pt1 = epi_geo.P1 * [pos3D.(boneName)(:,beadNum,fr);1];
   im_pt1 = im_pt1/im_pt1(3);
   im_pt2 = epi_geo.P2 * [pos3D.(boneName)(:,beadNum,fr);1];
   im_pt2 = im_pt2/im_pt2(3);
   
   numdata_new(fr,i:i+3) = [im_pt1(1:2)', im_pt2(1:2)'];
  
   
%         numdata(:,i) = pos2D.(boneName)(beadNum).(sprintf('cam%i',camNum))(:,pos_ind) ;

end
end

numdata_cell = num2cell(numdata_new);
writecell([headers;numdata_cell],outFile)

% 
% varargout{1} = pos2D;
% varargout{2} = nBones;
% varargout{3} = nBeads;