close all
clear
clc

% You will  need to download the mathworks function sphereFit: https://www.mathworks.com/matlabcentral/fileexchange/34129-sphere-fit-least-squared

% This code will take the rough exported iv files from Mimics, calculate a
% sphere fit, and then save the centroids of the beads to a text file. 
% It will also make a WristViz file in the specified iv directory with all
% the beads and their names that can be used to see all the beads in 3D. 

ivDir = 'E:\SOL001B\Models\IV\Beads\'; % the location of the exported bead iv files
modelsDir = 'E:\SOL001B\Models\'; % where the bead positions text file will be written to

list_files = dir([ivDir 'SOL*.iv']); % all the beads are iv files with the subject SOL letters in it

ivstring = createInventorHeader();

bone_list = {'tib','mt1','tal','cub','nav','mt5','cal','cmm','fib','ph1'}; % the bone codes for each bone with beads


%%
cmap = colormap('jet');
cm = cmap(1:6:end,:);

    ind = 1;
for i = 1:length(list_files) % for each bead
    surfaceFile = [ivDir list_files(i).name];
    [pts, cnt] = read_vrml_fast(surfaceFile); % load the bead
    cnt = cnt + 1;

    [Center,Radius] = sphereFit(pts); % fit a sphere
    
    if contains(surfaceFile,'out') % if it's a bead not in a bone -> some of the foot beads are not within a bone, designated with "out" in the file name
        col = [1 1 1];
    else
        
        bead{ind,1} = surfaceFile(end-6:end-3);
        centroid{ind,1} = Center;
        
       ci =  find(contains(bone_list,bead{ind,1}(1:3)));
       col = cm(ci,:);
        ind = ind+1;
        
       
    end
    
    ivstring = [ivstring createInventorSphere(Center,0.75,col,0.2)];
    ivstring = [ivstring createInventorText(bead{ind-1},4,Center + [0.5 0.5 0.5],col,0)];
   
end

  data_vals =  table(bead,centroid);
  
  writetable(data_vals,[ modelsDir 'bead_positions.txt']) % where to write the beads centres to
  
fid = fopen([ivDir 'SphereFitBeads.iv'],'w'); % write the full bead file
fprintf(fid,ivstring);
fclose(fid);
fprintf(['File written to: ' ivDir 'SphereFitBeads.iv\n'])
