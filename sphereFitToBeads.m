close all
clear
clc

ivDir = 'E:\SOL001_VISIT1\Models\IV\Beads\';

list_files = dir([ivDir 'SOL*.iv']);

ivstring = createInventorHeader();
cmap = colormap('jet');
close all;
bone_list = {'tib','mt1','tal','cub','nav','mt5','cal','cmm','fib','ph1'};
cm = cmap(1:6:end,:);

    ind = 1;
for i = 1:length(list_files)
    surfaceFile = [ivDir list_files(i).name];
    [pts cnt] = read_vrml_fast(surfaceFile);
    cnt = cnt + 1;

    [Center,Radius] = sphereFit(pts);
    
    if contains(surfaceFile,'out') % if it's a bead not in a bone
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
  writetable(data_vals,'E:\SOL001_VISIT1\Models\bead_positions.txt')
  
fid = fopen([ivDir 'SphereFitBeads.iv'],'w');
fprintf(fid,ivstring);
fclose(fid);
fprintf(['File written to: ' ivDir 'SphereFitBeads.iv\n'])
