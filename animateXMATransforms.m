function animateXMATransforms(varargin)
% animateXMATransforms(subj_dir,trial_name)
% animateXMATransforms(subj_dir,trial_name,fc)
% animateXMATransforms(subj_dir,trial_name,fc,anim_dir)
% animateXMATransforms(subj_dir,trial_name,fc,anim_dir,autoscoper_dir)
% 
%  The transforms are savedin a .mat file.
%  The autoscoper .tra files for each bone, as well as a combined "allBone"
%  trial is saved. The bones are also animated in Wrist Visualizer. 
% 
% ----------------------------------INPUTS--------------------------------
% subj_dir      =       The subject's folder location. It will be used to
%                       locate Model files, and trial directory.
% trial_name    =       The name of the trial, and subsequently how all the
%                       files saved in this code are named.
% fc            =       [OPTIONAL] The adaptive filtering cutoff frequency 
%                       in a 1x2 array. (Uses adaptiveLowPassButterworth.m) 
%                       If this argument is [], no filtering will be done.
% anim_dir      =       [OPTIONAL] Specify a specific directory under which
%                       to save the animations.
% autoscoper_dir =      [OPTIONAL] Specify a specific directory for the
%                       autoscoper.tra files.

% --------------------------------OUTPUTS---------------------------------
% None
% ------------------------------------------------------------------------
% L. Welte Oct 2019



subj_dir = varargin{1};
trial_name = varargin{2};

inFile = fullfile(subj_dir,trial_name,'XMA_csv',[trial_name 'TRANSFORMS.csv']);
subj_name = 'SOL001B';

ivdir = [subj_dir 'Models\IV\'];
trial_dir = [subj_dir trial_name '\'];

switch nargin
    case 3 
        if ~isempty(varargin{3})
            fc = varargin{3};
            filter_flag = 1;
        else
            filter_flag = 0;
        end
        anim_dir = strcat(trial_dir,'POS\Filtered\');
        autoscoper_dir = [trial_dir 'Autoscoper\'];
    case 4        
        if ~isempty(varargin{3})
            fc = varargin{3};
            filter_flag = 1;
        else
            filter_flag = 0;
        end
        anim_dir = varargin{4};
        autoscoper_dir = [trial_dir 'Autoscoper\'];
    case 5
        if ~isempty(varargin{3})
            fc = varargin{3};
            filter_flag = 1;
        else
            filter_flag = 0;
        end
        anim_dir = varargin{4};
        autoscoper_dir = varargin{5};
    otherwise
        filter_flag = 0;
        anim_dir = strcat(trial_dir,'POS\Unfiltered\');
        autoscoper_dir = [trial_dir 'Autoscoper\'];

end

Tauto_all = [];

if verLessThan('matlab','9.6')
    error('loadXMA2dPoints.m requires MATLAB 9.6 or later (R2019a).')
end

trial_dir = fullfile(subj_dir,trial_name,filesep);


numdata = readmatrix(inFile,'Delimiter',',','Range','A2');
[nfr,h] = size(numdata);
headers = readcell(inFile,'Delimiter',',','Range',[1 1 1 h]);
% headers_char = char(headers);
for hd = 1:16:h
    strs = strsplit(headers{hd},'_');
    T.(strs{1}) = convertRotation(numdata(:,hd:hd+15),'autoscoper','4x4xn');
end

bone_list = fields(T);
bone_list = sort(bone_list);
nBones = length(bone_list);
% tracked_frames = nan(nbones,nFr);
for bn = 1:nBones % for every bone, compute the transform between the beads and xma points
   
    % write the file with the bone transforms
    
    TQ.(bone_list{bn}) = convertRotation(T.(bone_list{bn}),'4x4xn','quaternion');
    tracked_frames{bn} = find(~isnan(TQ.(bone_list{bn})(:,1)));
     if filter_flag == 1
        TQfilt = adaptiveLowPassButterworth(TQ.(bone_list{bn})',fc,250,0)';
%         norm(TQ.(bone_list{bn})(140,1:4))
%         TQfilt(tracked_frames{bn},1:4) = unit(TQfilt(tracked_frames{bn},1:4));
        for fr = 1:nfr
            TQfilt(fr,1:4) = unit(TQfilt(fr,1:4));
        end
        auto_file = [autoscoper_dir trial_name '_' (bone_list{bn}) '_filt.tra'];
            
     else
        TQfilt =  TQ.(bone_list{bn});
        auto_file = [autoscoper_dir trial_name '_' (bone_list{bn}) '_unfilt.tra'];
     end
     
     
    Tauto.(bone_list{bn}) = convertRotation(TQfilt,'quaternion','autoscoper');
    TautoR.(bone_list{bn}) = convertRotation(TQfilt,'quaternion','autoscoperRows');
    Tanim.(bone_list{bn}) = convertRotation(TQfilt,'quaternion','4x4xn');
    
    ind_nan = isnan(Tanim.(bone_list{bn}));
    Tanim.(bone_list{bn})(ind_nan) = 1;
    if exist(autoscoper_dir,'dir') == 0
        mkdir(autoscoper_dir);
    end
    
    
    dlmwrite(auto_file,Tauto.(bone_list{bn}))
    
    Tauto_all = [Tauto_all TautoR.(bone_list{bn})];
end



% Save the ALL bones file - in rows
if filter_flag == 1
    auto_file = [autoscoper_dir 'allBones_' trial_name '_rows_filt.tra'];
else
    auto_file = [autoscoper_dir 'allBones_' trial_name '_rows_unfilt.tra'];
end
dlmwrite(auto_file,Tauto_all)



boneT_dir = fullfile(trial_dir,'BoneTransforms',filesep);
if exist(boneT_dir,'dir') == 0
        mkdir(boneT_dir);
end
% save the transforms in a nice MAT file
save([boneT_dir,trial_name '_transforms.mat'],'T')

fprintf('Autoscoper files written to %s.\n',autoscoper_dir)
fprintf('Bone transforms are written to %s.\n',boneT_dir)


%% Animate the trial



rigidiv_folder = fullfile(anim_dir,'rigidiv',filesep);

if exist(anim_dir,'dir')==0;     mkdir(anim_dir);     end
if exist(rigidiv_folder,'dir')==0;  mkdir(rigidiv_folder);  end



first_fr = min(cellfun(@min,tracked_frames)); 
end_fr   = max(cellfun(@max,tracked_frames));

for bn = 1:nBones
    
    ivstring = createInventorHeader();
    % make the linked iv file
    ivstring = [ivstring createInventorLink([ivdir subj_name '_'  bone_list{bn} '_aligned.iv'],eye(3,3),zeros(3,1),[0.7 0.7 0.7],0.5)];
    
    fid = fopen(fullfile(rigidiv_folder,[bone_list{bn} '.iv']),'w');
    fprintf(fid,ivstring);
    fclose(fid);
end


for bn = 1:nBones
    write_RTp(bone_list{bn} , Tanim.(bone_list{bn})(:,:,first_fr:end_fr) , anim_dir)
end

pos_text = write_pos(bone_list,anim_dir,trial_name);

    filename = fullfile(anim_dir, [trial_name '.pos']);



fid = fopen(filename,'w'); % open the file to write
fprintf(fid,pos_text);
fclose(fid);

fprintf('Animations are saved in %s.\n',anim_dir)