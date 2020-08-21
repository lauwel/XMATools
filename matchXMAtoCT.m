function varargout = matchXMAtoCT(varargin)
% matchXMAtoCT(pos3D,beadCT,subj_dir,trial_name)
% matchXMAtoCT(pos3D,beadCT,subj_dir,trial_name,filterOpts)
% matchXMAtoCT(pos3D,beadCT,subj_dir,trial_name,filterOpts,anim_dir,ivDir)
% matchXMAtoCT(pos3D,beadCT,subj_dir,trial_name,filterOpts,anim_dir,ivDir,autoscoper_dir)
% [T] = matchXMAtoCT(pos3D,beadCT,subj_dir,trial_name,...)
%
% This function will match the 3D bead positions to the beads in CT
% space, and will compute the transforms for each frame from CT to x-ray
% space. The transforms can be output as a variable 'T', but are also saved
% in a .mat file. The autoscoper .tra files for each bone, as well as a
% combined "allBone" trial is saved. The bones are also animated in Wrist
% Visualizer.
% ----------------------------------INPUTS--------------------------------
% pos3D         =       The 3D computed bead position structure from
%                       projectXMA2DPoints.m. Alternatively, the location
%                       to a 3D XMA points file can be used.
% beadCT        =       The 3D positions of the beads in CT space, created
%                       in projectXMA2DPoints. Alternatively, the location
%                       of the bead_positions.txt file can be used, and the
%                       structure will be recomputed in this function.
% subj_dir      =       The subject's folder location. It will be used to
%                       locate Model files, and trial directory.
% trial_name    =       The name of the trial, and subsequently how all the
%                       files saved in this code are named.
% filterOpts    =       [OPTIONAL] Uses adaptiveLowPassButterworth.m to low
%                       pass filter the 3D points. Is passed as a structure
%                       with the fields .Type, .Fs and .Fc (sampling and cutoff
%                       frequency respectively). Note: The adaptive
%                       filtering cutoff frequency must be in a 1x2 array.
%                       If this argument is [], no filtering will be done.
% anim_dir      =       [OPTIONAL] Specify a specific directory under which
%                       to save the animations. Must be specified WITH
%                       ivDir.
% ivDir         =       The iv file directory for the bones being animated.
% autoscoper_dir =      [OPTIONAL] Specify a specific directory for the
%                       autoscoper.tra files.
%
% --------------------------------OUTPUTS---------------------------------
%
% T             =       [OPTIONAL] The structure containing the transforms
%                       for all of the bones.
% ------------------------------------------------------------------------
% L. Welte Sept 2019
%



% assign all the variables
if isstruct(varargin{1})
    pos3D = varargin{1};
elseif isfile(varargin{1})
    
    [beadPos,names] = xlsread([varargin{1}]); % load the file
    beadPos = beadPos(:,2:end)'; % get rid of the first column and transpose
    names = names(1,2:end)'; % get rid of first column and transpose
    
    % save all the beads into the appropriate structure location
    for b = 1:size(beadPos,1)
        
        beadname = names{b}(2:4);
        beadnum = str2double(names{b}(5));
        tempBead = beadPos(b,:);
        
        switch names{b}(7)
            case 'X'
                pos3D.(beadname)(1,beadnum,:) = tempBead';
            case 'Y'
                pos3D.(beadname)(2,beadnum,:) = tempBead';
            case 'Z'
                pos3D.(beadname)(3,beadnum,:) = tempBead';
                
            otherwise
                error('Mislabelled bead.')
        end
    end
else
    error('Incorrect format for first input argument to matchXMAtoCT.')
end

if isstruct(varargin{2})
    beadCT = varargin{2};
elseif isfile(varargin{2})
    
    datBeadLocs = csvread(varargin{2},1,1);
    dat = readtable(varargin{2});
    beadStructTemp = table2struct(dat);
    nbeads = length(beadStructTemp);
    for b = 1:nbeads
        
        fieldname = beadStructTemp(b).bead;
        beadnum = str2double(fieldname(4));
        beadCT.(fieldname(1:3))(:,beadnum) = datBeadLocs(b,:);
        
    end
end

subjDir = varargin{3};
trial_name = varargin{4};

% subj_name = 'SOL001B';
trial_dir = [subjDir trial_name '\'];

switch nargin
    case 5
        if ~isempty(varargin{5})
            filterOpts  = varargin{5};
            filtType = filterOpts.Type;
            if strcmp(filtType,'none')
                filtType = 'Un';
                filter_flag = 0;
            else
                fc = filterOpts.Fc;
                fs = filterOpts.Fs;
                filter_flag = 1;
            end
        else
            filter_flag = 0;
            filtType = '';
        end
        anim_dir = strcat(trial_dir,['POS\' filtType 'filtered\']);
        autoscoper_dir = [trial_dir 'Autoscoper\'];
        ivDir = [subjDir 'Models/IV/'];
    case 7
        if ~isempty(varargin{5})
            filterOpts  = varargin{5};
            filtType = filterOpts.Type;
            if strcmp(filtType,'none')
                filtType = 'Un';
                filter_flag = 0;
            else
                fc = filterOpts.Fc;
                fs = filterOpts.Fs;
                filter_flag = 1;
            end
        else
            filter_flag = 0;
            filtType = 'Un';
        end
        anim_dir = [varargin{6}, filtType, 'filtered\'];
        autoscoper_dir = [trial_dir 'Autoscoper\'];
        ivDir = varargin{7};
        % arguments 6 and 7 are passed together
    case 8
        if ~isempty(varargin{5})
            filterOpts  = varargin{5};
            filtType = filterOpts.Type;
            if strcmp(filtType,'none')
                filtType = 'Un';
                filter_flag = 0;
            else
                fc = filterOpts.Fc;
                fs = filterOpts.Fs;
                filter_flag = 1;
            end
        else
            filter_flag = 0;
            filtType = 'Un';
        end
        anim_dir = [varargin{6}, filtType, 'filtered\'];
        ivDir = varargin{7};
        autoscoper_dir = varargin{8};
    otherwise
        ivDir = [subjDir 'Models/IV/'];
        filter_flag = 0;
        filtType = 'Un';
        anim_dir = strcat(trial_dir,['POS\' filtType 'filtered\']);
        autoscoper_dir = [trial_dir 'Autoscoper\'];
        
end



% subj_dir = 'P:\Data\2019-05-02 SOL001_Visit2\'; % if saving animation to
% server


bone_list = char(fields(pos3D));
bonesCell = fields(pos3D);
nBones = size(bone_list,1);
nFr = size(pos3D.(bone_list(1,:)),3);


%% Filter the points if necessary

if filter_flag == 1
    beadPos = [];
    
    for bn = 1:nBones
        nBeadsBone = size(pos3D.(bone_list(bn,:)),2);
        for bd = 1:nBeadsBone
            for dim = 1:3
                beadPos = squeeze(pos3D.(bone_list(bn,:))(dim,bd,:));
                
                if length(find(~isnan(beadPos)))<6
                    continue
                elseif contains(filtType,'3D')
                    filteredBead = adaptiveLowPassButterworth(beadPos',fc,fs);
                    pos3D.(bone_list(bn,:))(dim,bd,:) = filteredBead;
                else
                    pos3D.(bone_list(bn,:))(dim,bd,:) = beadPos';
                end
            end
        end
    end
end



%% match the bead positions to the CT to get the transform
T_save = [];

% tracked_frames = nan(nbones,nFr);
for bn = 1:nBones % for every bone, compute the transform between the beads and xma points
    bonebeadsCT = beadCT.(bone_list(bn,:));
    ind = 1;
    for fr = 1:nFr
        % get the transform from CT to x-ray space
        T.(bone_list(bn,:))(:,:,fr) = eye(4,4);
        Tanim.(bone_list(bn,:))(:,:,fr) = eye(4,4);
        bonebeadsXMA = pos3D.(bone_list(bn,:))(:,:,fr);
        ind_nonan = ~isnan(bonebeadsXMA(1,:));
        
        if sum(ind_nonan) >= 3 % if there are at least 3 beads with data
            
            [R,d,rms] = soder(bonebeadsCT(:,ind_nonan)', bonebeadsXMA(:,ind_nonan)');
            
            T.(bone_list(bn,:))(1:3,1:3,fr) = R;
            T.(bone_list(bn,:))(1:3,4,fr) = d;
            
            Tanim.(bone_list(bn,:))(1:3,1:3,fr) = R;
            Tanim.(bone_list(bn,:))(1:3,4,fr) = d;
            
            rms_save(bn,fr) = rms;
            
            tracked_frames(bn,ind) = fr;
            ind = ind+1;
        else
            T.(bone_list(bn,:))(:,:,fr) = nan(4,4);
            Tanim.(bone_list(bn,:))(:,:,fr) = ones(4,4);
        end
        
    end
    
    % write the file with the bone transforms
    
    Tstacked.(bone_list(bn,:)) = convertRotation(T.(bone_list(bn,:)),'4x4xn','autoscoper');
    
    TstackedR.(bone_list(bn,:)) = convertRotation(T.(bone_list(bn,:)),'4x4xn','autoscoperRows');
    %
    
    if exist(autoscoper_dir,'dir') == 0
        mkdir(autoscoper_dir);
    end
    
    
    auto_file = [autoscoper_dir trial_name '_' (bone_list(bn,:)) '_' filtType 'filt.tra'];
    
    dlmwrite(auto_file,Tstacked.(bone_list(bn,:)))
    
    T_save = [T_save TstackedR.(bone_list(bn,:))];
end



% Save the ALL bones file - in rows
auto_file = [autoscoper_dir 'allBones_' trial_name '_rows_' filtType 'filt.tra'];

dlmwrite(auto_file,T_save)




boneT_dir = fullfile(trial_dir,'BoneTransforms',filesep);
if exist(boneT_dir,'dir') == 0
    mkdir(boneT_dir);
end

% save the transforms in a nice MAT file

save([boneT_dir,trial_name '_transforms' upper(filtType) 'FILT.mat'],'T')


fprintf('Autoscoper files written to %s.\n',autoscoper_dir)
fprintf('Bone transforms are written to %s.\n',boneT_dir)

if nargout == 1
    varargout{1} = T;
end



figure;
h = bar(rms_save');
nbars = length(h);
cmap = colormap('jet');
indd = round(linspace(1,64,nbars));
cmap = cmap(indd,:);
for ii = 1:nbars
    h(ii).FaceColor = cmap(ii,:);
end
legend(bone_list)
ylabel('RMS of XMA to CT reconstruction')
xlabel('Frame');
%% Animate the trial



rigidiv_folder = fullfile(anim_dir,'rigidiv',filesep);

if exist(anim_dir,'dir')==0;     mkdir(anim_dir);     end
if exist(rigidiv_folder,'dir')==0;  mkdir(rigidiv_folder);  end

ind0 = tracked_frames == 0;
tracked_frames(ind0) = NaN;
first_fr = nanmin(nanmin(tracked_frames(:,:)));
if first_fr == 0
    first_fr = 1;
end
end_fr   = nanmax(nanmax(tracked_frames(:,:)));

ivFileListTemp = dir([ivDir '*.iv']);
ivFileList = {ivFileListTemp(:).name};
for bn = 1:nBones
    
    ivstring = createInventorHeader();
    % make the linked iv file
    boneFile = contains(ivFileList,bone_list(bn,:));
    ind_bone = find(boneFile);
    % verify the bone files found
    if length(ind_bone) > 1
        warning('There is more than one %s file in the iv directory. Using %s. ',bone_list(bn,:),ivFileList{ind_bone(1)} )
        
    elseif isempty(ind_bone)
          error('There is no %s file in the iv directory.',bone_list(bn,:))
    end
    newBoneFile = ivFileList{boneFile};
        
    ivstring = [ivstring createInventorLink([ivDir newBoneFile ],eye(3,3),zeros(3,1),[0.7 0.7 0.7],0.5)];
    
    fid = fopen(fullfile(rigidiv_folder,[bone_list(bn,:) '.iv']),'w');
    fprintf(fid,ivstring);
    fclose(fid);
end


for bn = 1:nBones
    write_RTp(bone_list(bn,:) , Tanim.(bone_list(bn,:))(:,:,first_fr:end_fr) , anim_dir)
end

pos_text = write_pos(bonesCell,anim_dir,trial_name);

filename = fullfile(anim_dir, [trial_name '.pos']);



fid = fopen(filename,'w'); % open the file to write
fprintf(fid,pos_text);
fclose(fid);

fprintf('Animations are saved in %s.\n',anim_dir)

%
