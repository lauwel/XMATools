


close all
clear 
clc
xma_dir = 'E:\SOL001_VISIT2\Testing\XMA\';
frs_test = 151:170;

camFolder = fullfile('E:\SOL001_VISIT2\Calibration\Set 1\Mayacam2\');
subj_dir = 'E:\SOL001_VISIT2\TESTING\';
[pos3Dsave(1),~] = projectXMA2DPoints(fullfile(xma_dir,'Data2D.csv'),camFolder,subj_dir,'Data2D',[],'optimize');
[pos3Dsave(2).cal,calres] = projectXMA2DPoints(fullfile(xma_dir,'Data2D_cal1missingcam1.csv'),camFolder,subj_dir,'Data2D_cal1missingcam1',[20 60],'optimize','cal');
[pos3Dsave(3).cub,cubres] = projectXMA2DPoints(fullfile(xma_dir,'CUB2.csv'),camFolder,subj_dir,'CUB2',[],'optimize','cub');
[pos3Dsave(4).tal1,tal1res] = projectXMA2DPoints(fullfile(xma_dir,'Data2D_tal2missingcam1.csv'),camFolder,subj_dir,'Data2D_tal2missingcam1',[],'optimize','tal');
[pos3Dsave(5).tal2,tal2res] = projectXMA2DPoints(fullfile(xma_dir,'Data2D_tal2missingcam2.csv'),camFolder,subj_dir,'Data2D_tal2missingcam2',[],'optimize','tal');
[pos3Dsave(6).mt1,~] = projectXMA2DPoints(fullfile(xma_dir,'Data2D_mt13missingcam2.csv'),camFolder,subj_dir,'Data2D_mt13missingcam2',[],'optimize','mt1');
[pos3Dsave(7).mt5,~] = projectXMA2DPoints(fullfile(xma_dir,'Data2D_mt53missingcam1.csv'),camFolder,subj_dir,'Data2D_mt53missingcam1',[],'optimize','mt5');


%%
% Look at calcaneus
close all
figure; hold on;
goldstd = squeeze(pos3Dsave(1).cal(:,1,frs_test))';
test = squeeze(pos3Dsave(2).cal.cal(:,1,frs_test))';
plot(goldstd,'-')
plot(test,'--')
legend('X-Hand-tracked','Y-Hand-tracked','Z-Hand-tracked',...
    'X-missing cam1','Y-missing cam1','Z-missing cam1')
xlabel('Frame')
ylabel('Position [mm]')
title('Calcaneus missing bead 1')

xyzStr = {'x','y','z'};
fprintf('Calcaneus results:\n')
% figure;hold on;

for bd = 1
    
    for dims = 1:3
        rm = rms(goldstd(:,dims)-test(:,dims));
        mm = max(abs(goldstd(:,dims)-test(:,dims)));
        fprintf('Bead %i : The %s RMS error for cam 1 missing bead 3 in MATLAB is %0.3fmm, max %0.3fmm.\n',bd,xyzStr{dims},rm,mm)
%         plot(calres.cal(bd).results(:,2:4),goldstd(:,dims)-test(:,dims),'.-')
    end
end
%%

figure;hold on;
goldstd = squeeze(pos3Dsave(1).tal(:,2,frs_test))';
test1 = squeeze(pos3Dsave(4).tal1.tal(:,2,frs_test))';
test2 = squeeze(pos3Dsave(5).tal2.tal(:,2,frs_test))';
plot(goldstd,'-')
plot(test1,'.-')
plot(test2,'--')
legend('X-Hand-tracked','Y-Hand-tracked','Z-Hand-tracked',...
    'X-missing cam1','Y-missing cam1','Z-missing cam1',...
    'X-missing cam2','Y-missing cam2','Z-missing cam2')
xlabel('Frame')
ylabel('Position [mm]')
title('Talus missing bead 2')

% figure; hold on;
fprintf('Talus results:\n')
for bd = 2
    
    for dims = 1:3
        rm = rms(goldstd(:,dims)-test1(:,dims));
        mm = max(abs(goldstd(:,dims)-test1(:,dims)));
        fprintf('Bead %i : The %s RMS error for cam 1 missing bead 2 in MATLAB is %0.3fmm, max %0.3fmm.\n',bd,xyzStr{dims},rm,mm)
        
%          plot(tal1res.tal(bd).results(:,2:4),goldstd(:,dims)-test1(:,dims),'.-')
        
        rm = rms(goldstd(:,dims)-test2(:,dims));
        mm = max(abs(goldstd(:,dims)-test2(:,dims)));
        fprintf('Bead %i : The %s RMS error for cam 2 missing bead 2 in MATLAB is %0.3fmm, max %0.3fmm.\n',bd,xyzStr{dims},rm,mm)
    
%          plot(tal2res.tal(bd).results(:,2:4),goldstd(:,dims)-test1(:,dims),'.-')
    end
end
%% cub
figure; hold on;
bd = 1
goldstd = squeeze(pos3Dsave(1).cub(:,bd,frs_test))';
test = squeeze(pos3Dsave(3).cub.cub(:,bd,frs_test))';
plot(goldstd,'-')
plot(test,'--')

legend('X-Hand-tracked','Y-Hand-tracked','Z-Hand-tracked',...
    'X-missing cam1','Y-missing cam1','Z-missing cam1')
xlabel('Frame')
ylabel('Position [mm]')
title('Cuboid missing bead 3')
%

% figure; hold on;
fprintf('Cuboid results:\n')
% for bd = 1
    
    for dims = 1:3
        rm = rms(goldstd(:,dims)-test(:,dims));
        mm = max(abs(goldstd(:,dims)-test(:,dims)));
        fprintf('Bead %i : The %s RMS error for cam 1 missing bead %i in MATLAB is %0.3fmm, max %0.3fmm.\n',bd,xyzStr{dims},bd,rm,mm)
        
%          plot(cubres.cub(1).results(:,2:4),goldstd(:,dims)-test(:,dims),'.-')
    end
% end

%% mt1
figure; hold on;
bd = 3;
goldstd = squeeze(pos3Dsave(1).mt1(:,bd,frs_test))';
test = squeeze(pos3Dsave(6).mt1.mt1(:,bd,frs_test))';
plot(goldstd,'-')
plot(test,'--')

legend('X-Hand-tracked','Y-Hand-tracked','Z-Hand-tracked',...
    'X-missing cam2','Y-missing cam2','Z-missing cam2')
xlabel('Frame')
ylabel('Position [mm]')
title('MT1 missing bead 3')
%

% figure; hold on;
fprintf('Mt1 results:\n')


for dims = 1:3
    rm = rms(goldstd(:,dims)-test(:,dims));
    mm = max(abs(goldstd(:,dims)-test(:,dims)));
    fprintf('Bead %i : The %s RMS error for cam 2 missing bead %i in MATLAB is %0.3fmm, max %0.3fmm.\n',bd,xyzStr{dims},bd,rm,mm)
    
    %          plot(cubres.cub(1).results(:,2:4),goldstd(:,dims)-test(:,dims),'.-')
end
%%
figure; hold on;
bd = 3;
goldstd = squeeze(pos3Dsave(1).mt5(:,bd,frs_test))';
test = squeeze(pos3Dsave(7).mt5.mt5(:,bd,frs_test))';
plot(goldstd,'-')
plot(test,'--')

legend('X-Hand-tracked','Y-Hand-tracked','Z-Hand-tracked',...
    'X-missing cam1','Y-missing cam1','Z-missing cam1')
xlabel('Frame')
ylabel('Position [mm]')
title('MT5 missing bead 3')
%

% figure; hold on;
fprintf('Mt5 results:\n')


for dims = 1:3
    rm = rms(goldstd(:,dims)-test(:,dims));
    mm = max(abs(goldstd(:,dims)-test(:,dims)));
    fprintf('Bead %i : The %s RMS error for cam 1 missing bead %i in MATLAB is %0.3fmm, max %0.3fmm.\n',bd,xyzStr{dims},bd,rm,mm)
    
    %          plot(cubres.cub(1).results(:,2:4),goldstd(:,dims)-test(:,dims),'.-')
end

