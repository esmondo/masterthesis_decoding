
close all
clear all
clc

%%

addpath('/export/data/reichert/toolbox/MD_tools')
addpath('/export/data/reichert/toolbox/Elekta')
addpath('/export/data/reichert/toolbox/fieldtrip')
addpath('/export/data/duerschm/allscripts')

addpath(genpath('/export/data/esmondo'))
savepath = '/export/data/esmondo/any_results/';

%%
if ~exist('topoplot','file')
    addpath /export/data/reichert/toolbox/ft4topoplot/
end

cfgfileMEG = '/export/data/reichert/capLayout/Neuromag_helmet.mat';
cfgfileEEG = '/export/data/reichert/capLayout/EEG_29_cartesianSorted.mat';

cfgMEG = topoprepare(cfgfileMEG);
cfgEEG = topoprepare(cfgfileEEG);
% ft_layoutplot(cfg);

%% Input dataset

whichdata = input('Select "2" for the old dataset (only MEG) | Select "4" for the new dataset (EEG & MEG) = ');

%%
switch whichdata
    
    case 1
    %% DATASET: srate = 2000Hz
    
    cd /export/data/wienke/data/Esmondo
    myFolder = '/export/data/wienke/data/Esmondo/';

    if ~isdir(myFolder)
        errorMessage = sprintf('Error: The following folder does not exists:\n%s', myFolder);
        uiwait(warndlg(errorMessage));
     return;
    end

    filePattern = fullfile(myFolder,'MW*meg.mat');
    filelist = dir(filePattern);
    endIter = length(filelist);

    
    case 2
%% DATASET: MAGNETOMETER DATA

%     cd /export/data/database/MEG/mindwandering_n2pc_sss
    myFolder = '/export/data/database/MEG/mindwandering_n2pc_sss';


    if ~isdir(myFolder)
        errorMessage = sprintf('Error: The following folder does not exists:\n%s', myFolder);
        uiwait(warndlg(errorMessage));
        return;
    end

    filePattern = fullfile(myFolder,'MW*megSSS.mat');
    filelist = dir(filePattern);
    endIter = length(filelist);
    
    
    case 3
 %% DATASET: GRADIOMETER DATA

    cd /export/data/database/MEG/mindwandering_n2pc_sss
    myFolder = '/export/data/database/MEG/mindwandering_n2pc_sss';


    if ~isdir(myFolder)
        errorMessage = sprintf('Error: The following folder does not exists:\n%s', myFolder);
        uiwait(warndlg(errorMessage));
        return;
    end

    filePattern = fullfile(myFolder,'MW*gradSSS.mat');
    filelist = dir(filePattern);
    endIter = length(filelist);

    
    case 4
%% NEW DATASET: MEG+EEG+EOG DATA

    cd '/export/data/esmondo/any_results/';
    myFolder = '/export/data/esmondo/any_results/';


    if ~isdir(myFolder)
        errorMessage = sprintf('Error: The following folder does not exists:\n%s', myFolder);
        uiwait(warndlg(errorMessage));
        return;
    end

    filePattern = fullfile(myFolder,'MW*newSSS.mat');
    filelist = dir(filePattern);
    endIter = length(filelist);
    
end    


%%
for i = 1:endIter

close all  

%% load data

% Load data option (1)
matFilename = fullfile(myFolder, filelist(i).name);
load(matFilename,'meg');
if whichdata == 4
    load(matFilename,'eeg');
end
load(matFilename,'eog');

fprintf('Analyzing %s , ', filelist(i).name);


%% define data

  
data_meg = double(meg.data);  % MEG  
    
if whichdata == 4          
    eeg.data = eeg.data( 1:29,:,: )-repmat( eeg.data( 30,:,: )./2,[29 1 1] ); % EEG re-referencing
    data_eeg = double(eeg.data);
end

data_eog = double(eog.data); % EOG

time = meg.time;

targetSide = meg.side;
subjectRes = meg.response;
targetPos = meg.position;
focusRate = meg.fokus;

%% Sensors Labelling

 % All Sensors
 chanALL = meg.header.label;
 
 % MEG
 chan_meg = meg.header.label(meg.channels); 
    % chanidx = strmatch('MEG',meg.header.label); % the same with meg.channels
    % chanidx = chanidx( 1:3:end );

 if whichdata == 4
     % EEG
    chan_eeg = meg.header.label(eeg.channels(1:29));
    chanEEGname = {'Fz','Cz','Pz','Oz','Iz','Fp1','Fp2','F3','F4','F7'...
                'F8','T7','T8','C3','C4','P3','P4','PO9','PO10',...
                'P7','P8','FC1','FC2','CP1','CP2','PO3','PO4','PO7','PO8'}';
 
 
 % Index of PO7 and PO8 in EEG
 PO7idx = find(strcmp(chanEEGname,'PO7')); % located in the left hemisphere
 PO8idx = find(strcmp(chanEEGname,'PO8')); % located in the right hemisphere
 
 % Index of posterior channels in EEG
 pterioreeg = {'P3','P4','PO9','PO10','Pz','Oz'...
                'P7','P8','PO3','PO4','PO7','PO8'};
 pterioreegidx = find(ismember(chanEEGname, pterioreeg));
 
 end

 % EOG
 chan_eog = meg.header.label(eog.channels); 

%% AOIs Index: Occipitotemporal region in MEG

[cL,cR,cAll] = getAOI_occitemp;
 
%% Size of necessary data

% indexing essential parameters
[x1_meg, x2_meg, x3_meg(i)] = size(data_meg);

if whichdata == 4
    [x1_eeg, x2_eeg, x3_eeg(i)] = size(data_eeg);
end

%% correct response
[~,corID] = corrResponse(meg.side,meg.response);

%% label
[cl_label,~,~] = hemisphaereTeilen(meg.position);
cl_left = (cl_label == 0);
cl_right = (cl_label == 1);

Xmeg_left_sum = sum(cl_left);
Xmeg_right_sum = sum(cl_right);

% Xmeg_L = data_meg(:,:,cl_left);
% Xmeg_R = data_meg(:,:,cl_right);
Xmeg_L = data_meg(cAll,:,cl_left);
Xmeg_R = data_meg(cAll,:,cl_right);

%%
% ERP_L = squeeze(mean(Xmeg_L, 3));
% ERP_R = squeeze(mean(Xmeg_R, 3));
% % fprintf('Size of ERP_attended: [%d %d]\n', size(ERP_L))
% % fprintf('Size of ERP_unattended: [%d %d]\n', size(ERP_R))
% 
% figure
% plot(time, ERP_L, 'r-')
% hold all, grid on
% plot(time, ERP_R, 'b-')
% 
% title('ERP at all channels (red=L, blue=R)')
% xlabel('Time [s]'), ylabel('Amplitude [T]')

% find the time points corresponding to 0-0.8 s
itime = find(time >= 0 & time <= 0.8);  

% Extract the mean activity in the interval as features
Xmeg2 = squeeze(mean(data_meg(cAll,itime,:),2));

%%
X = data_meg(cAll,itime,corID);
X = reshape(permute(X,[3,2,1]),[size(X,3),size(X,1)*size(X,2)]);
Y = cl_label(corID);

%%
cfg = [];
cfg.k = 10;
cfg.repeat = 1;
cfg.classifier = 'svm';
perf = mv_classify(cfg, X, Y);

end



