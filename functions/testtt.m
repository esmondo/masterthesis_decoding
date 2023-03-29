%%
close all
clear all
clc

%%
% cd /export/data/wienke/data/Esmondo
addpath('/export/data/wienke/data/Esmondo')
addpath('/export/data/esmondo')
addpath('/export/data/esmondo/functions')
addpath('/export/data/reichert/toolbox/MD_tools')
addpath('/export/data/reichert/toolbox/Elekta')
addpath('/export/data/database/MEG/mindwandering_n2pc_sss')
addpath('/export/data/duerschm/allscripts')
savePath = '/export/data/esmondo/';

%%
if ~exist('topoplot','file')
    addpath /export/data/reichert/toolbox/ft4topoplot/
end

cfgfile = '/export/data/reichert/capLayout/Neuromag_helmet.mat';
cfg = topoprepare(cfgfile);

%% Input dataset

whichdata = input('Select case = ');

%%
switch whichdata
    
    case 1
    %% DATASET: 1Hz-30Hz
    
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
%% DATASET: ARTIFACT REJECTED AND FILTERED BETWEEN 1Hz-200Hz

    cd /export/data/wienke/data/motor
    myFolder = '/export/data/wienke/data/motor';


    if ~isdir(myFolder)
        errorMessage = sprintf('Error: The following folder does not exists:\n%s', myFolder);
        uiwait(warndlg(errorMessage));
        return;
    end

    filePattern = fullfile(myFolder,'MW*meg.mat');
    filelist = dir(filePattern);
    endIter = length(filelist);

    case 3
%% DATASET: MAGNETOMETER DATA

    cd /export/data/database/MEG/mindwandering_n2pc_sss
    myFolder = '/export/data/database/MEG/mindwandering_n2pc_sss';


    if ~isdir(myFolder)
        errorMessage = sprintf('Error: The following folder does not exists:\n%s', myFolder);
        uiwait(warndlg(errorMessage));
        return;
    end

    filePattern = fullfile(myFolder,'MW*megSSS.mat');
    filelist = dir(filePattern);
    endIter = length(filelist);
    
    
    case 4
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
    
end    

%% Execute all subjects

for i = 1:endIter

close all  

%% load data

% Load data option (1)
matFilename = fullfile(myFolder, filelist(i).name);
load(matFilename,'meg');
load(matFilename,'eog');

fprintf('Analyzing %s , ', filelist(i).name);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 % PRE-PROCESSING AND ARTIFACT REJECTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% define data

data = double(meg.data);
data_decode = double(meg.data);
dataEOG = double(eog.data);

%% Time interval

time = meg.time;
timeEOG = eog.time;

%% MEG channels detail

chan = meg.header.label(meg.channels);
sensors = neuromagGetSensorAOIs;

%% all sensors

chansLeft = unique([sensors.frontal.left, sensors.parietal.left,...
    sensors.temporal.left, sensors.occipital.left])';

chansRight = unique([sensors.frontal.right, sensors.parietal.right,...
    sensors.temporal.right, sensors.occipital.right])';

chansAll = unique([sensors.frontal.left, sensors.parietal.left,...
    sensors.temporal.left, sensors.occipital.left,...
    sensors.frontal.right, sensors.parietal.right,...
    sensors.temporal.right, sensors.occipital.right])';

%% channel of interest (occipital)

chansOleft = unique([sensors.occipital.left])';

chansOright = unique([sensors.occipital.right])';

chansOall = unique([sensors.occipital.left, sensors.occipital.right])';

%% channel of interest (occipito-temporal)

chansOTleft = unique([sensors.temporal.left, sensors.occipital.left])';

chansOTright = unique([sensors.temporal.right, sensors.occipital.right])';

chansOTall = unique([sensors.temporal.left, sensors.temporal.right,...
    sensors.occipital.left, sensors.occipital.right])';

%% channel of interest (occipito-temporal-parietal)
chansOTPleft = unique([sensors.temporal.left, sensors.occipital.left, sensors.parietal.left])';

chansOTPright = unique([sensors.temporal.right, sensors.occipital.right, sensors.parietal.right])';

chansOTPall = unique([sensors.temporal.left, sensors.temporal.right,...
    sensors.occipital.left, sensors.occipital.right,...
    sensors.parietal.left,sensors.parietal.right])';

%% Size of necessary data

% indexing essential parameters
[idxChan, idxPoints, idxTrials(i)] = size(data);
[idxChanEOG, idxPointsEOG, idxTrialsEOG(i)] = size(dataEOG);


%%

targetSide = meg.side;
subjectRes = meg.response;
targetPos = meg.position;
focusRate = meg.fokus;


%% Filtering and Artifact Rejection


    if whichdata == 1
        
        % High Pass filter
        data = highpass(3.0,meg.srate,data);
    
        % Artficat Detection
        artFind = findArtifacts(data);
        data(:,:,artFind) = [];
    
        % removing artifacts in side of targets's tilt
        targetSide(artFind) = [];

        % removing artifacts in subject's response
        subjectRes(artFind) = [];

        % removing artifacts in position where the target appeared
        targetPos(artFind) = [];
        
        % removing artifacts in participant's focus ratings
        focusRate(artFind) = [];
        
    elseif whichdata == 3
        
        % Band Pass filter
        data = bandpass(1,200,meg.srate,data);
        dataEOG = bandpass(1,30,eog.srate,dataEOG);
             
        % Eye Movement Rejection
          data = rejectComponent(data,squeeze(dataEOG(1,:,:))); % HEOG
         
%         data = rejectComponent(data,squeeze(dataEOG(2,:,:))); % VEOG

%         % Exclude Blinks
%         eyeMove_1 = zeros(1,idxTrialsEOG(i));
%         eyeMove_2 = zeros(1,idxTrialsEOG(i));
%         for em = 1:idxTrialsEOG(i)
%             eyeMove_1(:,em) = max(abs(dataEOG(1,find(timeEOG>=-0.2 & timeEOG<=0.8),em)))';
%             eyeMove_2(:,em) = max(abs(dataEOG(2,find(timeEOG>=-0.2 & timeEOG<=0.8),em)))';
%         end
%         
%         % Indexes where blink occured (will be removed)
%         thBlinks1 = find(eyeMove_1 > mean(eyeMove_1)+(std(eyeMove_1)*2));
%         thBlinks1 = find(eyeMove_1 >  1e-4);
%         thBlinks2 = find(eyeMove_2 > mean(eyeMove_2)+(std(eyeMove_2)*2));
%         thBlinks2 = find(eyeMove_2 >  1e-4);
%         
%         % Combination of blink indexes VEOG and HEOG (use this if needed)
%         exclBlinks = cat(2,thBlinks1,thBlinks2); 
%         
%         data(:,:,thBlinks2) = [];
%         targetSide(thBlinks2) = [];
%         subjectRes(thBlinks2) = [];
%         targetPos(thBlinks2) = [];
%         focusRate(thBlinks2) = [];
         
%         % Artficat Detection
%         artFind = findArtifacts(data);
%         
%         % Artifact Rejection
%         data(:,:,artFind) = [];
%         targetSide(artFind) = []; % removing artifacts in side of targets's tilt
%         subjectRes(artFind) = []; % removing artifacts in subject's response
%         targetPos(artFind) = []; % removing artifacts in position where the target appeared
%         focusRate(artFind) = []; % removing artifacts in participant's focus ratings
                 
    end

%% Downsample data

    if whichdata == 1
        Fs = 100;
    elseif whichdata == 2
        Fs = 400;
    else
        Fs = meg.srate;
        FsEOG = eog.srate;
    end

[timeDS,dataDS] = dwsample(idxChan,Fs,meg.srate,data,time);

%% baseline correction

data = baseline(dataDS,timeDS);

%% time range -200 to 800ms

timeIdx = find(timeDS>=-0.2 & timeDS<=0.8);
timeNew = timeDS(timeIdx);

%% correct respose from subject

[corrResp, idxCorrResp] = corrResponse(targetSide,subjectRes);

% correct target presentation's side
corrSide = targetSide(idxCorrResp);

% correct focus ratings
corrFocus = focusRate(idxCorrResp);

% correct position
corrTargPos = targetPos(idxCorrResp);

%% index of trials when the focus question occured
idxMental = find(corrFocus > 0);

idxON = find(ismember(corrFocus,[4,5]));
idxMID = find(ismember(corrFocus,[3]));
idxOFF = find(ismember(corrFocus,[1,2]));


mentalStr = cell(size(corrFocus(idxMental)));
mentalStr(corrFocus(idxMental) == 1) = {'OFF'};
mentalStr(corrFocus(idxMental) == 2) = {'OFF'};
mentalStr(corrFocus(idxMental) == 3) = {'MID'};
mentalStr(corrFocus(idxMental) == 4) = {'ON'};
mentalStr(corrFocus(idxMental) == 5) = {'ON'};
% mentalStr = mentalStr';



%% Final and clean data

% data indicated by mental states
dataON = data(:,timeIdx,idxON);
dataMID = data(:,timeIdx,idxMID);
dataOFF = data(:,timeIdx,idxOFF);
dataMental = data(:,timeIdx,idxMental);

% the main data
data = data(:,timeIdx,idxCorrResp);



%% number of cleared trial:

numGutTrial(i) = size(data,3);
numBadTrial(i) = idxTrials(i)-numGutTrial(i);
numONtrial(i) = length(idxON);
numMIDtrial(i) = length(idxMID);
numOFFtrial(i) = length(idxOFF);
numMentalTrial(i) = length(idxMental);

%% CHANGE PARAMETERS HERE :)
% channel for revealing N2pc and ecoding N2pc 
% (can be changed, depends on which channel group you want to measure)

cL = chansOTleft;
cR = chansOTright;
cAll = chansOTall;

navg = 10;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       % STATISTICAL ANALYSIS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% data separation

% separate target by LVF and RVF 
[corrPos,targetLVF, targetRVF, idxTargetLVF, idxTargetRVF] = ...
    hemisphaereTeilen(data,corrTargPos);


numLeftTrial(i) = length(idxTargetLVF);
numRightTrial(i) = length(idxTargetRVF);


%% t-test

% compare targets in LVF vs targets in RVF
[h,pVal,ci,stats] = ttest2(targetLVF,targetRVF,'Dim',3);

% the t-values
tValues(:,:,i) = stats.tstat;

% the p-values
pValues(:,:,i) = pVal;

% t-values exceeding p-threshold (significantly different)
sigDiff(:,:,i) = tValues(:,:,i).*double(pValues(:,:,i)<0.05);


%% Max_positive and Max_negative sensors selection
% From "distributions of t-values", occipito-temporal sensors showing 
% Max(+) and Max(-) t-values in the time range 200-300ms were selected on
% each hemisphere


% time range 200 msec - 300 msec post stimulus
timeN2pcIdx = find(timeNew>=0.2 & timeNew<=0.3);
timeN2pcVal = timeNew(timeN2pcIdx);

[effluxLeft(i),influxLeft(i),effluxRight(i),influxRight(i)] = ...
    maxPmaxN(tValues(:,:,i),cL,cR,timeN2pcIdx);

effluxL(i) = cL(effluxLeft(i));
influxL(i) = cL(influxLeft(i));
effluxR(i) = cR(effluxRight(i));
influxR(i) = cR(influxRight(i));

maxchan = [effluxL(i),influxL(i),effluxR(i),influxR(i)];

%% combination of efflux channel & influx channel (substraction)


[avgLeftHem_LVF(:,:,i), avgLeftHem_RVF(:,:,i), avgRightHem_LVF(:,:,i), avgRightHem_RVF(:,:,i)] = effluxMinInflux...
    (data,effluxL(i),influxL(i),effluxR(i),influxR(i),idxTargetLVF,idxTargetRVF);


%% smoothering signals for plot

avgLeftHem_LVF(:,:,i) = bandpass(0.5,30,meg.srate,avgLeftHem_LVF(:,:,i));
avgLeftHem_RVF(:,:,i) = bandpass(0.5,30,meg.srate,avgLeftHem_RVF(:,:,i));
N2pcLH_sub(:,:,i) = avgLeftHem_LVF(:,:,i)-avgLeftHem_RVF(:,:,i);
avgRightHem_LVF(:,:,i) = bandpass(0.5,30,meg.srate,avgRightHem_LVF(:,:,i));
avgRightHem_RVF(:,:,i) = bandpass(0.5,30,meg.srate,avgRightHem_RVF(:,:,i));
N2pcRH_sub(:,:,i) = avgRightHem_LVF(:,:,i)-avgRightHem_RVF(:,:,i);

% % LH
% leftFig = figure('Name','Left Hemisphere');
% hold on
% plot(timeNew,avgLeftHem_LVF(:,:,i),'b');
% plot(timeNew,avgLeftHem_RVF(:,:,i),'r--');
% plot(timeNew,N2pcLH_sub(:,:,i),'k','LineWidth',2)
% hold off
% xlabel('time')
% ylabel('amplitude')
% title('Left Hemisphere')
% legend('LVF Target','RVF Target','LVF-RVF')
% ylim([-2.5e-13 2.5e-13])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');
% line([0.18 0.18],ylim,'Color','k','LineStyle','--');
% line([0.32 0.32],ylim,'Color','k','LineStyle','--');
% leftFigName = sprintf('leftHem_sub%d',i); 
% saveas(leftFig, fullfile(savePath,leftFigName),'png');
% 
% % RH
% rightFig = figure('Name','Right Hemisphere');
% hold on
% plot(timeNew,avgRightHem_LVF(:,:,i),'b');
% plot(timeNew,avgRightHem_RVF(:,:,i),'r--');
% plot(timeNew,N2pcRH_sub(:,:,i),'k','LineWidth',2)
% hold off
% xlabel('time')
% ylabel('amplitude')
% title('Right Hemisphere')
% legend('LVF Target','RVF Target','LVF-RVF')
% ylim([-2.5e-13 2.5e-13])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');
% line([0.18 0.18],ylim,'Color','k','LineStyle','--');
% line([0.32 0.32],ylim,'Color','k','LineStyle','--');
% rightFigName = sprintf('rightHem_sub%d',i);  
% saveas(rightFig, fullfile(savePath,rightFigName),'png');

%% check the N2pc amplitude
[N2pcLHpeak_val(i),N2pcLHpeak_id(i)] = max(N2pcLH_sub(:,timeN2pcIdx,i));
[N2pcRHpeak_val(i),N2pcRHpeak_id(i)] = max(N2pcRH_sub(:,timeN2pcIdx,i));

%% time point where the N2pc peaked
N2pcLHpeak_pos(i) = timeN2pcVal(N2pcLHpeak_id(i));
N2pcRHpeak_pos(i) = timeN2pcVal(N2pcRHpeak_id(i));

%% distribution plot for several time point


tPoint = 175:5:325;
% 
% for tp = 1:length(tPoint)
%     
%     if tp == 12
%         continue
%     else
%         tidx = timeNew>tPoint(tp)*1e-3 & timeNew<tPoint(tp+1)*1e-3;
%     end
%     
%     megframe = double(mean(tValues(:,tidx,i),2));  
%     topoAll = figure('Name','Topoplot');
%     hold on
%     toposhow(cfg, megframe);
%     clr = colorbar;
%     clr.Limits = [-2 2];
%     hold off
% 
%     titleFig = 'Distribution of N2pc at %d ms';
%     topoFilename = 'toposub%d_%d';
%     A1 = tPoint(tp)+5;
%     A2 = i;
%     
%     title(sprintf(titleFig,A1))
%     topoAllName = sprintf(topoFilename,A2,A1'); 
%     saveas(topoAll, fullfile(savePath,topoAllName),'png');
% 
% end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           % DECODING %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % ---------
% % time of interest: 0-800ms
% timeCutIdx = find(timeNew>=0 & timeNew(end));
% timeCut = timeNew(timeCutIdx);
% data = data(:,timeCutIdx,:);
% dataEOG = dataEOG(:,timeCutIdx,:);
% 
% % lowpass filter 12.5 Hz
% data = lowpass(12.5,Fs,data);
% 
% % Downsample 50 Hz
% FsCls = 50;
% [timeCls,dataCls] = dwsample(idxChan,FsCls,Fs,data,timeCut);
% 
% 
% %% Classification
% 
% % divTrial = fix(size(dataCls,3)/400);
% % remTrial = mod(size(dataCls,3),400);
% % trial_for_decoding = 1:divTrial:(size(dataCls,3)-remTrial);
% 
% X = double(meg.data(cAll,:,:)); 
% % X = dataCls(cAll,:,:) .* 1e12; % convert to pT
% % X = dataCls(cAll,:,trial_for_decoding) .* 1e12; % convert to pT
% 
% [~,corID] = corrResponse(meg.side,meg.response);
% Y = classPos(meg.position);
% Y = Y(corID);
% 
% % Y = corrPos; % initial decoding
% 
% % Y = corrPos(trial_for_decoding); % use the balanced trial
% 
% % for it = 1:1000 % performing permutation test
% %     Y = Y(randperm(length(Y)));
% % end
% 
% 
% N = 10; % N-fold
% 
% if navg <= 1
%     nreps = 1; 
% else
%     nreps = 5;
% end
% 
% 
% %%
% [predictions,predictionsSTF,alg1,alg2] = decodeLR(X,Y,N,navg,nreps,meg,eog,timeNew);
% 
% 
% %% precision & recall 
% 
% [precision_c1(i,:),recall_c1(i,:),precision_c2(i,:),recall_c2(i,:),acconly(i,:)] = precisionNrecall(Y,predictions);
% [precision_c1STF(i,:),recall_c1STF(i,:),precision_c2STF(i,:),recall_c2STF(i,:),accstf(i,:)] = precisionNrecall(Y,predictionsSTF); % STF
% 
% 
% %% accuracies
% 
% accuracy(i) = accCheck(Y,predictions,nreps);
% accuracySTF(i) = accCheck(Y,predictionsSTF,nreps);
% % accuracyL(i) = acc(Yl,predictionsL,nreps);
% % accuracyL_STF(i) = acc(Yl,predictionsL,nreps);
% % accuracyR(i) = acc(Yr,predictionsR,nreps);
% % accuracyR_STF(i) = acc(Yr,predictionsR,nreps);
% 
% %% 
% fprintf('Accuracy = %3.2f%% , ', accuracy(i));
% fprintf('Accuracy with spatial filter =  %3.2f%% \n', accuracySTF(i));


end



% %% used channel
% usedChan = cat(1,effluxL,influxL,effluxR,influxR);
% usedChanName = chan(usedChan);
% 
% %%  accuracy, precision, and recall of each subject | 
% accOnly = mean(acconly, 2)'*100;
% accSTF = mean(accstf, 2)'*100;
% 
% recall_1 = mean(recall_c1,2); %  Left Presentation Relevant Trials / TNR / specificity
% recallSTF_1 = mean(recall_c1STF,2);
% 
% recall_2 = mean(recall_c2,2); %  Right Presentation Relevant Trials / TPR / sensitivity
% recallSTF_2 = mean(recall_c2STF,2);
% 
% %% accuracy comparison
% accComp = cat(1,accuracy,accuracySTF);
% 
% %% Grand Average of Accuracy
% 
% accuracyMean = mean(accuracy);
% accuracySD = std(accuracy);
% 
% accuracyMeanSTF = mean(accuracySTF);
% accuracySD_STF = std(accuracySTF);


%% Grand Average N2pc on each hemisphere


% left hem
GA_avgLeftHem_LVF = mean(avgLeftHem_LVF,3);
GA_avgLeftHem_RVF = mean(avgLeftHem_RVF,3);
GA_N2pcLH = GA_avgLeftHem_LVF-GA_avgLeftHem_RVF;

% right hem
GA_avgRightHem_LVF = mean(avgRightHem_LVF,3);
GA_avgRightHem_RVF = mean(avgRightHem_RVF,3);
GA_N2pcRH = GA_avgRightHem_LVF-GA_avgRightHem_RVF;

GAleftFig = figure('Name','GA N2pc in Left Hemisphere');
hold on
plot(timeNew,GA_avgLeftHem_LVF,'b');
plot(timeNew,GA_avgLeftHem_RVF,'r--');
plot(timeNew,GA_N2pcLH,'k','LineWidth',2);
hold off
xlabel('time')
ylabel('amplitude')
title('GA N2pc in Left Hemisphere')
legend('LVF Target','RVF Target','LVF-RVF')
% yL = get(gca,'ylim');
ylim([-1e-13 1.5e-13])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');
line([0.18 0.18],ylim,'Color','k','LineStyle','--');
line([0.32 0.32],ylim,'Color','k','LineStyle','--');
GAleftFigName = sprintf('GAleftHem');
saveas(GAleftFig, fullfile(savePath,GAleftFigName),'png');

GArightFig = figure('Name','GA N2pc Right Hemisphere');
hold on
plot(timeNew,GA_avgRightHem_LVF,'b');
plot(timeNew,GA_avgRightHem_RVF,'r--');
plot(timeNew,GA_N2pcRH,'k','LineWidth',2);
hold off
xlabel('time')
ylabel('amplitude')
title('GA N2pc Right Hemisphere')
legend('LVF Target','RVF Target','LVF-RVF')
% yL = get(gca,'ylim');
ylim([-1e-13 1.5e-13])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');
line([0.18 0.18],ylim,'Color','k','LineStyle','--');
line([0.32 0.32],ylim,'Color','k','LineStyle','--');
GArightFigName = sprintf('GArightHem'); 
saveas(GArightFig, fullfile(savePath,GArightFigName),'png');
 
%% check the GA N2pc amplitude
[GA_N2pcLHpeak_val,GA_N2pcLHpeak_id] = max(GA_N2pcLH(timeN2pcIdx));
[GA_N2pcRHpeak_val,GA_N2pcRHpeak_id] = max(GA_N2pcRH(timeN2pcIdx));

%% time point where the GA N2pc peaked
GA_N2pcLHpeak_pos = timeN2pcVal(GA_N2pcLHpeak_id);
GA_N2pcRHpeakk_pos = timeN2pcVal(GA_N2pcRHpeak_id);

%% t-values exceeding p-threshold

tValExceedP = squeeze(mean(sigDiff,3));
tValExceedP_maxc = squeeze(mean(sigDiff(maxchan,:,:),3));


tvpvFig = figure('Name','t-values exceeding p-values < 0.05');
imagesc(tValExceedP);
colorbar;
title('t-values exceeding p-values < 0.05')
xlabel('sample points')
ylabel('occipitotemporal channels')
tvpvName = sprintf('tValExceedP'); 
saveas(tvpvFig, fullfile(savePath,tvpvName),'png');

tvpvFig_maxc = figure('Name','t-values exceeding p-values < 0.05 (efflux-influx channels)');
imagesc(tValExceedP_maxc);
colorbar;
title('t-values exceeding p-values < 0.05')
xlabel('sample points')
ylabel('efflux-influx channels')
yticks([1 2 3 4])
yticklabels({'eff L','inf L','eff R','inf R'})
tvpvName_maxc = sprintf('tValExceedP_maxc'); 
saveas(tvpvFig_maxc, fullfile(savePath,tvpvName_maxc),'png');

%% topoplot

   
% for tp = 1:length(tPoint)
%     
%     if tp == length(tPoint)
%         continue
%     else
%         tidx = timeNew>tPoint(tp)*1e-3 & timeNew<tPoint(tp+1)*1e-3;
%     end
%     
%     tValavg = double(mean(tValues(:,tidx,:),3));
%     megframe = double(mean(tValavg,2));    
%     topoAll = figure('Name','Topoplot');
%     hold on
%     toposhow(cfg, megframe);
%     colorbar
%     caxis([-2 2])
%     hold off
% 
%     titleFig = 'Distribution of overall N2pc at %d ms';
%     topoFilename = 'avgtopo_%d';
%     A1 = tPoint(tp)+5;
%     
%     title(sprintf(titleFig,A1))
%     topoAllName = sprintf(topoFilename,A1); 
%     saveas(topoAll, fullfile(savePath,topoAllName),'png');
% 
% end    

%%

% tidx = timeNew>0.200 & timeNew<0.220;
% tValavg = double(mean(tValues(maxchan,tidx,:),3));
% megframe = double(mean(tValavg,2));    
% topoAll = figure('Name','Topoplot');
% hold on
% toposhow(cfg, megframe);
% colorbar
% caxis([-2 2])
% hold off





