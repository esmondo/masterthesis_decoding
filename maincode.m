

close all
clear all
clc

% for k = 1:10 % looping all navg
    
%%
% cd /export/data/wienke/data/Esmondo
% addpath('/export/data/wienke/data/Esmondo')
% addpath('/export/data/esmondo')
% addpath('/export/data/esmondo/functions')
addpath('/export/data/reichert/toolbox/MD_tools')
addpath('/export/data/reichert/toolbox/Elekta')
addpath('/export/data/reichert/toolbox/fieldtrip')
% addpath('/export/data/database/MEG/mindwandering_n2pc_sss')
addpath('/export/data/duerschm/allscripts')

addpath(genpath('/export/data/esmondo'))
savepath = '/export/data/esmondo/any_results/13subs';
% savepath = '/export/data/esmondo/any_results/oldMEGdata';

%%
if ~exist('topoplot','file')
    addpath /export/data/reichert/toolbox/ft4topoplot/
end

cfgfileMEG = '/export/data/reichert/capLayout/Neuromag_helmet.mat';
cfgfileEEG = '/export/data/reichert/capLayout/EEG_29_cartesianSorted.mat';

cfgMEG = topoprepare(cfgfileMEG);
cfgEEG = topoprepare(cfgfileEEG);
% ft_layoutplot(cfg);

%% Input 
% 
whichdata = input('Select "2" for the old dataset (only MEG) | Select "4" for the new dataset (EEG & MEG) = ');
% whichnavg = input('navg = ');

% whichdata = 4;
% whichnavg = k;

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

    cd '/export/data/esmondo/any_results/exportnew';
    myFolder = '/export/data/esmondo/any_results/exportnew';


    if ~isdir(myFolder)
        errorMessage = sprintf('Error: The following folder does not exists:\n%s', myFolder);
        uiwait(warndlg(errorMessage));
        return;
    end

    filePattern = fullfile(myFolder,'MW*newSSS.mat');
    filelist = dir(filePattern);
    endIter = length(filelist);
    
end    

%% opening display
% % resname = ['navg_' num2str(whichnavg)];
% % fprintf('----- Processing %s ----- \n', resname);

%% Execute all subjects

cd '/export/data/esmondo/any_results/13subs';
% cd '/export/data/esmondo/any_results/oldMEGdata';

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

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                               % DATA SET-UP %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define data
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

    
    % Sensors Labelling    
    chanALL = meg.header.label; % All Sensors
    
    chan_meg = meg.header.label(meg.channels); % MEG
    % chanidx = strmatch('MEG',meg.header.label); % the same with meg.channels
    % chanidx = chanidx( 1:3:end );
    
    [cL,cR,cAll] = getAOI_occitemp; % AOIs Index: Occipitotemporal region in MEG

    if whichdata == 4
        chan_eeg = meg.header.label(eeg.channels(1:29)); % EEG
        chanEEGname = {'Fz','Cz','Pz','Oz','Iz','Fp1','Fp2','F3','F4','F7'...
                   'F8','T7','T8','C3','C4','P3','P4','O9','O10',...
                   'P7','P8','FC1','FC2','CP1','CP2','PO3','PO4','PO7','PO8'}'; 
 
        % Index of PO7 and PO8 in EEG
        PO7idx = find(strcmp(chanEEGname,'PO7')); % located in the left hemisphere
        PO8idx = find(strcmp(chanEEGname,'PO8')); % located in the right hemisphere
 
        % Index of posterior channels in EEG
        pterioreeg = {'P3','P7','PO3','PO7',...
               'P4','P8','PO4','PO8',...
               'Pz','Oz'};
        pterioreegidx = find(ismember(chanEEGname, pterioreeg));
    end
   
    chan_eog = meg.header.label(eog.channels); % EOG
 

    % indexing essential parameters
    [x1_meg, x2_meg, x3_meg(i)] = size(data_meg);

    if whichdata == 4
        [x1_eeg, x2_eeg, x3_eeg(i)] = size(data_eeg);
    end


    %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % ERP PREPROCESSING %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [data_meg,data_eeg,data_eog,corrTargPos,num(i),timeIdx,timeNew] = preprocess(meg,eeg,eog,whichdata);
    
    if whichdata == 2
        eeg = [];
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             % ERP ANALYSIS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % MEG %
    erpdata_MEG(i) = analyzeN2pcMEG(meg,data_meg,cL,cR,corrTargPos,timeNew);
  
    % EEG %
    if whichdata == 4 
        erpdata_EEG(i) = analyzeN2pcEEG(meg,data_eeg,chanEEGname,corrTargPos,timeNew);
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % DECODING %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% involve correct trials only
    [~,corID] = corrResponse(meg.side,meg.response);
    numCorID(i) = numel(corID);
    numAll(i) = numel(meg.side);
    badID = setdiff([1:numAll(i)],corID); % use this for TESTING
    numBadID(i) = numel(badID);

    %% LABEL
    [Y_ori,~,~] = hemisphaereTeilen(meg.position);
    [Y_corID,idxLVF,idxRVF] = hemisphaereTeilen(meg.position(corID));
    numLvor(i) = numel(idxLVF);
    numRvor(i) = numel(idxRVF);

    [Y,idxLVF,idxRVF,idxCorr] = balancetrial(Y_corID); % sample balancing
    numLnach(i) = numel(idxLVF);
    numRnach(i) = numel(idxRVF);

    %% EOG that needs to be rejected
    Z = eog.data(:,:,corID);
    Z = Z(:,:,idxCorr);

    %% PREDICTOR
    Xmeg = meg.data(cAll,:,corID); %% processing epoched data (involving correct response only)
    Xmeg = Xmeg(:,:,idxCorr); %% processing epoched data (from the correct response, involve balanced trials only)
    if whichdata == 4
        Xeeg = eeg.data(pterioreegidx,:,corID);
        Xeeg = Xeeg(:,:,idxCorr);
    end

    % % for it = 1:1000 % Activate this to perform a permutation test
    % %       Y = Y(randperm(length(Y)));
    % % end
 
    N = 10; % N-fold / k-fold
    navg = whichnavg; % navg: average trial
    nreps = 10;

    %% struct for decoding
    dcodest.megAOI = cAll;
    dcodest.navg = navg;
    dcodest.nreps = nreps;
    dcodest.srate = meg.srate;
    dcodest.time = meg.time;
    dcodest.Y_ori = Y_ori;
    dcodest.Y_corID = Y_corID;
    dcodest.Y_corbalanced = Y;
    dcodest.eog = eog.data;
    dcodest.k_fold = N;
    dcodest.corID = corID;
    dcodest.badID = badID;

    %% Classifications
    [pred_meg,predSTF_meg,alg1_meg,alg2_meg] = decode_meg(Xmeg,Y,Z,N,navg,nreps,meg);
    % [pred_meg_LDA,predSTF_meg_LDA,alg1_meg_LDA,alg2_meg_LDA] = decode_meg_LDA(Xmeg,Y,Z,N,navg,nreps,meg);
    [pred_meg2,predSTF_meg2,alg1_meg2,alg2_meg2] = decode_meg_2(Xmeg,Y,Z,dcodest,meg); 

    if whichdata == 4
        [pred_eeg,predSTF_eeg,alg1_eeg,alg2_eeg] = decode_eeg(Xeeg,Y,Z,N,navg,nreps,meg);
%       [pred_eeg_LDA,predSTF_eeg_LDA,alg1_eeg_LDA,alg2_eeg_LDA] = decode_eeg_LDA(Xeeg,Y,Z,N,navg,nreps,meg);
    end

    %% Predictions
    predictions.meg = pred_meg;
    predictions.megSTF = predSTF_meg;     
%     predictions.meg2 = pred_meg2;
%     predictions.megSTF2 = predSTF_meg2;
    
%     predictions.meg_LDA = pred_meg_LDA;
%     predictions.megSTF_LDA = predSTF_meg_LDA;

    if whichdata == 4
        predictions.eeg = pred_eeg;
        predictions.eegSTF = predSTF_eeg;
    
    %     predictions.eeg_LDA = pred_eeg_LDA;
    %     predictions.eegSTF_LDA = predSTF_eeg_LDA;
    end


    %% Results
    results(i).MEG = calculate_performance(predictions.meg,Y);
    results(i).MEGstf = calculate_performance(predictions.megSTF,Y);    
%     results(i).MEG2 = calculate_performance(predictions.meg2,Y);
%     results(i).MEGstf2 = calculate_performance(predictions.megSTF2,Y);
    
%     results(i).MEG_LDA = calculate_performance(predictions.meg_LDA,Y);
%     results(i).MEGstf_LDA = calculate_performance(predictions.megSTF_LDA,Y);
    
    if whichdata == 4
        results(i).EEG = calculate_performance(predictions.eeg,Y);
        results(i).EEGstf = calculate_performance(predictions.eegSTF,Y);
    
    %     results(i).EEG_LDA = calculate_performance(predictions.eeg_LDA,Y);
    %     results(i).EEGstf_LDA = calculate_performance(predictions.eegSTF_LDA,Y);
    end

    %% Confidence Interval (only for permutation test or single-trial (navg=1) )

    % interv(i) = 1.96.*sqrt((missclass(i,:).*(1-missclass(i,:)))./length(Y));
    % interv_stf(i) = 1.96.*sqrt((missclass_stf(i,:).*(1-missclass_stf(i,:)))./length(Y));


    %% Accuracy per subject to be displayed

    if whichdata == 2
        fprintf('MEG Accuracy = %3.2f%% , ', results(i).MEG.perf);
        fprintf('MEG Accuracy with STF =  %3.2f%% \n', results(i).MEGstf.perf);
    
    %     fprintf('MEG Accuracy (LDA) = %3.2f%% , ', results(i).MEG_LDA.perf);
    %     fprintf('MEG Accuracy with STF (LDA) =  %3.2f%% \n', results(i).MEGstf_LDA.perf);
    
    elseif whichdata == 4
        fprintf('MEG Accuracy = %3.2f%% , ', results(i).MEG.perf);
        fprintf('MEG Accuracy with STF =  %3.2f%% | ', results(i).MEGstf.perf);
        fprintf('EEG Accuracy = %3.2f%% , ', results(i).EEG.perf);
        fprintf('EEG Accuracy with STF =  %3.2f%% \n', results(i).EEGstf.perf);
    
    %     fprintf('MEG Accuracy (LDA) = %3.2f%% , ', results(i).MEG_LDA.perf);
    %     fprintf('MEG Accuracy with STF (LDA) =  %3.2f%% | ', results(i).MEGstf_LDA.perf);
    %     fprintf('EEG Accuracy (LDA) = %3.2f%% , ', results(i).EEG_LDA.perf);
    %     fprintf('EEG Accuracy with STF (LDA) =  %3.2f%% \n', results(i).EEGstf_LDA.perf);
    end



end

 

%% Summary of performance

for i=1:endIter
    
    % accuracy
    allPerf_MEG(i) = cat(2,results(i).MEG.perf);
    allPerf_MEGstf(i) = cat(2,results(i).MEGstf.perf);   
%     allPerf_MEG_LDA(i) = cat(2,results(i).MEG_LDA.perf);
%     allPerf_MEGstf_LDA(i) = cat(2,results(i).MEGstf_LDA.perf);
    
    if whichdata == 4
        allPerf_EEG(i) = cat(2,results(i).EEG.perf);
        allPerf_EEGstf(i) = cat(2,results(i).EEGstf.perf);        
%         allPerf_EEG_LDA(i) = cat(2,results(i).EEG_LDA.perf);
%         allPerf_EEGstf_LDA(i) = cat(2,results(i).EEGstf_LDA.perf);
    end
    
        
    % recall for class 0 
    recall_c1_MEG(i,:) = cat(1,results(i).MEG.recall_c1);
    recall_c1_MEGstf(i,:) = cat(1,results(i).MEGstf.recall_c1);    
%     recall_c1_MEG_LDA(i,:) = cat(1,results(i).MEG_LDA.recall_c1);
%     recall_c1_MEGstf_LDA(i,:) = cat(1,results(i).MEGstf_LDA.recall_c1);
    
    if whichdata == 4
        recall_c1_EEG(i,:) = cat(1,results(i).EEG.recall_c1);
        recall_c1_EEGstf(i,:) = cat(1,results(i).EEGstf.recall_c1);        
%         recall_c1_EEG_LDA(i,:) = cat(1,results(i).EEG_LDA.recall_c1);
%         recall_c1_EEGstf_LDA(i,:) = cat(1,results(i).EEGstf_LDA.recall_c1);
    end
    
    
    % recall class 1 
    recall_c2_MEG(i,:) = cat(1,results(i).MEG.recall_c2);
    recall_c2_MEGstf(i,:) = cat(1,results(i).MEGstf.recall_c2);    
%     recall_c2_MEG_LDA(i,:) = cat(1,results(i).MEG_LDA.recall_c2);
%     recall_c2_MEGstf_LDA(i,:) = cat(1,results(i).MEGstf_LDA.recall_c2);
    
    if whichdata == 4
        recall_c2_EEG(i,:) = cat(1,results(i).EEG.recall_c2);
        recall_c2_EEGstf(i,:) = cat(1,results(i).EEGstf.recall_c2);        
%         recall_c2_EEG_LDA(i,:) = cat(1,results(i).EEG_LDA.recall_c2);
%         recall_c2_EEGstf_LDA(i,:) = cat(1,results(i).EEGstf_LDA.recall_c2);
    end   
    
end

%% Table sum
sum_navg = zeros(numel(results),numel(fieldnames(results)));

if whichdata == 2
    sum_navg(:,1) = allPerf_MEG';
    sum_navg(:,2) = allPerf_MEGstf';
elseif whichdata == 4
    sum_navg(:,1) = allPerf_EEG';
    sum_navg(:,2) = allPerf_EEGstf';
    sum_navg(:,3) = allPerf_MEG';
    sum_navg(:,4) = allPerf_MEGstf';
    sum_navg(:,5) = allPerf_EEG_LDA';
    sum_navg(:,6) = allPerf_EEGstf_LDA';
    sum_navg(:,7) = allPerf_MEG_LDA';
    sum_navg(:,8) = allPerf_MEGstf_LDA';
end


%% Results struct
results_eachsubject.perf.MEG = allPerf_MEG;
results_eachsubject.perf.MEGstf = allPerf_MEGstf;
results_eachsubject.recall_c1.MEG = recall_c1_MEG;
results_eachsubject.recall_c1.MEGstf = recall_c1_MEGstf;
results_eachsubject.recall_c2.MEG = recall_c2_MEG;
results_eachsubject.recall_c2.MEGstf = recall_c2_MEGstf;
% results_eachsubject.perf.MEG_LDA = allPerf_MEG_LDA;
% results_eachsubject.perf.MEGstf_LDA = allPerf_MEGstf_LDA;
% results_eachsubject.recall_c1.MEG_LDA = recall_c1_MEG_LDA;
% results_eachsubject.recall_c1.MEGstf_LDA = recall_c1_MEGstf_LDA;
% results_eachsubject.recall_c2.MEG_LDA = recall_c2_MEG_LDA;
% results_eachsubject.recall_c2.MEGstf_LDA = recall_c2_MEGstf_LDA;

if whichdata == 4
    results_eachsubject.perf.EEG = allPerf_EEG;
    results_eachsubject.perf.EEGstf = allPerf_EEGstf;
    results_eachsubject.recall_c1.EEG = recall_c1_EEG;
    results_eachsubject.recall_c1.EEGstf =recall_c1_EEGstf;
    results_eachsubject.recall_c2.EEG = recall_c2_EEG;
    results_eachsubject.recall_c2.EEGstf =recall_c2_EEGstf;   
%     results_eachsubject.perf.EEG_LDA = allPerf_EEG_LDA;
%     results_eachsubject.perf.EEGstf_LDA = allPerf_EEGstf_LDA;
%     results_eachsubject.recall_c1.EEG_LDA = recall_c1_EEG_LDA;
%     results_eachsubject.recall_c1.EEGstf_LDA =recall_c1_EEGstf_LDA;
%     results_eachsubject.recall_c2.EEG_LDA = recall_c2_EEG_LDA;
%     results_eachsubject.recall_c2.EEGstf_LDA =recall_c2_EEGstf_LDA;
end

 %% Accuracy comparison
perfComparison_MEG = cat(1,results_eachsubject.perf.MEG,results_eachsubject.perf.MEGstf);
% perfComparison_MEG_LDA = cat(1,results_eachsubject.perf.MEG_LDA,results_eachsubject.perf.MEGstf_LDA);

if whichdata == 4
    perfComparison_EEG = cat(1,results_eachsubject.perf.EEG, results_eachsubject.perf.EEGstf);
%     perfComparison_EEG_LDA = cat(1,results_eachsubject.perf.EEG_LDA, results_eachsubject.perf.EEGstf_LDA);
end


%% Grand Results
results_navg.perf.MEG = mean(allPerf_MEG);
results_navg.perf.MEGstf = mean(allPerf_MEGstf);
results_navg.perf_std.MEG = std(allPerf_MEG);
results_navg.perf_std.MEGstf = std(allPerf_MEGstf);
results_navg.recall_c1.MEG = mean(recall_c1_MEG);
results_navg.recall_c1.MEGstf = mean(recall_c1_MEGstf);
results_navg.recall_c2.MEG = mean(recall_c2_MEG);
results_navg.recall_c2.MEGstf = mean(recall_c2_MEGstf);
% results_navg.perf.MEG_LDA = mean(allPerf_MEG_LDA);
% results_navg.perf.MEGstf_LDA = mean(allPerf_MEGstf_LDA);
% results_navg.perf_std.MEG_LDA = std(allPerf_MEG_LDA);
% results_navg.perf_std.MEGstf_LDA = std(allPerf_MEGstf_LDA);
% results_navg.recall_c1.MEG_LDA = mean(recall_c1_MEG_LDA);
% results_navg.recall_c1.MEGstf_LDA = mean(recall_c1_MEGstf_LDA);
% results_navg.recall_c2.MEG_LDA = mean(recall_c2_MEG_LDA);
% results_navg.recall_c2.MEGstf_LDA = mean(recall_c2_MEGstf_LDA);

if whichdata == 4
    results_navg.perf.EEG = mean(allPerf_EEG);
    results_navg.perf.EEGstf = mean(allPerf_EEGstf);
    results_navg.perf_std.EEG = std(allPerf_EEG);
    results_navg.perf_std.EEGstf = std(allPerf_EEGstf);
    results_navg.recall_c1.EEG = mean(recall_c1_EEG);
    results_navg.recall_c1.EEGstf = mean(recall_c1_EEGstf);
    results_navg.recall_c2.EEG = mean(recall_c2_EEG);
    results_navg.recall_c2.EEGstf = mean(recall_c2_EEGstf);   
%     results_navg.perf.EEG_LDA = mean(allPerf_EEG_LDA);
%     results_navg.perf.EEGstf_LDA = mean(allPerf_EEGstf_LDA);
%     results_navg.perf_std.EEG_LDA = std(allPerf_EEG_LDA);
%     results_navg.perf_std.EEGstf_LDA = std(allPerf_EEGstf_LDA);
%     results_navg.recall_c1.EEG_LDA = mean(recall_c1_EEG_LDA);
%     results_navg.recall_c1.EEGstf_LDA = mean(recall_c1_EEGstf_LDA);
%     results_navg.recall_c2.EEG_LDA = mean(recall_c2_EEG_LDA);
%     results_navg.recall_c2.EEGstf_LDA = mean(recall_c2_EEGstf_LDA);
end



%% Saving Results

resname = ['navg_' num2str(navg)];
fprintf('----- Saving the results of %s ----- \n', resname);
save(fullfile(savepath,resname),'sum_navg','results_navg','results_eachsubject');

close all
clear all

% end % for looping all navg (k = 1:10)

 %% used channel
% usedChan = cat(1,effluxL,influxL,effluxR,influxR);
% usedChanName = chanMEG(usedChan);

%% ERP PLOT

% time range 200 msec - 300 msec post stimulus
timeN2pcIdx = find(timeNew>=0.2 & timeNew<=0.4);
timeN2pcVal = timeNew(timeN2pcIdx);

% %t-values exceeding p-threshold (significantly different)
% sigDiff_LH = tVal_LH.*double(pVal_LH<0.05);
% avgSigDiff_LH = mean(sigDiff_LH);
% [tVn2pc_LH_val,tVn2pc_LH_pos] = max(avgSigDiff_LH(timeN2pcIdx));
% 
% %t-values exceeding p-threshold (significantly different)
% sigDiff_RH = tVal_RH.*double(pVal_RH<0.05);
% avgSigDiff_RH = mean(sigDiff_RH);
% [tVn2pc_RH_val,tVn2pc_RH_pos] = min(avgSigDiff_RH(timeN2pcIdx));

%% Grand Average N2pc on each hemisphere: MEG

GA_MEG.avgLeftHem_LVF = zeros(endIter,length(timeNew));
GA_MEG.avgLeftHem_RVF = zeros(endIter,length(timeNew));
GA_MEG.avgRightHem_LVF = zeros(endIter,length(timeNew));
GA_MEG.avgRightHem_RVF = zeros(endIter,length(timeNew));

for k = 1:endIter
    GA_MEG.avgLeftHem_LVF(k,:) = erpdata_MEG(k).avgLeftHem_LVF;
    GA_MEG.avgLeftHem_RVF(k,:) = erpdata_MEG(k).avgLeftHem_RVF;
    GA_MEG.avgRightHem_LVF(k,:) = erpdata_MEG(k).avgRightHem_LVF;
    GA_MEG.avgRightHem_RVF(k,:) = erpdata_MEG(k).avgRightHem_RVF;    
end


%left hem
GA_MEG.GA_avgLeftHem_LVF = mean(GA_MEG.avgLeftHem_LVF);
GA_MEG.GA_avgLeftHem_RVF = mean(GA_MEG.avgLeftHem_RVF);
GA_MEG.GA_N2pcLH = GA_MEG.GA_avgLeftHem_LVF - GA_MEG.GA_avgLeftHem_RVF;

%right hem
GA_MEG.GA_avgRightHem_LVF = mean(GA_MEG.avgRightHem_LVF);
GA_MEG.GA_avgRightHem_RVF = mean(GA_MEG.avgRightHem_RVF);
GA_MEG.GA_N2pcRH = GA_MEG.GA_avgRightHem_LVF - GA_MEG.GA_avgRightHem_RVF;


%% GA plot MEG
GAleftFig = figure('Name','GA N2pc in Left Hemisphere');
hold on
plot(timeNew,GA_MEG.GA_avgLeftHem_LVF,'b');
plot(timeNew,GA_MEG.GA_avgLeftHem_RVF,'r--');
plot(timeNew,GA_MEG.GA_N2pcLH,'k','LineWidth',2);
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA N2pc in Left Hemisphere')
legend('LVF Target','RVF Target','LVF-RVF')
yL = get(gca,'ylim');
ylim([-1.5e-13 1.5e-13])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');
% line([0.18 0.18],ylim,'Color','k','LineStyle','--');
% line([0.3 0.3],ylim,'Color','k','LineStyle','--');
% GAleftFigName = sprintf('GAleftHem');
% saveas(GAleftFig, fullfile(savePath,GAleftFigName),'png');

GArightFig = figure('Name','GA N2pc Right Hemisphere');
hold on
plot(timeNew,GA_MEG.GA_avgRightHem_LVF,'b');
plot(timeNew,GA_MEG.GA_avgRightHem_RVF,'r--');
plot(timeNew,GA_MEG.GA_N2pcRH,'k','LineWidth',2);
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA N2pc Right Hemisphere')
legend('LVF Target','RVF Target','LVF-RVF')
yL = get(gca,'ylim');
ylim([-1.5e-13 1.5e-13])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');
% line([0.18 0.18],ylim,'Color','k','LineStyle','--');
% line([0.3 0.3],ylim,'Color','k','LineStyle','--');
% GArightFigName = sprintf('GArightHem'); 
% saveas(GArightFig, fullfile(savePath,GArightFigName),'png');

%% PEAK MEG

% check the GA N2pc amplitude
[GA_MEG.GA_N2pcLHpeak_val,GA_MEG.GA_N2pcLHpeak_id] = max(GA_MEG.GA_N2pcLH(timeN2pcIdx));
[GA_MEG.GA_N2pcRHpeak_val,GA_MEG.GA_N2pcRHpeak_id] = min(GA_MEG.GA_N2pcRH(timeN2pcIdx));

% time point where the GA N2pc peaked
GA_MEG.GA_N2pcLHpeak_pos = timeN2pcVal(GA_MEG.GA_N2pcLHpeak_id);
GA_MEG.GA_N2pcRHpeak_pos = timeN2pcVal(GA_MEG.GA_N2pcRHpeak_id);



%% Grand Average N2pc from contralateral and ipsilateral signals: EEG

GA_EEG.avgPO7_LVF = zeros(endIter,length(timeNew));
GA_EEG.avgPO7_RVF = zeros(endIter,length(timeNew));
GA_EEG.avgPO8_LVF = zeros(endIter,length(timeNew));
GA_EEG.avgPO8_RVF = zeros(endIter,length(timeNew));
GA_EEG.avgContra = zeros(endIter,length(timeNew));
GA_EEG.avgIpsi = zeros(endIter,length(timeNew));

for k = 1:endIter
    GA_EEG.avgPO7_LVF(k,:) = erpdata_EEG(k).avgPO7_LVF;
    GA_EEG.avgPO7_RVF(k,:) = erpdata_EEG(k).avgPO7_RVF;
    GA_EEG.avgPO8_LVF(k,:) = erpdata_EEG(k).avgPO8_LVF;
    GA_EEG.avgPO8_RVF(k,:) = erpdata_EEG(k).avgPO8_RVF;
    GA_EEG.avgContra(k,:) = erpdata_EEG(k).contra;
    GA_EEG.avgIpsi(k,:) = erpdata_EEG(k).ipsi;
end

GA_EEG.GA_avgPO7_LVF = mean(GA_EEG.avgPO7_LVF);
GA_EEG.GA_avgPO7_RVF = mean(GA_EEG.avgPO7_RVF);

GA_EEG.GA_avgPO8_LVF = mean(GA_EEG.avgPO8_LVF);
GA_EEG.GA_avgPO8_RVF = mean(GA_EEG.avgPO8_RVF);

GA_EEG.GA_avgContra = mean(GA_EEG.avgContra);
GA_EEG.GA_avgIpsi = mean(GA_EEG.avgIpsi);
GA_EEG.GA_avgConIpsi = GA_EEG.GA_avgContra - GA_EEG.GA_avgIpsi;

%% GA plot EEG
GA_PO7 = figure('Name','GA PO7 / Left Hemisphere');
hold on
plot(timeNew,GA_EEG.GA_avgPO7_LVF)
plot(timeNew,GA_EEG.GA_avgPO7_RVF)
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA PO7 / Left Hemisphere')
legend('PO7 LVF','PO7 RVF')
yL = get(gca,'ylim');
ylim([-5.5e-6 4e-6])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');
% line([0.2 0.2],ylim,'Color','k','LineStyle','--');
% line([0.3 0.3],ylim,'Color','k','LineStyle','--');

GA_PO8 = figure('Name','GA PO8 / Right Hemisphere');
hold on
plot(timeNew,GA_EEG.GA_avgPO8_LVF)
plot(timeNew,GA_EEG.GA_avgPO8_RVF)
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA PO8 / Right Hemisphere')
legend('PO8 LVF','PO8 RVF')
yL = get(gca,'ylim');
ylim([-5.5e-6 4e-6])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');
% line([0.2 0.2],ylim,'Color','k','LineStyle','--');
% line([0.3 0.3],ylim,'Color','k','LineStyle','--');

GAconipsi = figure('Name','GA Contralateral and Ipsilateral');
hold on
plot(timeNew,GA_EEG.GA_avgContra) % contralateral
plot(timeNew,GA_EEG.GA_avgIpsi) % ipsilateral
plot(timeNew,(GA_EEG.GA_avgContra - GA_EEG.GA_avgIpsi))
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA Contralateral and Ipsilateral')
legend('Contralateral Target','Ipsilateral Target','Contralateral-Ipsilateral')
yL = get(gca,'ylim');
ylim([-5.5e-6 4e-6])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');
% line([0.2 0.2],ylim,'Color','k','LineStyle','--');
% line([0.3 0.3],ylim,'Color','k','LineStyle','--');

%% PEAK EEG

% check the GA contra-ipsi amplitude
[GA_EEG.GA_contrapeak_val,GA_EEG.GA_contrapeak_id] = max(GA_EEG.GA_avgContra(timeN2pcIdx));
[GA_EEG.GA_ipsipeak_val,GA_EEG.GA_ipsipeak_id] = max(GA_EEG.GA_avgIpsi(timeN2pcIdx));
[GA_EEG.GA_conipsipeak_val,GA_EEG.GA_conipsipeak_id] = max(GA_EEG.GA_avgConIpsi(timeN2pcIdx));

% time point where the GA N2pc peaked
GA_EEG.GA_conipsipeak_pos = timeN2pcVal(GA_EEG.GA_conipsipeak_id);


%% t-values exceeding p-threshold


% tPoint = 175:5:325;


% tValExceedP = squeeze(mean(sigDiff,3));
% tValExceedP_maxc = squeeze(mean(sigDiff(maxchan,:,:),3));



% tvpvFig = figure('Name','t-values exceeding p-values < 0.05');
% imagesc(sigDiff_RH);
% colorbar;
% title('t-values exceeding p-values < 0.05')
% xlabel('sample points')
% ylabel('occipitotemporal channels')
% tvpvName = sprintf('tValExceedP'); 
% saveas(tvpvFig, fullfile(savePath,tvpvName),'png');
 
% tvpvFig_maxc = figure('Name','t-values exceeding p-values < 0.05 (efflux-influx channels)');
% imagesc(tValExceedP_maxc);
% colorbar;
% title('t-values exceeding p-values < 0.05')
% xlabel('sample points')
% ylabel('efflux-influx channels')
% yticks([1 2 3 4])
% yticklabels({'eff L','inf L','eff R','inf R'})
% tvpvName_maxc = sprintf('tValExceedP_maxc'); 
% saveas(tvpvFig_maxc, fullfile(savePath,tvpvName_maxc),'png');

%% topoplot
% alltValues = zeros(endIter,length(timeNew));
%    
% for k = 1:endIter
%     alltValues(k,:) = erpdata_EEG(k).tValues;
% end
% 
% avg_tValues = mean(alltValues);

% for tp = 1:length(tPoint)
%     
%     if tp == length(tPoint)
%         continue
%     else
%         tidx = timeNew>tPoint(tp)*1e-3 & timeNew<tPoint(tp+1)*1e-3;
%     end
%     
%     ft_layoutplot(cfgfileEEG);
%     tValavg = double(mean(alltValues(:,tidx)));
%     megframe = double(mean(tValavg,2));    
%     topoAll = figure('Name','Topoplot');
%     hold on
%     toposhow(cfg, megframe);
%     colorbar
%     caxis([-2 2])
%     hold off
% 
% %     titleFig = 'Distribution of overall N2pc at %d ms';
% %     topoFilename = 'avgtopo_%d';
% %     A1 = tPoint(tp)+5;
% %     
% %     title(sprintf(titleFig,A1))
% %     topoAllName = sprintf(topoFilename,A1); 
% %     saveas(topoAll, fullfile(savePath,topoAllName),'png');
% 
% end    

%%

% tidx = timeNew>0.200 & timeNew<0.400;
% tValavg = alltValues(tidx);
% eegframe = mean(tValavg);    
% topoAll = figure('Name','Topoplot');
% hold on
% toposhow(cfgEEG, tValavg);
% colorbar
% caxis([-2 2])
% hold off

%%
% tidx = timeNew>0.175 & timeNew<0.185; % average range around 180 ms
% eegframe = double(mean(avg_tValues(:,tidx),2)); % convert to double
% 
% figure;
% toposhow(cfgEEG, eegframe);


