function [X_tr,X_ts,Y,Y_badID,num] = getXY(meg,data,sensors,whichtool,whichresponse)

 

%% EOG that needs to be rejected
    
% Xeog = eog(:,:,corID);
% Xeog = Xeog(:,:,idxCorrBal);

%% Predictor

if whichresponse == 'corID'
    X = data(sensors,:,corID); %% processing epoched data (involving correct response only)
    X = X(:,:,idxCorrBal); %% processing epoched data (from the correct response, involve balanced trials only) 

    %% heog to be rejected (if necessary)
    % heog = bandpass(1,30,meg.srate,Xeog(1,:,:));

    %% Bandpass filter X
    X_tr = bandpass(1,200,meg.srate,X);
    X_ts = bandpass(1,200,meg.srate,X);

    %% Reject Eye Movement Component
    % X_tr = rejectComponent(X_tr,squeeze(heog)); % HEOG

    %% Artifact rejection X
    if whichtool == 'meg'
        xcl = findArti_meg(X_tr);
%         X_tr(:,:,xcl) = [];
%         Y(xcl) = [];
        num.artifactIdx = xcl;
    elseif whichtool == 'eeg'
        xcl = findArti_eeg(X_tr);
%         X_tr(:,:,xcl) = [];
%         Y(xcl) = [];
        num.artifactIdx = xcl;
    end

    %% Convert Unit
    if whichtool == 'meg'
        X_tr = X_tr.*1e12;
        X_ts = X_ts.*1e12;
    elseif whichtool == 'eeg'
        X_tr = X_tr.*1e6;
        X_ts = X_ts.*1e6;
    end

    %% lowpass filter 12.5 Hz
    X_tr = lowpass(12.5,meg.srate,X_tr);
    X_ts = lowpass(12.5,meg.srate,X_ts);

    %% DownSample 50 Hz
    FsNew = 50; 
    [ds_Xrej,t_Xrej] = dwsample(size(X_tr,1),FsNew,meg.srate,X_tr,meg.time); 
    [ds_Xori,t_Xori] = dwsample(size(X_ts,1),FsNew,meg.srate,X_ts,meg.time);

    %% Baseline Correction
    X_tr = baseline(ds_Xrej,t_Xrej);
    X_ts = baseline(ds_Xori,t_Xori);

    %% time of interest: 0 - 800ms
    timeCutIdx = find(t_Xori>=0 & t_Xori<=0.8);
    X_tr = X_tr(:,timeCutIdx,:); 
    X_ts = X_ts(:,timeCutIdx,:);

    
elseif whichresponse == 'badID'
    
    X = data(sensors,:,badID);
    
    X = bandpass(1,200,meg.srate,X);
    
    if whichtool == 'meg'
        X = X.*1e12;
    elseif whichtool == 'eeg'
        X = X.*1e6;
    end
    
    X = lowpass(12.5,meg.srate,X);
    
    FsNew = 50; 
    [ds_X,t_X] = dwsample(size(X,1),FsNew,meg.srate,X,meg.time);
    
    X = baseline(ds_X,t_X);
    
    timeCutIdx = t_X>=0 & t_X<=0.8;
    X_ts = X(:,timeCutIdx,:);
    X_tr = NaN;
    
end




%% struct number

num.numAll = numAll;
num.numCorID = numCorID;
num.idxCorrBal = idxCorrBal;
num.numBadID = numBadID;
num.numLall = numL;
num.numRall = numR;
num.numLvor = numLvor;
num.numRvor = numRvor;
num.numLvor_x = numLvor_x;
num.numRvor_x = numRvor_x;
num.numLnach = numLnach;
num.numRnach = numRnach;


end