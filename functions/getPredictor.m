function [X_tr,xcl] = getPredictor(meg,data,sensors,whichtool,whichresponse,idx)
 

%% Predictor

if whichresponse == 'corID'
    X = data(sensors,:,idx.corID); %% processing epoched data (involving correct response only)
 
    %% Bandpass filter X
    X_tr = bandpass(1,200,meg.srate,X);
    
    %% Baseline Correction
    X_tr = baseline(X_tr,meg.time);

    %% Reject Eye Movement Component (if necessary)
    % heog = bandpass(1,30,meg.srate,Xeog(1,:,:));
    % X_tr = rejectComponent(X_tr,squeeze(heog)); % HEOG

    %% Artifact rejection X
    if whichtool == 'meg'
        xcl = findArti_meg(X_tr);
    elseif whichtool == 'eeg'
        xcl = findArti_eeg(X_tr);
    end

    %% Convert Unit
    if whichtool == 'meg'
        X_tr = X_tr.*1e12;
    elseif whichtool == 'eeg'
        X_tr = X_tr.*1e6;
    end

    %% lowpass filter 12.5 Hz
    X_tr = lowpass(12.5,meg.srate,X_tr);

    %% DownSample 50 Hz
    FsNew = 50; 
    [ds_Xrej,t_Xrej] = dwsample(size(X_tr,1),FsNew,meg.srate,X_tr,meg.time); 

    %% time of interest: 0 - 800ms
    timeCutIdx = find(t_Xrej>=0 & t_Xrej<=0.8);
%     X_tr = X_tr(:,timeCutIdx,:); 
    X_tr = ds_Xrej(:,timeCutIdx,:); 
    
elseif whichresponse == 'badID'
    
    X = data(sensors,:,idx.badID);
       
    X = bandpass(1,200,meg.srate,X);
    
    X = baseline(X,meg.time);
    
    if whichtool == 'meg'
        X = X.*1e12;
    elseif whichtool == 'eeg'
        X = X.*1e6;
    end
    
    X = lowpass(12.5,meg.srate,X);
    
    FsNew = 50; 
    [ds_X,t_X] = dwsample(size(X,1),FsNew,meg.srate,X,meg.time);
    
%     X = baseline(ds_X,t_X);
    
    timeCutIdx = t_X>=0 & t_X<=0.8;
%     X_tr = X(:,timeCutIdx,:);
    X_tr = ds_X(:,timeCutIdx,:);
    
    
    xcl = 'there was no artifact rejection';
end


end