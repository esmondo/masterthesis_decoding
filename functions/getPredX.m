% function [X_tr,X_ts,xcl] = getPredictor(meg,data,sensors,whichtool,whichresponse,idx)
function [X_tr,xcl] = getPredX(meg,data,sensors,whichtool,whichresponse,idx)
 

%% Predictor

if whichresponse == 'corID'
    X = data(sensors,:,idx.corID); %% processing epoched data (involving correct response only)
%     X = X(:,:,idx.idxCorrBal); %% processing epoched data (from the correct response, involve balanced trials only) 

    %% heog to be rejected (if necessary)
    % heog = bandpass(1,30,meg.srate,Xeog(1,:,:));

    %% Bandpass filter X
    
%     X = baseline(X,meg.time); % baseline
    
    X_tr = bandpass(1,200,meg.srate,X);
%     X_ts = bandpass(1,200,meg.srate,X);
%     X_tr = bandpass(1,12.5,meg.srate,X);

    %% Reject Eye Movement Component
    % X_tr = rejectComponent(X_tr,squeeze(heog)); % HEOG

    %% Artifact rejection X
    if whichtool == 'meg'
        xcl = findArti_meg(X_tr);
%         X_tr(:,:,xcl) = [];
%         Y(xcl) = [];
    elseif whichtool == 'eeg'
        xcl = findArti_eeg(X_tr);
%         X_tr(:,:,xcl) = [];
%         Y(xcl) = [];
    end

    %% Convert Unit
    if whichtool == 'meg'
        X_tr = X_tr.*1e12;
%         X_ts = X_ts.*1e12;
    elseif whichtool == 'eeg'
        X_tr = X_tr.*1e6;
%         X_ts = X_ts.*1e6;
    end

    %% lowpass filter 12.5 Hz
    X_tr = lowpass(12.5,meg.srate,X_tr);
%     X_ts = lowpass(12.5,meg.srate,X_ts);

    %% DownSample 50 Hz
    FsNew = 50; 
    [ds_Xrej,t_Xrej] = dwsample(size(X_tr,1),FsNew,meg.srate,X_tr,meg.time); 
%     [ds_Xori,t_Xori] = dwsample(size(X_ts,1),FsNew,meg.srate,X_ts,meg.time);

    %% Baseline Correction
    X_tr = baseline(ds_Xrej,t_Xrej);
%     X_ts = baseline(ds_Xori,t_Xori);

    %% time of interest: 0 - 800ms
    timeCutIdx = find(t_Xrej>=0 & t_Xrej<=0.8);
    X_tr = X_tr(:,timeCutIdx,:); 
%     X_ts = X_ts(:,timeCutIdx,:);

    
elseif whichresponse == 'badID'
    
    X = data(sensors,:,idx.badID);
    
%     X = baseline(X,meg.time); % baseline
    
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
    X_tr = X(:,timeCutIdx,:);
%     X_tr = NaN;
    
    xcl = 'there was no artifact rejection';
end




%% struct number

% num.numAll = numAll;
% num.numCorID = numCorID;
% num.idxCorrBal = idxCorrBal;
% num.numBadID = numBadID;
% num.numLall = numL;
% num.numRall = numR;
% num.numLvor = numLvor;
% num.numRvor = numRvor;
% num.numLvor_x = numLvor_x;
% num.numRvor_x = numRvor_x;
% num.numLnach = numLnach;
% num.numRnach = numRnach;

    %% struct for decoding
%     dcodest.megAOI = cAll;
%     dcodest.navg = navg;
%     dcodest.nreps = nreps;
%     dcodest.srate = meg.srate;
%     dcodest.time = meg.time;
%     dcodest.Y_ori = Y_ori;
%     dcodest.Y_corID = Y_corID;
%     dcodest.Y_corbalanced = Y;
%     dcodest.eog = eog.data;
%     dcodest.k_fold = N;
%     dcodest.corID = corID;
%     dcodest.badID = badID;


%     %% EOG that needs to be rejected
%     
%     Z = eog.data(:,:,corID);
%     Z = Z(:,:,idxCorrBal);
%   
%     Xmeg = meg.data(cAll,:,corID); %% processing epoched data (involving correct response only)
%     Xmeg = Xmeg(:,:,idxCorrBal); %% processing epoched data (from the correct response, involve balanced trials only)    
%     
% %     Xmeg_x = meg.data(cAll,:,badID);
% %     Xmeg_x = Xmeg_x(:,:,idxBadBal);
%     
%     if whichdata == 4
%         Xeeg = eeg.data(pterioreegidx,:,corID);
%         Xeeg = Xeeg(:,:,idxCorrBal);
%         
% %         Xeeg_x = eeg.data(pterioreegidx,:,badID);
% %         Xeeg_x = Xeeg_x(:,:,idxBadBal);
%     end
% % 
% %     % % for it = 1:1000 % Activate this to perform a permutation test
% %     % %       Y = Y(randperm(length(Y)));
% %     % % end




%     [pred_meg,predSTF_meg,~,~] = decode_meg(Xmeg,Y,Z,N,navg,nreps,meg);
% %     [pred_meg_x,predSTF_meg_x,~,~] = decode_meg_x(Xmeg,Y,Xmeg_x,Y_x,Z,N,navg,nreps,meg);
% 
%     if whichdata == 4
%         [pred_eeg,predSTF_eeg,~,~] = decode_eeg(Xeeg,Y,Z,N,navg,nreps,meg);
% %         [pred_eeg_x,predSTF_eeg_x,~,~] = decode_eeg_x(Xeeg,Y,Xeeg_x,Y_x,Z,N,navg,nreps,meg);
%     end


%     %% Predictions
%     predictions.meg = pred_meg;
%     predictions.megSTF = predSTF_meg;    
% 
% %     predictions.meg_x = pred_meg_x;
% %     predictions.megSTF_x = predSTF_meg_x; 
% 
%     if whichdata == 4
%         predictions.eeg = pred_eeg;
%         predictions.eegSTF = predSTF_eeg;
% 
% %         predictions.eeg_x = pred_eeg_x;
% %         predictions.eegSTF_x = predSTF_eeg_x;
%     end

end