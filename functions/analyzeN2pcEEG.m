function erpdata_EEG = analyzeN2pcEEG(meg,data_eeg,chanEEGname,corrTargPos,timeNew)

disp('Begin analyzing ERP...')

%% separate target by LVF and RVF 
[~,idxTargetLVF, idxTargetRVF] = hemisphaereTeilen(corrTargPos);


% EEG left and EEG right
targetLVF_eeg = data_eeg(:,:,idxTargetLVF);
targetRVF_eeg = data_eeg(:,:,idxTargetRVF);


%% Electrode selection: parietal/occipital electrodes on the Left and Right Hemishphere

% parietal/occipital electrodes located in the left hemisphere
CP1idx = find(strcmp(chanEEGname,'CP1'));
P7idx = find(strcmp(chanEEGname,'P7'));
P3idx = find(strcmp(chanEEGname,'P3'));
PO7idx = find(strcmp(chanEEGname,'PO7'));
PO3idx = find(strcmp(chanEEGname,'PO3'));
O9idx = find(strcmp(chanEEGname,'O9'));
LHidx = sort([CP1idx P7idx P3idx PO7idx PO3idx O9idx]);

% parietal/occipital electrodes located in the right hemisphere
CP2idx = find(strcmp(chanEEGname,'CP2'));
P4idx = find(strcmp(chanEEGname,'P4'));
P8idx = find(strcmp(chanEEGname,'P8'));
PO4idx = find(strcmp(chanEEGname,'PO4'));
PO8idx = find(strcmp(chanEEGname,'PO8'));
O10idx = find(strcmp(chanEEGname,'O10'));
RHidx = sort([CP2idx P4idx P8idx PO4idx PO8idx O10idx]);

% brain response in both visual field with the respected electrode: all
% parieto-occipital electrodes located in the left and right hemishpere
chanLH_LVF = targetLVF_eeg(LHidx,:,:); % left hemisphere
chanLH_RVF = targetRVF_eeg(LHidx,:,:); % left hemisphere
chanRH_LVF = targetLVF_eeg(RHidx,:,:); % right hemisphere
chanRH_RVF = targetRVF_eeg(RHidx,:,:); % right hemisphere

% concatenate parietal/occipital LH and RH electrodes
eegcontra = cat(3,chanLH_RVF,chanRH_LVF); % contralateral
eegipsi = cat(3,chanRH_RVF,chanLH_LVF);   % ipsilateral

%%%%% row 1 --> P3 + P4
%%%%% row 2 --> O9 + O10
%%%%% row 3 --> P7 + P8
%%%%% row 4 --> CP1 + CP2
%%%%% row 5 --> PO3 + PO4
%%%%% row 6 --> PO7 + PO8


%% t-test: contralateral vs ipsilateral

[~,pValues,~,stats] = ttest2(eegcontra,eegipsi,'Dim',3);
tValues = stats.tstat;

%% N2pc

sigInterv = {};
for i = 1:size(pValues,1)
    sigInterv{i} = find(pValues(i,:) < 0.05);
end

sigDiff = tValues.*double(pValues<0.05);
sigDiffx = mean(sigDiff);

% figure(1)
% hold on
% plot(timeNew,sigDiff)
% hold off

%% time range 200 msec - 300 msec post stimulus
% 
timeN2pcIdx = find(timeNew>=0.18 & timeNew<=0.35);
timeN2pcVal = timeNew(timeN2pcIdx);

[maxch1,idx_maxch1] = min(sigDiff(1,timeN2pcIdx));
[maxch2,idx_maxch2] = min(sigDiff(2,timeN2pcIdx));
[maxch3,idx_maxch3] = min(sigDiff(3,timeN2pcIdx));
[maxch4,idx_maxch4] = min(sigDiff(4,timeN2pcIdx));
[maxch5,idx_maxch5] = min(sigDiff(5,timeN2pcIdx));
[maxch6,idx_maxch6] = min(sigDiff(6,timeN2pcIdx));

tpos_maxch1 = timeN2pcVal(idx_maxch1);
tpos_maxch2 = timeN2pcVal(idx_maxch2);
tpos_maxch3 = timeN2pcVal(idx_maxch3);
tpos_maxch4 = timeN2pcVal(idx_maxch4);
tpos_maxch5 = timeN2pcVal(idx_maxch5);
tpos_maxch6 = timeN2pcVal(idx_maxch6);

%% Average of all electrodes for every target position across trials

avgLH_LVF = mean(chanLH_LVF,3);
avgLH_RVF = mean(chanLH_RVF,3);
avgRH_LVF = mean(chanRH_LVF,3);
avgRH_RVF = mean(chanRH_RVF,3);

avgContra = mean(eegcontra,3);
avgIpsi = mean(eegipsi,3);


%% lowpass filtering
avgLH_LVF = lowpass(30,meg.srate,avgLH_LVF);
avgLH_RVF = lowpass(30,meg.srate,avgLH_RVF);
avgRH_LVF = lowpass(30,meg.srate,avgRH_LVF);
avgRH_RVF = lowpass(30,meg.srate,avgRH_RVF);

avgContra = lowpass(30,meg.srate,avgContra);
avgIpsi = lowpass(30,meg.srate,avgIpsi);
avg_conminipsi = avgContra - avgIpsi; % contra-minus-ipsi


% figure(1)
% hold on
% plot(timeNew,avgContra(1,:))
% plot(timeNew,avgIpsi(1,:))
% plot(timeNew,avg_conminipsi(1,:))
% hold off
% 
% figure(2)
% hold on
% plot(timeNew,avgContra(2,:))
% plot(timeNew,avgIpsi(2,:))
% plot(timeNew,avg_conminipsi(2,:))
% hold off
% 
% figure(3)
% hold on
% plot(timeNew,avgContra(3,:))
% plot(timeNew,avgIpsi(3,:))
% plot(timeNew,avg_conminipsi(3,:))
% hold off
% 
% figure(4)
% hold on
% plot(timeNew,avgContra(4,:))
% plot(timeNew,avgIpsi(4,:))
% plot(timeNew,avg_conminipsi(4,:))
% hold off
% 
% figure(5)
% hold on
% plot(timeNew,avgContra(5,:))
% plot(timeNew,avgIpsi(5,:))
% plot(timeNew,avg_conminipsi(5,:))
% hold off
% 
% figure(6)
% hold on
% plot(timeNew,avgContra(6,:))
% plot(timeNew,avgIpsi(6,:))
% plot(timeNew,avg_conminipsi(6,:))
% hold off


%% Data Summary

erpdata_EEG.avgLH_LVF = avgLH_LVF;
erpdata_EEG.avgLH_RVF = avgLH_RVF;
erpdata_EEG.avgRH_LVF = avgRH_LVF;
erpdata_EEG.avgRH_RVF = avgRH_RVF;

erpdata_EEG.avgcontra = avgContra;
erpdata_EEG.avgipsi = avgIpsi;
erpdata_EEG.avg_conminipsi = avg_conminipsi;

erpdata_EEG.tValues = tValues;
erpdata_EEG.pValues = pValues;
erpdata_EEG.sigInterv = sigInterv;
erpdata_EEG.sigDiff = sigDiff;
erpdata_EEG.sigDiff_mean = sigDiffx;

erpdata_EEG.maxval.maxch1 = maxch1;
erpdata_EEG.maxval.maxch2 = maxch2;
erpdata_EEG.maxval.maxch3 = maxch3;
erpdata_EEG.maxval.maxch4 = maxch4;
erpdata_EEG.maxval.maxch5 = maxch5;
erpdata_EEG.maxval.maxch6 = maxch6;

erpdata_EEG.maxval.idx_maxch1 = idx_maxch1;
erpdata_EEG.maxval.idx_maxch2 = idx_maxch2;
erpdata_EEG.maxval.idx_maxch3 = idx_maxch3;
erpdata_EEG.maxval.idx_maxch4 = idx_maxch4;
erpdata_EEG.maxval.idx_maxch5 = idx_maxch5;
erpdata_EEG.maxval.idx_maxch6 = idx_maxch6;

erpdata_EEG.maxval.tpos_maxch1 = tpos_maxch1;
erpdata_EEG.maxval.tpos_maxch2 = tpos_maxch2;
erpdata_EEG.maxval.tpos_maxch3 = tpos_maxch3;
erpdata_EEG.maxval.tpos_maxch4 = tpos_maxch4;
erpdata_EEG.maxval.tpos_maxch5 = tpos_maxch5;
erpdata_EEG.maxval.tpos_maxch6 = tpos_maxch6;


end