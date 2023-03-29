function erpdata_MEG = analyzeN2pcMEG(meg,data_meg,chanOT_LH,chanOT_RH,corrTargPos,timeNew)

disp('Begin analyzing ERF...')

% time range 200 msec - 300 msec post stimulus
timeN2pcIdx = find(timeNew>=0.2 & timeNew<=0.3);
timeN2pcVal = timeNew(timeN2pcIdx);

%% data separation

% separate target by LVF and RVF 
[~,idxTargetLVF, idxTargetRVF] = hemisphaereTeilen(corrTargPos);
numLeftTrial = length(idxTargetLVF);
numRightTrial = length(idxTargetRVF);

% MEG left and MEG right
targetLVF_meg = data_meg(:,:,idxTargetLVF);
targetRVF_meg = data_meg(:,:,idxTargetRVF);

%%
targetLVF_meg = bandpass(0.5,30,meg.srate, targetLVF_meg);
targetRVF_meg = bandpass(0.5,30,meg.srate, targetRVF_meg);

avg_targetLVF_meg = mean(targetLVF_meg,3);
avg_targetRVF_meg = mean(targetRVF_meg,3);

avg_LVFminRVF_meg = avg_targetLVF_meg - avg_targetRVF_meg;

% figure(1)
% hold on
% plot(timeNew,avg_targetLVF_meg)
% hold off
% 
% figure(2)
% hold on
% plot(timeNew,avg_targetRVF_meg)
% hold off
% 
% figure(3)
% hold on
% plot(timeNew,avg_LVFminRVF_meg(23,:))
% hold off


% % LH
% % xl_LHchan = avg_targetLVF_meg(chanOT_LH,timeN2pcIdx);
% xl_LHchan = avg_LVFminRVF_meg(chanOT_LH,:);
% xl_LHchan(xl_LHchan <= 0) = NaN;
% [xl_LHchan,efflux_LHt] = max(xl_LHchan,[],2);
% [efflux_LHval,efflux_LHChan] = max(xl_LHchan); 
% efflux_LHt = efflux_LHt(efflux_LHChan);
% % efflux_LHt = timeN2pcVal(efflux_LHt);
% efflux_LHt = timeNew(efflux_LHt);
% 
% % xl_LHchan = avg_targetLVF_meg(chanOT_LH,timeN2pcIdx);
% xl_LHchan = avg_LVFminRVF_meg(chanOT_LH,:);
% xl_LHchan(xl_LHchan >= 0) = NaN; 
% [xl_LHchan,influx_LHt] = min(xl_LHchan,[],2);
% [influx_LHval,influx_LHChan] = min(xl_LHchan); 
% influx_LHt = influx_LHt(influx_LHChan);
% % influx_LHt = timeN2pcVal(influx_LHt);
% influx_LHt = timeNew(influx_LHt);
% 
% % RH
% % xr_RHchan = avg_targetRVF_meg(chanOT_RH,timeN2pcIdx); 
% xr_RHchan = avg_LVFminRVF_meg(chanOT_RH,:);
% xr_RHchan(xr_RHchan <= 0) = NaN;
% [xr_RHchan,efflux_RHt] = max(xr_RHchan,[],2);
% [efflux_RHval,efflux_RHChan] = max(xr_RHchan); 
% efflux_RHt = efflux_RHt(efflux_RHChan);
% % efflux_RHt = timeN2pcVal(efflux_RHt);
% efflux_RHt = timeNew(efflux_RHt);
% 
% % xr_RHchan = avg_targetRVF_meg(chanOT_RH,timeN2pcIdx);
% xr_RHchan = avg_LVFminRVF_meg(chanOT_RH,:);
% xr_RHchan(xr_RHchan >= 0) = NaN; 
% [xr_RHchan,influx_RHt] = min(xr_RHchan,[],2);
% [influx_RHval,influx_RHChan] = min(xr_RHchan); 
% influx_RHt = influx_RHt(influx_RHChan);
% % influx_RHt = timeN2pcVal(influx_RHt);
% influx_RHt = timeNew(influx_RHt);

%% t-test

% compare targets in LVF vs targets in RVF
[~,pValues,~,stats] = ttest2(targetLVF_meg,targetRVF_meg,'Dim',3);

% the t-values
tValues = stats.tstat;

% t-values exceeding p-threshold (significantly different)
sigDiff = tValues.*double(pValues<0.05);


% [~,pValues2,~,stats2] = ttest2(targetLVF_meg(:,timeN2pcIdx,:),targetRVF_meg(:,timeN2pcIdx,:),'Dim',3);
% tValues2 = stats2.tstat;
% sigDiff2 = tValues2.*double(pValues2<0.05);


% tVal_LH = tValues(chanOT_LH,timeN2pcIdx);
% tVal_RH = tValues(chanOT_RH,timeN2pcIdx);
% pVal_LH = pValues(chanOT_LH,timeN2pcIdx);
% pVal_RH = pValues(chanOT_RH,timeN2pcIdx);
tVal_LH = tValues(chanOT_LH,:);
tVal_RH = tValues(chanOT_RH,:);
pVal_LH = pValues(chanOT_LH,:);
pVal_RH = pValues(chanOT_RH,:);
sigDiff_LH = tVal_LH.*double(pVal_LH<0.05);
sigDiff_RH = tVal_RH.*double(pVal_RH<0.05);
avg_sigDiff_LH = mean(sigDiff_LH);
avg_sigDiff_RH = mean(sigDiff_RH);



% figure(1)
% hold on
% plot(timeNew,avg_sigDiff_LH)
% plot(timeNew,avg_sigDiff_RH)
% plot(timeNew,(avg_sigDiff_LH-avg_sigDiff_RH))
% hold off



%% Max_positive and Max_negative sensors selection
% From "distributions of t-values", occipito-temporal sensors showing 
% Max(+) and Max(-) t-values in the time range 200-300ms were selected on
% each hemisphere


% LH
xl_LHchan = tVal_LH(:,timeN2pcIdx); % taking only the signals from LEFT occipitotemporal sensors (bcs this section is LH)
% xl_LHchan = tVal_LH;
xl_LHchan(xl_LHchan <= 0) = NaN;
[xl_LHchan,efflux_LHt] = max(xl_LHchan,[],2);
[efflux_LHval,efflux_LHChan] = max(xl_LHchan); 
efflux_LHt = efflux_LHt(efflux_LHChan);
efflux_LHt = timeN2pcVal(efflux_LHt);
% efflux_LHt = timeNew(efflux_LHt);

xl_LHchan = tVal_LH(:,timeN2pcIdx);
% xl_LHchan = tVal_LH;
xl_LHchan(xl_LHchan >= 0) = NaN; 
[xl_LHchan,influx_LHt] = min(xl_LHchan,[],2);
[influx_LHval,influx_LHChan] = min(xl_LHchan); 
influx_LHt = influx_LHt(influx_LHChan);
influx_LHt = timeN2pcVal(influx_LHt);
% influx_LHt = timeNew(influx_LHt);

% RH
xr_RHchan = tVal_RH(:,timeN2pcIdx); % taking only the signals from RIGHT occipitotemporal sensors (bcs this section is RH)
% xr_RHchan = tVal_RH;
xr_RHchan(xr_RHchan <= 0) = NaN;
[xr_RHchan,efflux_RHt] = max(xr_RHchan,[],2);
[efflux_RHval,efflux_RHChan] = max(xr_RHchan); 
efflux_RHt = efflux_RHt(efflux_RHChan);
efflux_RHt = timeN2pcVal(efflux_RHt);
% efflux_RHt = timeNew(efflux_RHt);

xr_RHchan = tVal_RH(:,timeN2pcIdx);
% xr_RHchan = tVal_RH;
xr_RHchan(xr_RHchan >= 0) = NaN; 
[xr_RHchan,influx_RHt] = min(xr_RHchan,[],2);
[influx_RHval,influx_RHChan] = min(xr_RHchan); 
influx_RHt = influx_RHt(influx_RHChan);
influx_RHt = timeN2pcVal(influx_RHt);
% influx_RHt = timeNew(influx_RHt);

effluxL = chanOT_LH(efflux_LHChan); 
influxL = chanOT_LH(influx_LHChan);
effluxR = chanOT_RH(efflux_RHChan);
influxR = chanOT_RH(influx_RHChan);
maxchan = [effluxL,influxL,effluxR,influxR];

%% 

% [effluxLeft,influxLeft,effluxRight,influxRight] = maxPmaxN(tValues,chanOT_LH,chanOT_RH,timeN2pcIdx);
% effluxL = chanOT_LH(effluxLeft); 
% influxL = chanOT_LH(influxLeft);
% effluxR = chanOT_RH(effluxRight);
% influxR = chanOT_RH(influxRight);



%% combination of efflux channel & influx channel (substraction)

avgTargetLVF = mean(targetLVF_meg,3);
avgTargetRVF = mean(targetRVF_meg,3);

LH_LVF = -avgTargetRVF(effluxL,:) - -avgTargetRVF(influxL,:);
LH_RVF = -avgTargetLVF(effluxL,:) - -avgTargetLVF(influxL,:);
RH_LVF = avgTargetRVF(effluxR,:) - avgTargetRVF(influxR,:);
RH_RVF = avgTargetLVF(effluxR,:) - avgTargetLVF(influxR,:);


% [~,pValues_LH,~,stats_LH] = ttest2(LH_LVF,LH_RVF,'Dim',1);
% tValues_LH = stats_LH.tstat;
% sigDiff_LH = tValues_LH.*double(pValues_LH<0.05);

% figure(3)
% hold on
% plot(timeNew,-avgTargetRVF(effluxL,:))
% plot(timeNew,-avgTargetRVF(influxL,:))
% plot(timeNew,(-avgTargetRVF(effluxL,:) - -avgTargetRVF(influxL,:)),'LineWidth',2)
% hold off
% 
% figure(4)
% hold on
% plot(timeNew,-avgTargetLVF(effluxL,:))
% plot(timeNew,-avgTargetLVF(influxL,:))
% plot(timeNew,(-avgTargetLVF(effluxL,:) - -avgTargetLVF(influxL,:)),'LineWidth',2)
% hold off
% 
% figure(7)
% hold on
% plot(timeNew,(avgTargetRVF(effluxL,:) - -avgTargetRVF(influxL,:)))
% plot(timeNew,(avgTargetLVF(effluxL,:) - -avgTargetLVF(influxL,:)))
% hold off
% 
% figure(8)
% hold on
% plot(timeNew,(avgTargetRVF(effluxR,:) - -avgTargetRVF(influxR,:)))
% plot(timeNew,(avgTargetLVF(effluxR,:) - -avgTargetLVF(influxR,:)))
% hold off
% 
% figure(5)
% hold on
% plot(timeNew,avgTargetRVF(effluxR,:))
% plot(timeNew,avgTargetRVF(influxR,:))
% hold off
% 
% figure(6)
% hold on
% plot(timeNew,avgTargetLVF(effluxR,:))
% plot(timeNew,avgTargetLVF(influxR,:))
% hold off

%%
% data_meg = bandpass(0.5,30,meg.srate, data_meg);
% LH_LVF = -data_meg(effluxL,:,idxTargetRVF) - -data_meg(influxL,:,idxTargetRVF);
% LH_RVF =  -data_meg(effluxL,:,idxTargetLVF) - -data_meg(influxL,:,idxTargetLVF);
% RH_LVF = data_meg(effluxR,:,idxTargetRVF) - data_meg(influxR,:,idxTargetRVF);
% RH_RVF = data_meg(effluxR,:,idxTargetLVF) - data_meg(influxR,:,idxTargetLVF);
% 
% LH_LVF = mean(LH_LVF,3);
% LH_RVF = mean(LH_RVF,3);
% RH_LVF = mean(RH_LVF,3);
% RH_RVF = mean(RH_RVF,3);

%% N2pc
N2pcLH = LH_LVF - LH_RVF;
N2pcRH = RH_LVF - RH_RVF;


figure(1)
hold on
plot(timeNew,LH_LVF)
plot(timeNew,LH_RVF)
plot(timeNew,N2pcLH,'k','LineWidth',2);
hold off

figure(2)
hold on
plot(timeNew,RH_LVF)
plot(timeNew,RH_RVF)
plot(timeNew,N2pcRH,'k','LineWidth',2);
hold off






% effinf = effluxMinInflux(data_meg,effluxL,influxL,effluxR,influxR,idxTargetLVF,idxTargetRVF);





%% check the N2pc amplitude
[N2pcLHpeak_val,N2pcLHpeak_id] = max(N2pcLH(:,timeN2pcIdx));
[N2pcRHpeak_val,N2pcRHpeak_id] = min(N2pcRH(:,timeN2pcIdx));

% time point where the N2pc peaked
N2pcLHpeak_pos = timeN2pcVal(N2pcLHpeak_id);
N2pcRHpeak_pos = timeN2pcVal(N2pcRHpeak_id);

%%
[max_positive,idx_max_positive] = max(sigDiff_LH(:,timeN2pcIdx));
[max_negative,idx_max_negative] = min(sigDiff_RH(:,timeN2pcIdx));

tpos_max_positive = timeN2pcVal(idx_max_positive);
tpos_max_negative = timeN2pcVal(idx_max_negative);

%%
erpdata_MEG.numLeftTrial = numLeftTrial;
erpdata_MEG.numRightTrial = numRightTrial;
erpdata_MEG.tValues = tValues;
erpdata_MEG.pValues = pValues;
erpdata_MEG.sigDiff = sigDiff;
erpdata_MEG.timeN2pcIdx = timeN2pcIdx;
erpdata_MEG.timeN2pcVal = timeN2pcVal;
erpdata_MEG.effluxL = effluxL;
erpdata_MEG.influxL = influxL;
erpdata_MEG.effluxR = effluxR;
erpdata_MEG.influxR = influxR;
erpdata_MEG.maxchan = maxchan;
erpdata_MEG.avgLeftHem_LVF = LH_LVF;
erpdata_MEG.avgLeftHem_RVF = LH_RVF;
erpdata_MEG.N2pcLH = N2pcLH;
erpdata_MEG.avgRightHem_LVF = RH_LVF;
erpdata_MEG.avgRightHem_RVF = RH_RVF;
erpdata_MEG.N2pcRH = N2pcRH;
erpdata_MEG.tVal_LH = tVal_LH;
erpdata_MEG.pVal_LH = pVal_LH;
erpdata_MEG.tVal_RH = tVal_RH;
erpdata_MEG.pVal_RH = pVal_RH;
erpdata_MEG.sigDiff_LH = sigDiff_LH;
erpdata_MEG.sigDiff_RH = sigDiff_RH;
erpdata_MEG.N2pcLHpeak_pos = N2pcLHpeak_pos;
erpdata_MEG.N2pcRHpeak_pos = N2pcRHpeak_pos;
erpdata_MEG.N2pcLHpeak_val = N2pcLHpeak_val;
erpdata_MEG.N2pcRHpeak_val = N2pcRHpeak_val;
erpdata_MEG.max_positive = max_positive;
erpdata_MEG.max_negative = max_negative;
erpdata_MEG.tpos_max_positive = tpos_max_positive;
erpdata_MEG.tpos_max_negative = tpos_max_negative;


% % LH
% xl = mean(targetLVF_meg,3); % averaging MEG signals when target is in LVF
% % xl = xl(:,timeN2pcIdx);
% 
% xl_LHchan = xl(chanOT_LH,:); % taking only the signals from LEFT occipitotemporal sensors (bcs this section is LH)
% xl_LHchan(xl_LHchan <= 0) = NaN;
% xl_LHchan = max(xl_LHchan,[],2);
% [effluxLHval,efflux_LHChan] = max(xl_LHchan); 
% 
% xl_LHchan = xl(chanOT_LH,:);
% xl_LHchan(xl_LHchan >= 0) = NaN; 
% xl_LHchan = min(xl_LHchan,[],2);
% [influxLHval,influx_LHChan] = min(xl_LHchan); 
% 
% 
% % RH
% xr = mean(targetRVF_meg,3); % averaging MEG signals when target is in RVF
% % xr = xr(:,timeN2pcIdx);
% 
% xr_RHchan = xr(chanOT_R,:); % taking only the signals from RIGHT occipitotemporal sensors (bcs this section is RH)
% xr_RHchan(xr_RHchan <= 0) = NaN;
% xr_RHchan = max(xr_RHchan,[],2);
% [effluxRHval,efflux_RHChan] = max(xr_RHchan); 
% 
% xr_RHchan = xr(chanOT_R,:);
% xr_RHchan(xr_RHchan >= 0) = NaN; 
% xr_RHchan = min(xr_RHchan,[],2);
% [influxRHval,influx_RHChan] = min(xr_RHchan); 


% figure(1)
% hold on
% plot(timeNew,(mean(targetLVF_meg(cL,:,:),3)))
% hold off
% 
% figure(2)
% hold on
% plot(timeNew,(mean(targetLVF_meg(9,:,:),3)))
% plot(timeNew,(mean(targetLVF_meg(19,:,:),3)))
% hold off
end