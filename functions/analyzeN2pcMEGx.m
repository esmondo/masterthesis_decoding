function erpdata_MEG = analyzeN2pcMEGx(meg,data_meg,chanOT_LH,chanOT_RH,corrTargPos,timeNew)

disp('Begin analyzing ERF...')

% time range 200 msec - 300 msec post stimulus
timeN2pcIdx = find(timeNew>=0.2 & timeNew<=0.3);
timeN2pcVal = timeNew(timeN2pcIdx);

%% data separation

% separate target by LVF and RVF 
[~,idxTargetLVF, idxTargetRVF] = hemisphaereTeilen(corrTargPos);
numLeftTrial = length(idxTargetLVF);
numRightTrial = length(idxTargetRVF);
erpdata_MEG.numLeftTrial = numLeftTrial;
erpdata_MEG.numRightTrial = numRightTrial;

% MEG left and MEG right
targetLVF_meg = data_meg(:,:,idxTargetLVF);
targetRVF_meg = data_meg(:,:,idxTargetRVF);


%% t-test

% compare targets in LVF vs targets in RVF
[~,pValues,~,stats] = ttest2(targetLVF_meg,targetRVF_meg,'Dim',3);

% the t-values
tValues = stats.tstat;

% t-values exceeding p-threshold (significantly different)
sigDiff = tValues.*double(pValues<0.05);

erpdata_MEG.tValues = tValues;
erpdata_MEG.pValues = pValues;
erpdata_MEG.sigDiff = sigDiff;

%% t-values of targets in both VF in each hemisphere

tVal_LH = tValues(chanOT_LH,:);
tVal_RH = tValues(chanOT_RH,:);

pVal_LH = pValues(chanOT_LH,:);
pVal_RH = pValues(chanOT_RH,:);

sigDiff_LH = tVal_LH.*double(pVal_LH<0.05);
sigDiff_RH = tVal_RH.*double(pVal_RH<0.05);


erpdata_MEG.tVal_LH = tVal_LH;
erpdata_MEG.pVal_LH = pVal_LH;
erpdata_MEG.tVal_RH = tVal_RH;
erpdata_MEG.pVal_RH = pVal_RH;
erpdata_MEG.sigDiff_LH = sigDiff_LH;
erpdata_MEG.sigDiff_RH = sigDiff_RH;



%% Max_positive and Max_negative sensors selection
% From "distributions of t-values", occipito-temporal sensors showing 
% Max(+) and Max(-) t-values in the time range 200-300ms were selected on
% each hemisphere

% LH
% maxpositive_LH = tValues(chanOT_LH,:);
maxpositive_LH = tValues(chanOT_LH,timeN2pcIdx); 
maxpositive_LH(maxpositive_LH <= 0) = NaN;
[maxpositive_LH,efflux_LHt] = max(maxpositive_LH,[],2);
[efflux_LHval,efflux_LHChan] = max(maxpositive_LH); 
efflux_LHt = efflux_LHt(efflux_LHChan);
% efflux_LHt = timeNew(efflux_LHt);
efflux_LHt = timeN2pcVal(efflux_LHt);
effluxL = chanOT_LH(efflux_LHChan); 

% maxnegative_LH = tValues(chanOT_LH,:);
maxnegative_LH = tValues(chanOT_LH,timeN2pcIdx);
maxnegative_LH(efflux_LHChan,:,:) = [];
maxnegative_LH(maxnegative_LH >= 0) = NaN; 
[maxnegative_LH,influx_LHt] = min(maxnegative_LH,[],2);
[influx_LHval,influx_LHChan] = min(maxnegative_LH); 
influx_LHt = influx_LHt(influx_LHChan);
% influx_LHt = timeNew(influx_LHt);
influx_LHt = timeN2pcVal(influx_LHt);
chanOT_LH(efflux_LHChan) = [];
influxL = chanOT_LH(influx_LHChan);

% RH
% maxpositive_RH = tValues(chanOT_RH,:); 
maxpositive_RH = tValues(chanOT_RH,timeN2pcIdx); 
maxpositive_RH(maxpositive_RH <= 0) = NaN;
[maxpositive_RH,efflux_RHt] = max(maxpositive_RH,[],2);
[efflux_RHval,efflux_RHChan] = max(maxpositive_RH); 
efflux_RHt = efflux_RHt(efflux_RHChan);
% efflux_RHt = timeNew(efflux_RHt);
efflux_RHt = timeN2pcVal(efflux_RHt);
effluxR = chanOT_RH(efflux_RHChan);

% maxnegative_RH = tValues(chanOT_RH,:);
maxnegative_RH = tValues(chanOT_RH,timeN2pcIdx);
maxnegative_RH(efflux_RHChan,:,:) = [];
maxnegative_RH(maxnegative_RH >= 0) = NaN; 
[maxnegative_RH,influx_RHt] = min(maxnegative_RH,[],2);
[influx_RHval,influx_RHChan] = min(maxnegative_RH); 
influx_RHt = influx_RHt(influx_RHChan);
% influx_RHt = timeNew(influx_RHt);
influx_RHt = timeN2pcVal(influx_RHt);
chanOT_RH(efflux_RHChan) = [];
influxR = chanOT_RH(influx_RHChan);


maxchan = [effluxL,influxL,effluxR,influxR];
erpdata_MEG.maxchan = maxchan;

maxchan_x = [efflux_LHChan, influx_LHChan, efflux_RHChan, influx_RHChan];
maxchan_xL = [efflux_LHChan, influx_LHChan];
maxchan_xR = [efflux_RHChan, influx_RHChan];
erpdata_MEG.maxchan_x = maxchan_x;

erpdata_MEG.efflux_LHChan = efflux_LHChan;
erpdata_MEG.influx_LHChan = influx_LHChan;
erpdata_MEG.efflux_RHChan = efflux_RHChan;
erpdata_MEG.influx_RHChan = influx_RHChan;

erpdata_MEG.effluxL = effluxL;
erpdata_MEG.influxL = influxL;
erpdata_MEG.effluxR = effluxR;
erpdata_MEG.influxR = influxR;

%%
avg_sigDiff_LH = mean(sigDiff_LH(maxchan_xL,:));
avg_sigDiff_RH = mean(sigDiff_RH(maxchan_xR,:));

erpdata_MEG.avg_sigDiff_LH = avg_sigDiff_LH;
erpdata_MEG.avg_sigDiff_RH = avg_sigDiff_RH;

%% combination of efflux channel & influx channel (substraction)

% LH_LVF = -targetLVF_meg(effluxL,:,:) - -targetLVF_meg(influxL,:,:);
% LH_RVF = -targetRVF_meg(effluxL,:,:) - -targetRVF_meg(influxL,:,:); %co
% RH_LVF = targetLVF_meg(effluxR,:,:) - targetLVF_meg(influxR,:,:); %co
% RH_RVF = targetRVF_meg(effluxR,:,:) - targetRVF_meg(influxR,:,:);

LH_LVF = targetLVF_meg(influxL,:,:) - targetLVF_meg(effluxL,:,:);
LH_RVF = targetRVF_meg(influxL,:,:) - targetRVF_meg(effluxL,:,:); %co
RH_LVF = targetLVF_meg(effluxR,:,:) - targetLVF_meg(influxR,:,:); %co
RH_RVF = targetRVF_meg(effluxR,:,:) - targetRVF_meg(influxR,:,:);

megcontra = cat(3,LH_RVF,RH_LVF);
megipsi = cat(3,LH_LVF,RH_RVF);

%%
[~,pVal_LH,~,stats_LH] = ttest2(LH_LVF,LH_RVF,'Dim',3);
tVal_LH = stats_LH.tstat;
sigDiff_LH = tVal_LH.*double(pVal_LH<0.05);
erpdata_MEG.tVal_LH = tVal_LH;
erpdata_MEG.sigDiff_LH = sigDiff_LH;

[max_LH,idx_max_LH] = min(sigDiff_LH(:,timeN2pcIdx));
tpos_max_LH = timeN2pcVal(idx_max_LH);
erpdata_MEG.max_LH_sdiff = max_LH;
erpdata_MEG.tpos_max_LH_sdiff = tpos_max_LH;

%
[~,pVal_RH,~,stats_RH] = ttest2(RH_LVF,RH_RVF,'Dim',3);
tVal_RH = stats_RH.tstat;
sigDiff_RH = tVal_RH.*double(pVal_RH<0.05);
erpdata_MEG.tVal_RH = tVal_RH;
erpdata_MEG.sigDiff_RH = sigDiff_RH;

[max_RH,idx_max_RH] = max(sigDiff_RH(:,timeN2pcIdx));
tpos_max_RH = timeN2pcVal(idx_max_RH);
erpdata_MEG.max_RH_sdiff = max_RH;
erpdata_MEG.tpos_max_RH_sdiff = tpos_max_RH;

%
[~,pVal_comb,~,stats_comb] = ttest2(megcontra,megipsi,'Dim',3);
tVal_comb = stats_comb.tstat;
sigDiff_comb = tVal_comb.*double(pVal_comb<0.05);
erpdata_MEG.tVal_comb = tVal_comb;
erpdata_MEG.sigDiff_comb = sigDiff_comb;

[max_comb,idx_max_comb] = max(sigDiff_comb(:,timeN2pcIdx));
tpos_max_comb = timeN2pcVal(idx_max_comb);
erpdata_MEG.max_conip_sdiff = max_comb;
erpdata_MEG.tpos_max_conip_sdiff = tpos_max_comb;

% figure(1);
% hold on
% plot(timeNew,sigDiff_LH,'b');
% plot(timeNew,sigDiff_RH,'r--');
% plot(timeNew,sigDiff_comb,'k','LineWidth',2);
% hold off

%% N2pc (for plotting)

avgLH_LVF = mean(LH_LVF,3);
avgLH_RVF = mean(LH_RVF,3);
avgRH_LVF = mean(RH_LVF,3);
avgRH_RVF = mean(RH_RVF,3);

avgLH_LVF = lowpass(30,meg.srate, avgLH_LVF);
avgLH_RVF = lowpass(30,meg.srate, avgLH_RVF);
avgRH_LVF = lowpass(30,meg.srate, avgRH_LVF);
avgRH_RVF = lowpass(30,meg.srate, avgRH_RVF);

erpdata_MEG.avgLH_LVF = avgLH_LVF;
erpdata_MEG.avgLH_RVF = avgLH_RVF;
erpdata_MEG.avgRH_LVF = avgRH_LVF;
erpdata_MEG.avgRH_RVF = avgRH_RVF;

%
avgmegcontra = mean(megcontra,3);
avgmegipsi = mean(megipsi,3);

avgmegcontra = lowpass(30,meg.srate, avgmegcontra);
avgmegipsi = lowpass(30,meg.srate, avgmegipsi);

erpdata_MEG.avgmegcontra = avgmegcontra;
erpdata_MEG.avgmegipsi = avgmegipsi;


% collapse LVF-RVF (N2pc)
N2pcLH = avgLH_LVF - avgLH_RVF;
N2pcRH = avgRH_LVF - avgRH_RVF;
N2pc_conip = avgmegcontra - avgmegipsi;
erpdata_MEG.N2pcLH = N2pcLH;
erpdata_MEG.N2pcRH = N2pcRH;
erpdata_MEG.N2pc_conip = N2pc_conip;


% figure(2);
% hold on
% plot(timeNew,avgLH_LVF,'b');
% plot(timeNew,avgLH_RVF,'r--');
% plot(timeNew,N2pcLH,'k','LineWidth',2);
% hold off
% 
% figure(3);
% hold on
% plot(timeNew,avgRH_LVF,'b');
% plot(timeNew,avgRH_RVF,'r--');
% plot(timeNew,N2pcRH,'k','LineWidth',2);
% hold off
% 
% figure(4);
% hold on
% plot(timeNew,avgmegcontra,'b');
% plot(timeNew,avgmegipsi,'r--');
% plot(timeNew,N2pc_conip,'k','LineWidth',2);
% hold off


%% N2pc peak (amplitude)
[max_N2pcLH,idx_max_N2pcLH] = max(N2pcLH(:,timeN2pcIdx));
tpos_max_N2pcLH = timeN2pcVal(idx_max_N2pcLH);
erpdata_MEG.max_LH_N2pcpeak = max_N2pcLH;
erpdata_MEG.tpos_max_LH_N2pcpeak = tpos_max_N2pcLH;

[max_N2pcRH,idx_max_N2pcRH] = max(N2pcRH(:,timeN2pcIdx));
tpos_max_N2pcRH = timeN2pcVal(idx_max_N2pcRH);
erpdata_MEG.max_RH_N2pcpeak = max_N2pcRH;
erpdata_MEG.tpos_max_RH_N2pcpeak = tpos_max_N2pcRH;

[max_N2pc_conip,idx_max_N2pc_conip] = max(N2pc_conip(:,timeN2pcIdx));
tpos_max_N2pc_conip = timeN2pcVal(idx_max_N2pc_conip);
erpdata_MEG.max_conip_N2pcpeak = max_N2pc_conip;
erpdata_MEG.tpos_max_conip_N2pcpeak = tpos_max_N2pc_conip;


end