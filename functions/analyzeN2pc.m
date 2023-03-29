function [tValues,pValues,sigDiff,...
    numLeftTrial,numRightTrial,...
    timeN2pcIdx,timeN2pcVal,...
    effluxL,influxL,effluxR,influxR,maxchan,...
    avgLeftHem_LVF,avgLeftHem_RVF,N2pcLH_sub,...
    avgRightHem_LVF,avgRightHem_RVF,N2pcRH_sub,...
    tVal_LH,pVal_LH, tVal_RH,pVal_RH,...
    N2pcLHpeak_pos,N2pcRHpeak_pos,...
    N2pcLHpeak_val,N2pcRHpeak_val] = analyzeN2pcMEG(meg,data,cL,cR,corrTargPos,timeNew)



%% data separation

% separate target by LVF and RVF 
[~,idxTargetLVF, idxTargetRVF] = hemisphaereTeilen(corrTargPos);
numLeftTrial = length(idxTargetLVF);
numRightTrial = length(idxTargetRVF);

% MEG left and MEG right
targetLVF_meg = data(:,:,idxTargetLVF);
targetRVF_meg = data(:,:,idxTargetRVF);





%% t-test

% compare targets in LVF vs targets in RVF
[~,pVal,~,stats] = ttest2(targetLVF_meg,targetRVF_meg,'Dim',3);

% the t-values
tValues = stats.tstat;

% the p-values
pValues = pVal;

% t-values exceeding p-threshold (significantly different)
sigDiff = tValues.*double(pValues<0.05);


%% Max_positive and Max_negative sensors selection
% From "distributions of t-values", occipito-temporal sensors showing 
% Max(+) and Max(-) t-values in the time range 200-300ms were selected on
% each hemisphere


% time range 200 msec - 300 msec post stimulus
timeN2pcIdx = find(timeNew>=0.198 & timeNew<=0.302);
timeN2pcVal = timeNew(timeN2pcIdx);

[effluxLeft,influxLeft,effluxRight,influxRight] = ...
    maxPmaxN(tValues,cL,cR,timeN2pcIdx);

effluxL = cL(effluxLeft); 
influxL = cL(influxLeft);
effluxR = cR(effluxRight);
influxR = cR(influxRight);

maxchan = [effluxL,influxL,effluxR,influxR];



%% combination of efflux channel & influx channel (substraction)


[avgLeftHem_LVF, avgLeftHem_RVF, avgRightHem_LVF, avgRightHem_RVF,...
   tVal_LH,pVal_LH, tVal_RH,pVal_RH] = effluxMinInflux...
    (data,effluxL,influxL,effluxR,influxR,idxTargetLVF,idxTargetRVF);



%% smoothering signals for plot

avgLeftHem_LVF = bandpass(0.5,30,meg.srate,avgLeftHem_LVF);
avgLeftHem_RVF = bandpass(0.5,30,meg.srate,avgLeftHem_RVF);
N2pcLH_sub = avgLeftHem_LVF-avgLeftHem_RVF;
avgRightHem_LVF = bandpass(0.5,30,meg.srate,avgRightHem_LVF);
avgRightHem_RVF = bandpass(0.5,30,meg.srate,avgRightHem_RVF);
N2pcRH_sub = avgRightHem_LVF-avgRightHem_RVF;



%%

% sub = 7;
% 
% % LH
% leftFig = figure('Name','Left Hemisphere');
% hold on
% plot(timeNew,avgLeftHem_LVF(:,:,sub),'b');
% plot(timeNew,avgLeftHem_RVF(:,:,sub),'r--');
% plot(timeNew,N2pcLH_sub(:,:,sub),'k','LineWidth',2)
% hold off
% xlabel('time')
% ylabel('amplitude')
% title('Left Hemisphere')
% legend('LVF Target','RVF Target','LVF-RVF')
% ylim([-2.5e-13 2.5e-13])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');
% line([0.18 0.18],ylim,'Color','k','LineStyle','--');
% line([0.3 0.3],ylim,'Color','k','LineStyle','--');
% 
% 
% % RH
% rightFig = figure('Name','Right Hemisphere');
% hold on
% plot(timeNew,avgRightHem_LVF(:,:,sub),'b');
% plot(timeNew,avgRightHem_RVF(:,:,sub),'r--');
% plot(timeNew,N2pcRH_sub(:,:,sub),'k','LineWidth',2)
% hold off
% xlabel('time')
% ylabel('amplitude')
% title('Right Hemisphere')
% legend('LVF Target','RVF Target','LVF-RVF')
% ylim([-2.5e-13 2.5e-13])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');
% line([0.18 0.18],ylim,'Color','k','LineStyle','--');
% line([0.3 0.3],ylim,'Color','k','LineStyle','--');


%% check the N2pc amplitude
[N2pcLHpeak_val,N2pcLHpeak_id] = max(N2pcLH_sub(:,timeN2pcIdx));
[N2pcRHpeak_val,N2pcRHpeak_id] = min(N2pcRH_sub(:,timeN2pcIdx));



%% time point where the N2pc peaked
N2pcLHpeak_pos = timeN2pcVal(N2pcLHpeak_id);
N2pcRHpeak_pos = timeN2pcVal(N2pcRHpeak_id);



end