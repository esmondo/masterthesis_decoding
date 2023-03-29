% function [avgLeftHem_LVF, avgLeftHem_RVF, avgRightHem_LVF, avgRightHem_RVF,...
%     tVal_LH,pVal_LH, tVal_RH,pVal_RH] = effluxMinInflux...
%     (data,effluxLeft,influxLeft,effluxRight,influxRight,idxTargetLVF,idxTargetRVF)
    
function effinf = effluxMinInflux...
    (data_meg,effluxLeft,influxLeft,effluxRight,influxRight,idxTargetLVF,idxTargetRVF)


%% LEFT HEMISPHERE

% (Max_positive - Max_negative) for targets in LVF
leftHem_LVF = -data_meg(effluxLeft,:,idxTargetRVF) - -data_meg(influxLeft,:,idxTargetRVF); %influx-effkux
avgLeftHem_LVF = mean(leftHem_LVF,3);

% (Max_positive - Max_negative) for targets in RVF
leftHem_RVF = -data_meg(effluxLeft,:,idxTargetLVF) - -data_meg(influxLeft,:,idxTargetLVF);
avgLeftHem_RVF = mean(leftHem_RVF,3);

%% RIGHT HEMISPHERE

% (Max_positive - Max_negative) for targets in LVF
rightHem_LVF = data_meg(effluxRight,:,idxTargetRVF) - data_meg(influxRight,:,idxTargetRVF);
avgRightHem_LVF = mean(rightHem_LVF,3);

% (Max_positive - Max_negative) for targets in RVF
rightHem_RVF = data_meg(effluxRight,:,idxTargetLVF) - data_meg(influxRight,:,idxTargetLVF);
avgRightHem_RVF = mean(rightHem_RVF,3);

%% compare targets in LVF vs targets in RVF
[~,pVal_LH,~,stats_LH] = ttest2(leftHem_LVF,leftHem_RVF,'Dim',3);

% the t-values
tVal_LH = stats_LH.tstat;

% t-values exceeding p-threshold (significantly different)
sigDiff_LH = tVal_LH.*double(pVal_LH<0.05);


% 
[~,pVal_RH,~,stats_RH] = ttest2(rightHem_LVF,rightHem_RVF,'Dim',3);

% the t-values
tVal_RH = stats_RH.tstat;

% t-values exceeding p-threshold (significantly different)
sigDiff_RH = tVal_RH.*double(pVal_RH<0.05);

%%
effinf.avgLeftHem_LVF = avgLeftHem_LVF;
effinf.avgLeftHem_RVF = avgLeftHem_RVF;
effinf.avgRightHem_LVF = avgRightHem_LVF;
effinf.avgRightHem_RVF = avgRightHem_RVF;
effinf.tVal_LH = tVal_LH;
effinf.pVal_LH = pVal_LH;
effinf.tVal_RH = tVal_RH;
effinf.pVal_RH = pVal_RH;
effinf.sigDiff_LH = sigDiff_LH;
effinf.sigDiff_RH = sigDiff_RH;

%%
% % LEFT HEMISPHERE
% 
% % (Max_positive - Max_negative) for targets in LVF
% leftHem_LVF = data(influxLeft,:,idxTargetLVF) - data(effluxLeft,:,idxTargetLVF);
% avgLeftHem_LVF = mean(leftHem_LVF,3);
% 
% % (Max_positive - Max_negative) for targets in RVF
% leftHem_RVF = data(influxLeft,:,idxTargetRVF) - data(effluxLeft,:,idxTargetRVF);
% avgLeftHem_RVF = mean(leftHem_RVF,3);



% % RIGHT HEMISPHERE
% 
% % (Max_positive - Max_negative) for targets in LVF
% rightHem_LVF = data(effluxRight,:,idxTargetLVF) - data(influxRight,:,idxTargetLVF);
% avgRightHem_LVF = mean(rightHem_LVF,3);
% 
% % (Max_positive - Max_negative) for targets in RVF
% rightHem_RVF = data(effluxRight,:,idxTargetRVF) - data(influxRight,:,idxTargetRVF);
% avgRightHem_RVF = mean(rightHem_RVF,3);






% % LEFT HEMISPHERE
% 
% % (Max_positive - Max_negative) for targets in LVF
% leftHem_LVF = -(data(effluxLeft,:,idxTargetRVF) - data(influxLeft,:,idxTargetRVF)); %influx-effkux
% avgLeftHem_LVF = mean(leftHem_LVF,3);
% 
% % (Max_positive - Max_negative) for targets in RVF
% leftHem_RVF = -(data(effluxLeft,:,idxTargetLVF) - data(influxLeft,:,idxTargetLVF));
% avgLeftHem_RVF = mean(leftHem_RVF,3);
% % 
% % 
% % 
% % RIGHT HEMISPHERE
% 
% % (Max_positive - Max_negative) for targets in LVF
% rightHem_LVF = data(effluxRight,:,idxTargetRVF) - data(influxRight,:,idxTargetRVF);
% avgRightHem_LVF = mean(rightHem_LVF,3);
% 
% % (Max_positive - Max_negative) for targets in RVF
% rightHem_RVF = data(effluxRight,:,idxTargetLVF) - data(influxRight,:,idxTargetLVF);
% avgRightHem_RVF = mean(rightHem_RVF,3);
end