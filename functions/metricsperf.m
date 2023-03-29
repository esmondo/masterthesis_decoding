function metrics = metricsperf(predictions,Y)

if size(Y,1) == 1
    Y = Y';
end

numOfSamplePred = length(Y);

Ytest = repmat(Y,1,size(predictions,2));

%% ACCURACY
% perf = sum(predictions == repmat(Y,1,nreps))./length(Y)*100;
% perf = sum(predictions == Ybig)./length(Y);
perf =  mean(sum(predictions == Ytest)./length(Y));

    

%% CONFUSION MATRIX
confmat = strings(length(Y),size(predictions,2));

confmat(Ytest==1 & predictions==1) = 'TP';
confmat(Ytest==1 & predictions==0) = 'FN';
confmat(Ytest==0 & predictions==1) = 'FP';
confmat(Ytest==0 & predictions==0) = 'TN';

TP = sum(confmat == 'TP',1);
TN = sum(confmat == 'TN',1);
FP = sum(confmat == 'FP',1);
FN = sum(confmat == 'FN',1);


%% Accuracy
accuracy = (TP+TN)./(TP+TN+FP+FN);
[bestacc,bestmodel] = max(accuracy);

interval = 1.96.*sqrt((accuracy.*(1-accuracy))./length(Y)); % CI 95%

% se = std(accuracy) ./ sqrt(length(accuracy));
% ts = tinv([0.025 0.975],length(accuracy)-1);
% CI = mean(accuracy) + ts.*se;

%% Missclassification
missclass = (FP+FN)./(TP+TN+FP+FN);

%% Recall 
recall_c1 = TN./(TN+FP); % c1=0 / specificity / TNR
recall_c2 = TP./(TP+FN); % c2=1 / sensitivity / TPR

%% Precision
precision_c1 = TN./(TN+FN); % c1=0
precision_c2 = TP./(TP+FP); % c2=1

%% F1 measure
f1_score = 2.*TP ./ (2.*TP + FP + FN);

%% ROC Curve
tpr = TP./(TP+FN); % TPR
fpr = FP./(FP+TN); % FPR


%%
metrics.bestmodel = bestmodel;
metrics.bestacc = bestacc;
metrics.perf = perf;
metrics.accuracy = mean(accuracy);
metrics.missclass = mean(missclass);
metrics.recall_c1 = mean(recall_c1);
metrics.recall_c2 = mean(recall_c2);
metrics.precision_c2 = mean(precision_c2);
metrics.precision_c1 = mean(precision_c1);
metrics.f1 = mean(f1_score);
metrics.tpr = mean(tpr);
metrics.fpr = mean(fpr);
metrics.interval = interval;
metrics.TP = TP;
metrics.TN = TN;
metrics.FP = FP;
metrics.FN = FN;
metrics.numOfSamplePred = numOfSamplePred;




%%

% perf = mean(perf,2);
% perf_std = std(accuracy);
% precisionC1 = mean(precisionC1,2);
% recallC1 = mean(recallC1,2);
% precisionC2 = mean(precisionC2,2);
% recallC2 = mean(recallC2,2);
% accuracy = mean(accuracy,2);
% missclass = mean(missclass,2);

%% combine all in a result struct
% metrics.perf = perf .*100;
% % metrics.perf_std = perf_std .* 100;
% metrics.accuracy = accuracy .*100;
% metrics.precision_c1 = precisionC1 .*100;
% metrics.recall_c1 = recallC1 .*100;
% metrics.precision_c2 = precisionC2 .*100;
% metrics.recall_c2 = recallC2 .*100;
% metrics.missclassification = missclass .*100;



end