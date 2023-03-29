function results = calculate_performance(predictions,Y)

if size(Y,1) == 1
    Y = Y';
end

Ybig = repmat(Y,1,size(predictions,2));

%% ACCURACY
% perf = sum(predictions == repmat(Y,1,nreps))./length(Y)*100;
perf = sum(predictions == Ybig)./length(Y);

%% CONFUSION MATRIX
confmat = strings(length(Y),size(predictions,2));

confmat(Ybig==1 & predictions==1) = 'TP';
confmat(Ybig==1 & predictions==0) = 'FN';
confmat(Ybig==0 & predictions==1) = 'FP';
confmat(Ybig==0 & predictions==0) = 'TN';

TP = sum(confmat == 'TP',1);
TN = sum(confmat == 'TN',1);
FP = sum(confmat == 'FP',1);
FN = sum(confmat == 'FN',1);

% Left Presentation Relevant Trials = 0
precisionC1 = TN./(TN+FN);
recallC1 = TN./(TN+FP);

% Right Presentation Relevant Trials = 1
precisionC2 = TP./(TP+FP);
recallC2 = TP./(TP+FN); 

% Accuracy
accuracy = (TP+TN)./(TP+TN+FP+FN);

% Missclassification
missclass = (FP+FN)./(TP+TN+FP+FN);

%%

perf = mean(perf,2);
perf_std = std(accuracy);
precisionC1 = mean(precisionC1,2);
recallC1 = mean(recallC1,2);
precisionC2 = mean(precisionC2,2);
recallC2 = mean(recallC2,2);
accuracy = mean(accuracy,2);
missclass = mean(missclass,2);

%% combine all in a result struct
results.perf = perf .*100;
results.perf_std = perf_std .* 100;
results.precision_c1 = precisionC1 .*100;
results.recall_c1 = recallC1 .*100;
results.precision_c2 = precisionC2 .*100;
results.recall_c2 = recallC2 .*100;
results.accuracy = accuracy .*100;
results.missclassification = missclass .*100;


% results.std = perf_std;

end