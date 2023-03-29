function [precisionC2,recallC2,precisionC1,recallC1] = confmat(Y,predictions)

confMat = confusionmat(Y',predictions,'Order',[1 0]);
TP = confMat(1,1);
FN = confMat(1,2);
FP = confMat(2,1);
TN = confMat(2,2);

%% Right Presentation Relevant Trials = 1

precisionC2 = TP/(TP+FP); % Precision / Positive Predictive Value
recallC2 = TP/(TP+FN); % Recall / True Positive Rate / Sensitivity

%% Left Presentation Relevant Trials = 0

precisionC1 = TN/(TN+FN); % Negative Predictive Value
recallC1 = TN/(TN+FP); % True Negative Rate / Specificity / Selectivity

end