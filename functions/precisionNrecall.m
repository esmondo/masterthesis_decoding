function [precisionC1,recallC1,precisionC2,recallC2,accuracy,missclass] = precisionNrecall(Y,predictions)
% Y = label / actual
% P = predictions
% confComp = 


if size(Y,1) == 1
    Y = Y';
end

%%
confmat1 = strings(length(Y),size(predictions,2));
% confmat2 = strings(length(Y),size(predictions,2));

Ybig = repmat(Y,1,size(predictions,2));

confComp = predictions == Ybig;
% confComp = predictions == repmat(Y',1,size(predictions,2));



confmat1(confComp==1 & predictions==1) = 'TP';
confmat1(confComp==1 & predictions==0) = 'TN';
confmat1(confComp==0 & predictions==1) = 'FP';
confmat1(confComp==0 & predictions==0) = 'FN';

% confmat2(Ybig==1 & predictions==1) = 'TP';
% confmat2(Ybig==1 & predictions==0) = 'FN';
% confmat2(Ybig==0 & predictions==1) = 'FP';
% confmat2(Ybig==0 & predictions==0) = 'TN';



TP = sum(confmat1 == 'TP',1);
TN = sum(confmat1 == 'TN',1);
FP = sum(confmat1 == 'FP',1);
FN = sum(confmat1 == 'FN',1);

% TPx = sum(confmat2 == 'TP',1);
% TNx = sum(confmat2 == 'TN',1);
% FPx = sum(confmat2 == 'FP',1);
% FNx = sum(confmat2 == 'FN',1);

%% Left Presentation Relevant Trials = 0

precisionC1 = TN./(TN+FN);
recallC1 = TN./(TN+FP);

% precisionC1 = TNx./(TNx+FNx);
% recallC1 = TNx./(TNx+FPx);

%% Right Presentation Relevant Trials = 1

precisionC2 = TP./(TP+FP);
recallC2 = TP./(TP+FN); 

% precisionC2 = TPx./(TPx+FPx);
% recallC2 = TPx./(TPx+FNx); 

%% Accuracy

accuracy = (TP+TN)./(TP+TN+FP+FN);

% accuracy = (TPx+TNx)./(TPx+TNx+FPx+FNx);

%% Missclassification

missclass = (FP+FN)./(TP+TN+FP+FN);

% missclass = (FPx+FNx)./(TPx+TNx+FPx+FNx);



end