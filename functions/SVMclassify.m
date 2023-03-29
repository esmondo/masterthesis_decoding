function [errorTrain,accuracyTrain,errorTest,accuracyTest,dataset] = SVMclassify(allData,corrPos)
    
%
X = permute(allData,[3 2 1]);
y = corrPos';
C = 1/mean(diag(X*X'));

%
dataset = array2table(X);
dataset.labels = y;

predictors = dataset(:,1:end-1);
response = dataset.labels;

%
rng('default')
n = length(response);
hpartition = cvpartition(response,'Holdout',0.1);
idxTrain = training(hpartition);
tblTrain = dataset(idxTrain,:);

idxTest = test(hpartition);
tblTest = dataset(idxTest,:);
% testData = table2array(tblTest);


%
Mdl = fitcsvm(tblTrain,'labels',...
    'KernelFunction', 'linear',...
    'KernelScale', 1,...
    'BoxConstraint',C);

%
errorTrain = resubLoss(Mdl);
accuracyTrain = 100*(1-errorTrain);

cvMdl = crossval(Mdl);
cvTrainError = kfoldLoss(cvMdl);
cvTrainAccuracy = 1-cvTrainError;

errorTest = loss(Mdl,tblTest, 'labels');
accuracyTest = 100*(1-errorTest);

end