function [accuracyTrain,accuracyTest,errorTrain,errorTest,...
    dataTrain,trainLabel,dataTest,testLabel] = ...
    crossvalidate(corrPos,dataCls,chanOcciTemp)

rng('default')
% partition = cvpartition(label','KFold',10);
partition = cvpartition(corrPos,'Holdout',0.1);
% partition = cvpartition(label','LeaveOut');

for ii = 1:partition.NumTestSets
    trainIdx = find(training(partition,ii) == 1);
    testIdx = find(test(partition,ii) == 1);
    
    trainLabel = corrPos(:,trainIdx)';
    testLabel = corrPos(:,testIdx)';
end

% estimate spatio-temporal filters for training set
stf = estimateSTF(dataCls(chanOcciTemp,:,trainIdx), 2*(corrPos(trainIdx)'-0.5), 0.05, 'impulse');

% Spatially filtered train data
filteredDataTrain = applySTF(dataCls(chanOcciTemp,:,trainIdx),stf);  

dataTrain = permute(filteredDataTrain,[3 2 1]);

    if size(dataTrain,3) > 1
        dataTrain = mean(dataTrain,3);
    end
    
C = 1/mean(diag(dataTrain*dataTrain'));


% Spatially filtered test data
filteredDataTest = applySTF(dataCls(chanOcciTemp,:,testIdx),stf);
dataTest = permute(filteredDataTest,[3 2 1]);

    if size(dataTest,3) > 1
        dataTest = mean(dataTest,3);
    end

% SVM Model
Mdl = fitcsvm(dataTrain,trainLabel,...
    'KernelFunction', 'linear',...
    'KernelScale', 1,...
    'BoxConstraint',C);

% Train Accuracy
errorTrain = resubLoss(Mdl);
accuracyTrain = 100*(1-errorTrain);

% CV Accuracy
cvMdl = crossval(Mdl);
cvTrainError = kfoldLoss(cvMdl);
cvTrainAccuracy = 1-cvTrainError;

% Test Accuracy
errorTest = loss(Mdl,dataTest,testLabel);
accuracyTest = 100*(1-errorTest);


end

