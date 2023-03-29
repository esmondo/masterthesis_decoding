function [predictions,predictions_OT,accuracy,accuracy_OT,error,error_OT] = decodeN2pc(foldIdx,foldIDs,y)

for k = 1:length(foldIDs)
    
    trainIdx = foldIdx ~= foldIDs(k);
    testIdx = foldIdx == foldIDs(k);
    
    %% applying spatial filter
    
    % estimate spatio-temporal filter
    stf = estimateSTF(X(:,:,trainIdx), 2*(y(trainIdx)-0.5), 0.05, 'impulse');
    
    % apply spatial filter for training
    filtDataTrain = applySTF(X(:,:,trainIdx),stf);
    
    % train with SVM
    trainMdl = trainClassifier(filtDataTrain,y(trainIdx));
    
    % apply spatial filter for validation
    filtDataTest = applySTF(X(:,:,testIdx),stf);

    % predict the result of validation with the trained SVM
    predictions(testIdx) = testClassifier(trainMdl,filtDataTest);
    
%     % train
%     predictionsTrain = testClassifier(trainMdl,filtDataTrain);
%     predictionsTr(k) = sum(predictionsTrain == y(trainIdx)')/length(y(trainIdx)')*100;
    
    %% without spatial filter
    
    % train with SVM
    trainMdl_OT = trainClassifier(X(:,:,trainIdx),y(trainIdx));
    
    % predict the result of validation with the trained SVM
    predictions_OT(testIdx) = testClassifier(trainMdl_OT,X(:,:,testIdx));
    
end

accuracy = sum(predictions == y)/length(y)*100;
accuracy_OT = sum(predictions_OT == y)/length(y)*100;
error = (1-accuracy/100)*100;
error_OT = (1-accuracy_OT/100)*100;


end