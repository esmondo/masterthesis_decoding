function [metrics,metrics_stf] = goDecode_x(X_cor,X_bad,Y_corID,Y_badID,N,nreps,navg,artifactIdx,idx)

corID = idx.corID;
badID = idx.badID;

% N*2 for 80:20 partitions
% N*3 for 70:30 partitions
foldIdx = mod(1:length(Y_corID),N);
foldIDs = unique(foldIdx);

foldIdx_x = mod(1:length(Y_badID),N);
foldIDs_x = unique(foldIdx_x);

% if whichresponse == 'corID'
%     predictions = zeros(size(Y_badID,2),nreps);
%     predictionsSTF = zeros(size(Y_badID,2),nreps);
% elseif whichresponse == 'badID'
%     predictions = zeros(size(Xori,3),nreps);
%     predictionsSTF = zeros(size(Xori,3),nreps);
% end

% classes
c1 = 0;
c2 = 1;


%% 10-fold Cross-Validation

    for k = 1:length(foldIDs)
        
        disp(['Begin 10-fold cross validation: Fold ', num2str(k)])
 
    %% CV Index     
    
    if isempty(artifactIdx) == 0
        trainIdx = setdiff(find(foldIdx~=foldIDs(k)),artifactIdx);
    else
        trainIdx = find(foldIdx~=foldIDs(k)); % you could opt to balance the class sizes if quite unbalanced
    end
     
    testIdx = find(foldIdx==foldIDs(k));   
    testIdx_x = find(foldIdx_x==foldIDs_x(k));
    
    predictions = zeros(size(testIdx_x,2),nreps);
    predictionsSTF = zeros(size(testIdx_x,2),nreps);
    
    %%
        if navg > 1
        %% with averaging the trials
            
        classSize = round(length(Y_corID(trainIdx))/2); % just one possibility to set class sizes
                                                  % the length of label for training divided by 2
        
        Xtrain = zeros(size(X_cor,1), size(X_cor,2), 2*classSize);
        Ytrain = zeros(2*classSize,1);
        Xtest = zeros(size(X_bad,1), size(X_bad,2), length(testIdx_x));
        Ytest = Y_badID(testIdx_x);
        
        t1 = trainIdx(Y_corID(trainIdx) == c1); % c1 1st class label
        t2 = trainIdx(Y_corID(trainIdx) == c2); % c2 2nd class label
        
            for r = 1:nreps % repetitions            
                %% Training Dataset
                
                   for c = 1:classSize                                     
                        ridx = randperm(length(t1));  % randomize number with the length of c1 trial for training 
                        Xtrain(:,:,c) = mean(X_cor(:,:,t1(ridx(1:navg))),3);  % it's taking the average of c1 trial (for training) from the randomized c1 trial index (ridx) as much as the predefined number of trial that will be averaged (navg)
                        Ytrain(c) = c1;
                
                        ridx = randperm(length(t2));
                        Xtrain(:,:,c + classSize) = mean(X_cor(:,:,t2(ridx(1:navg))),3);  % it's taking the average of c2 trial (for training) from the randomized c2 trial index (ridx) as much as the predefined number of trial that will be averaged (navg)
                        Ytrain(c + classSize) = c2;
                    end
                           
                %% Validation Dataset
                    for c = 1:length(testIdx_x)                
                        t3 = testIdx_x(Y_badID(testIdx_x)==Y_badID(testIdx_x(c))); % take only trials of current class label
                        
                         if length(t3) < navg                            
                            xxx = Y_badID == Y_badID(testIdx_x(c));
                            ygbeda = setdiff(t3,foldIdx_x(xxx));
                            xtra = randsample(ygbeda,navg-length(t3));                           
%                             xtra = randsample(t3,navg-length(t3));
                            t3 = cat(2,t3,xtra);
                            ridx = randperm(length(t3));  
                         end
                        
                        ridx = randperm(length(t3));
                        Xtest(:,:,c) = mean(X_bad(:,:,t3(ridx(1:navg))),3);               
                    end
            
                        %  SVM Classification without STF
                        alg1 = trainClassifier(Xtrain,Ytrain); % involves feature extraction and C calculation
                        predictions(1:length(testIdx_x),r) = testClassifier(alg1,Xtest); 
                       
                %% applying STF
     
                    % estimate spatio-temporal filter
                    stf = estimateSTF(Xtrain, 2*(Ytrain-0.5), 0.05, 'impulse');
   
                    % apply spatial filter for training
                    filtDataTrain = applySTF(Xtrain,stf);
            
                    % apply spatial filter for validation
                    filtDataTest = applySTF(Xtest,stf);
            
                    % SVM Classification applying STF
                    alg2 = trainClassifier(filtDataTrain,Ytrain); % involves feature extraction and C calculation 
                    predictionsSTF(1:length(testIdx_x),r) = testClassifier(alg2,filtDataTest);  % predicting the test set using the SVM model
            
            end
        
        else                
               %% Single-trial: Training Dataset
               
               for r = 1:nreps % repetitions
      
                    Xtrain = X_cor(:,:,trainIdx);   
                    Ytrain = Y_corID(trainIdx);  
                    
                    Xtest = X_bad(:,:,testIdx_x);
                    Ytest =  Y_badID(testIdx_x);

                    %% Without spatial filter                 
     
                    % SVM Classification without STF
                    alg1 = trainClassifier(Xtrain,Ytrain); % involves feature extraction and C calculation
                    predictions(1:length(testIdx_x),r) = testClassifier(alg1, Xtest);
            
                    %% With spatial filter
                    stf = estimateSTF(Xtrain, 2*(Ytrain-0.5), 0.05, 'impulse');
                    filtDataTrain = applySTF(Xtrain,stf);
                    filtDataTest = applySTF(Xtest,stf);
            
                    % SVM Classification applying STF
                    alg2 = trainClassifier(filtDataTrain,Ytrain); % involves feature extraction and C calculation
                    predictionsSTF(1:length(testIdx_x),r) = testClassifier(alg2,filtDataTest);  
                    
               end       
        end  
        
        %% Evaluating Classification Model 
       
        evalmetrics(k) = metricsperf(predictions,Ytest);
        evalmetrics_stf(k) = metricsperf(predictionsSTF,Ytest);
    
    end
    
    perf = mean([evalmetrics.perf]).*100;
    accuracy = mean([evalmetrics.accuracy]).*100;
    perf_std = std([evalmetrics.accuracy]).*100;
    missclass = mean([evalmetrics.missclass]).*100;
    recall_c1 = mean([evalmetrics.recall_c1]).*100;
    recall_c2 = mean([evalmetrics.recall_c2]).*100;
    precision_c2 = mean([evalmetrics.precision_c2]).*100;
    precision_c1 = mean([evalmetrics.precision_c1]).*100;
    f1 = mean([evalmetrics.f1]).*100;
    tpr = mean([evalmetrics.tpr]).*100;
    fpr = mean([evalmetrics.fpr]).*100;

    perf_stf = mean([evalmetrics_stf.perf]).*100;
    accuracy_stf = mean([evalmetrics_stf.accuracy]).*100;
    perf_std_stf = std([evalmetrics_stf.accuracy]).*100;
    missclass_stf = mean([evalmetrics_stf.missclass]).*100;
    recall_c1_stf = mean([evalmetrics_stf.recall_c1]).*100;
    recall_c2_stf = mean([evalmetrics_stf.recall_c2]).*100;
    precision_c2_stf = mean([evalmetrics_stf.precision_c2]).*100;
    precision_c1_stf = mean([evalmetrics_stf.precision_c1]).*100;
    f1_stf = mean([evalmetrics_stf.f1]).*100;
    tpr_stf = mean([evalmetrics_stf.tpr]).*100;
    fpr_stf = mean([evalmetrics_stf.fpr]).*100;
    
    
    
    metrics.perf = perf;
    metrics.accuracy = accuracy;
    metrics.perf_std = perf_std;
    metrics.missclass = missclass;
    metrics.recall_c1 = recall_c1;
    metrics.recall_c2 = recall_c2;
    metrics.precision_c2 = precision_c2;
    metrics.precision_c1 = precision_c1;
    metrics.f1 = f1;
    metrics.tpr = tpr;
    metrics.fpr = fpr;

    metrics_stf.perf = perf_stf;
    metrics_stf.accuracy = accuracy_stf;
    metrics_stf.perf_std = perf_std_stf;
    metrics_stf.missclass = missclass_stf;
    metrics_stf.recall_c1 = recall_c1_stf;
    metrics_stf.recall_c2 = recall_c2_stf;
    metrics_stf.precision_c2 = precision_c2_stf;
    metrics_stf.precision_c1 = precision_c1_stf;
    metrics_stf.f1 = f1_stf;
    metrics_stf.tpr = tpr_stf;
    metrics_stf.fpr = fpr_stf;
    
end