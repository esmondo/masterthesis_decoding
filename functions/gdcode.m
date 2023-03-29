% function [valmetrics,valmetrics_stf,metrics,metrics_stf,metrics_x,metrics_stf_x] = gdcode(X_corID,Y_corID,X_badID,Y_badID,N,nreps,navg,artifactIdx,idx)
function [valmetrics,valmetrics_stf,metrics_x,metrics_stf_x] = gdcode(X_corID,Y_corID,X_badID,Y_badID,N,nreps,navg,artifactIdx,idx)

%% 
% corID = idx.corID;
% badID = idx.badID;
% 
% m = numel(Y_corID);
% P = 0.8;
% 
% idxpart = randperm(m);
% 
% % Train and Test data
% idxTrain = idxpart(1:round(P*m));
% idxTest = idxpart(round(P*m)+1:end);
% 
% X4train = X_corID(:,:,idxTrain);
% Y4train = Y_corID(idxTrain);
% 
% X4test = X_corID(:,:,idxTest);
% Y4test = Y_corID(idxTest);

%% classes
c1 = 0;
c2 = 1;

% foldIdx = mod(1:length(Y4train),N); % 80:20
foldIdx = mod(1:length(Y_corID),N);
foldIDs = unique(foldIdx);

for k = 1:length(foldIDs)

    
    %% Subset fold indices  
    
    if isempty(artifactIdx) == 0
        trainIdx = setdiff(find(foldIdx~=foldIDs(k)),artifactIdx);
    else
        trainIdx = find(foldIdx~=foldIDs(k));
    end
        
    testIdx = find(foldIdx==foldIDs(k));
    
    predictions = zeros(size(testIdx,2),nreps);
    predictionsSTF = zeros(size(testIdx,2),nreps);
    
    
    %% Single-trial: Training Dataset
    
    if navg == 1
        
 
        for r = 1:nreps % repetitions
%             X_train = X4train(:,:,trainIdx); % 80:20   
%             Y_train = Y4train(trainIdx); 
            X_train = X_corID(:,:,trainIdx);   
            Y_train = Y_corID(trainIdx);  
            
            % balance trial
            numc1 = sum(Y_train==c1);
            numc2 = sum(Y_train==c2);
            if numc1 > numc2
                d = numc1 - numc2;
                bb = find(Y_train==c1);
                f = randsample(bb,d); 
                Y_train(f) = [];
                X_train(:,:,f) = [] ;
            elseif numc2 > numc1
                d = numc2 - numc1;
                bb = find(Y_train==c2);
                f = randsample(bb,d);
                Y_train(f) = [];
                X_train(:,:,f) = [] ;
            end
            
            
            numc1_balanced = sum(Y_train==c1);
            numc2_balanced = sum(Y_train==c2);
                              
                    
%             X_test = X4train(:,:,testIdx); % 80:20
%             Y_test = Y4train(testIdx);
            X_test = X_corID(:,:,testIdx);
            Y_test = Y_corID(testIdx);

                % Without spatial filter                     
                Mdl = trainClassifier(X_train,Y_train);
                predictions(1:length(testIdx),r) = testClassifier(Mdl, X_test);
            
                % With spatial filter 
                stf = estimateSTF(X_train, 2*(Y_train - 0.5), 0.05, 'impulse');
                filtDataTrain = applySTF(X_train,stf);
                filtDataTest = applySTF(X_test,stf);
                
                Mdl_stf = trainClassifier(filtDataTrain,Y_train); 
                predictionsSTF(1:length(testIdx),r) = testClassifier(Mdl_stf,filtDataTest);                                         
        end            
    else
        
    %% with averaging the trials
            
%         classSize = round(length(Y4train(trainIdx))/2); % so that the classes will be balanced (80:20)
        classSize = round(length(Y_corID(trainIdx))/2); % so that the classes will be balanced

%         X_train = zeros(size(X4train,1), size(X4train,2), 2*classSize); % 80:20
%         Y_train = zeros(2*classSize,1);       
%         X_test = zeros(size(X4train,1), size(X4train,2), length(testIdx));
%         Y_test = Y4train(testIdx);
        X_train = zeros(size(X_corID,1), size(X_corID,2), 2*classSize);
        Y_train = zeros(2*classSize,1);       
        X_test = zeros(size(X_corID,1), size(X_corID,2), length(testIdx));
        Y_test = Y_corID(testIdx);
        
%         t1 = trainIdx(Y4train(trainIdx) == c1); % c1 1st class label (80:20)
%         t2 = trainIdx(Y4train(trainIdx) == c2); % c2 2nd class label
        t1 = trainIdx(Y_corID(trainIdx) == c1); % c1 1st class label
        t2 = trainIdx(Y_corID(trainIdx) == c2); % c2 2nd class label
        
            for r = 1:nreps % repetitions            
                %% Training Dataset
                                      
                    for c = 1:classSize  
                        
                        % first half classes are c1
                        ridx = randperm(length(t1));  % randomize number with the length of c1 trial for training 
                        X_train(:,:,c) = mean(X_corID(:,:,t1(ridx(1:navg))),3);  % it's taking the average of c1 trial (for training) from the randomized c1 trial index (ridx) as much as the predefined number of trial that will be averaged (navg)
%                         X_train(:,:,c) = mean(X4train(:,:,t1(ridx(1:navg))),3); % 80:20
                        Y_train(c) = c1;
                        
                        % second half classes are c2
                        ridx = randperm(length(t2));
                        X_train(:,:,c + classSize) = mean(X_corID(:,:,t2(ridx(1:navg))),3);  % it's taking the average of c2 trial (for training) from the randomized c2 trial index (ridx) as much as the predefined number of trial that will be averaged (navg)
%                         X_train(:,:,c + classSize) =
%                         mean(X4train(:,:,t2(ridx(1:navg))),3); % 80:20
                        Y_train(c + classSize) = c2;
                        
                    end
                           
                %% Validation Dataset
                Y_test = Y_corID(testIdx);
%                 Y_test = Y4train(testIdx); % 80:20
%                 numc1 = sum(Y_test==c1);
%                 numc2 = sum(Y_test==c2);
                        
                    for c = 1:length(testIdx) 
                        t3 = testIdx( Y_test == Y_test(c) ); % take only trials of current class from testIdx / testIdx will contain only trial with current classes (according to c)
                        
                        
%                         if length(t3) < navg
%                             oo = Y_corID == Y_test(c);
%                             
%                             xxx = Y_badID == Y_badID(testIdx_x(c));
%                             ygbeda = setdiff(t3,foldIdx_x(xxx));
%                             xtra = randsample(ygbeda,navg-length(t3));                           
% %                             xtra = randsample(t3,navg-length(t3));
%                             t3 = cat(2,t3,xtra);
%                             ridx = randperm(length(t3));  
%                         end
                        
                        
                        ridx = randperm(length(t3));
                        X_test(:,:,c) = mean(X_corID(:,:,t3(ridx(1:navg))),3);
%                         X_test(:,:,c) = mean(X4train(:,:,t3(ridx(1:navg))),3); % 80:20
                    end
                        
                %  SVM Classification without STF
                Mdl = trainClassifier(X_train,Y_train); % involves feature extraction and C calculation
                predictions(1:length(testIdx),r) = testClassifier(Mdl,X_test); 
                       
                %% applying STF
     
                    % estimate spatio-temporal filter
                    stf = estimateSTF(X_train, 2*(Y_train-0.5), 0.05, 'impulse');
   
                    % apply spatial filter for training
                    filtDataTrain = applySTF(X_train,stf);
            
                    % apply spatial filter for validation
                    filtDataTest = applySTF(X_test,stf);
            
                    % SVM Classification applying STF
                    Mdl_stf = trainClassifier(filtDataTrain,Y_train); % involves feature extraction and C calculation 
                    predictionsSTF(1:length(testIdx),r) = testClassifier(Mdl_stf,filtDataTest);  % predicting the test set using the SVM model
   
            end
                        
    end
    
    %% Testing 20% data that has been trained for validation
%     pred_testsplit(:,k) = testClassifier(Mdl,X4test);
%     
%     filtDataTestSplit = applySTF(X4test,stf);
%     pred_testsplit_stf(:,k) = testClassifier(Mdl_stf,filtDataTestSplit);
    
     %% Testing incorrect responses
    pred_bad(:,k) = testClassifier(Mdl,X_badID);
    
    filtDataTest_x = applySTF(X_badID,stf);
    pred_bad_stf(:,k) = testClassifier(Mdl_stf,filtDataTest_x);
    
%% Evaluating Classification Model 
       
validation(k) = metricsperf(predictions,Y_test);
validation_stf(k) = metricsperf(predictionsSTF,Y_test);

% testmetrics(k) = metricsperf(pred_testsplit,Y4test);
% testmetrics_stf(k) = metricsperf(pred_testsplit_stf,Y4test);
    
testmetrics_x(k) = metricsperf(pred_bad,Y_badID);
testmetrics_stf_x(k) = metricsperf(pred_bad_stf,Y_badID);
    
end

valmetrics = sumperf(validation);
valmetrics_stf = sumperf(validation_stf);

% metrics = sumperf(testmetrics);
% metrics_stf = sumperf(testmetrics_stf);

metrics_x = sumperf(testmetrics_x);
metrics_stf_x = sumperf(testmetrics_stf_x);



%     perf = mean([evalmetrics.perf]).*100;
%     accuracy = mean([evalmetrics.accuracy]).*100;
%     perf_std = std([evalmetrics.accuracy]).*100;
%     missclass = mean([evalmetrics.missclass]).*100;
%     recall_c1 = mean([evalmetrics.recall_c1]).*100;
%     recall_c2 = mean([evalmetrics.recall_c2]).*100;
%     precision_c2 = mean([evalmetrics.precision_c2]).*100;
%     precision_c1 = mean([evalmetrics.precision_c1]).*100;
%     f1 = mean([evalmetrics.f1]).*100;
%     tpr = mean([evalmetrics.tpr]).*100;
%     fpr = mean([evalmetrics.fpr]).*100;
% 
%     perf_stf = mean([evalmetrics_stf.perf]).*100;
%     accuracy_stf = mean([evalmetrics_stf.accuracy]).*100;
%     perf_std_stf = std([evalmetrics_stf.accuracy]).*100;
%     missclass_stf = mean([evalmetrics_stf.missclass]).*100;
%     recall_c1_stf = mean([evalmetrics_stf.recall_c1]).*100;
%     recall_c2_stf = mean([evalmetrics_stf.recall_c2]).*100;
%     precision_c2_stf = mean([evalmetrics_stf.precision_c2]).*100;
%     precision_c1_stf = mean([evalmetrics_stf.precision_c1]).*100;
%     f1_stf = mean([evalmetrics_stf.f1]).*100;
%     tpr_stf = mean([evalmetrics_stf.tpr]).*100;
%     fpr_stf = mean([evalmetrics_stf.fpr]).*100;
% 
%  
%     metrics.perf = perf;
%     metrics.accuracy = accuracy;
%     metrics.perf_std = perf_std;
%     metrics.missclass = missclass;
%     metrics.recall_c1 = recall_c1;
%     metrics.recall_c2 = recall_c2;
%     metrics.precision_c2 = precision_c2;
%     metrics.precision_c1 = precision_c1;
%     metrics.f1 = f1;
%     metrics.tpr = tpr;
%     metrics.fpr = fpr;
% 
%     metrics_stf.perf = perf_stf;
%     metrics_stf.accuracy = accuracy_stf;
%     metrics_stf.perf_std = perf_std_stf;
%     metrics_stf.missclass = missclass_stf;
%     metrics_stf.recall_c1 = recall_c1_stf;
%     metrics_stf.recall_c2 = recall_c2_stf;
%     metrics_stf.precision_c2 = precision_c2_stf;
%     metrics_stf.precision_c1 = precision_c1_stf;
%     metrics_stf.f1 = f1_stf;
%     metrics_stf.tpr = tpr_stf;
%     metrics_stf.fpr = fpr_stf;

    
    %%
%     perf_x = mean([evalmetrics_x.perf]).*100;
%     accuracy_x = mean([evalmetrics_x.accuracy]).*100;
%     perf_std_x = std([evalmetrics_x.accuracy]).*100;
%     missclass_x = mean([evalmetrics_x.missclass]).*100;
%     recall_c1_x = mean([evalmetrics_x.recall_c1]).*100;
%     recall_c2_x = mean([evalmetrics_x.recall_c2]).*100;
%     precision_c2_x = mean([evalmetrics_x.precision_c2]).*100;
%     precision_c1_x = mean([evalmetrics_x.precision_c1]).*100;
%     f1_x = mean([evalmetrics_x.f1]).*100;
%     tpr_x = mean([evalmetrics_x.tpr]).*100;
%     fpr_x = mean([evalmetrics_x.fpr]).*100;
% 
%     perf_stf_x = mean([evalmetrics_stf_x.perf]).*100;
%     accuracy_stf_x = mean([evalmetrics_stf_x.accuracy]).*100;
%     perf_std_stf_x = std([evalmetrics_stf_x.accuracy]).*100;
%     missclass_stf_x = mean([evalmetrics_stf_x.missclass]).*100;
%     recall_c1_stf_x = mean([evalmetrics_stf_x.recall_c1]).*100;
%     recall_c2_stf_x = mean([evalmetrics_stf_x.recall_c2]).*100;
%     precision_c2_stf_x = mean([evalmetrics_stf_x.precision_c2]).*100;
%     precision_c1_stf_x = mean([evalmetrics_stf_x.precision_c1]).*100;
%     f1_stf_x = mean([evalmetrics_stf_x.f1]).*100;
%     tpr_stf_x = mean([evalmetrics_stf_x.tpr]).*100;
%     fpr_stf_x = mean([evalmetrics_stf_x.fpr]).*100;
% 
%  
%     metrics_x.perf = perf_x;
%     metrics_x.accuracy = accuracy_x;
%     metrics_x.perf_std = perf_std_x;
%     metrics_x.missclass = missclass_x;
%     metrics_x.recall_c1 = recall_c1_x;
%     metrics_x.recall_c2 = recall_c2_x;
%     metrics_x.precision_c2 = precision_c2_x;
%     metrics_x.precision_c1 = precision_c1_x;
%     metrics_x.f1 = f1_x;
%     metrics_x.tpr = tpr_x;
%     metrics_x.fpr = fpr_x;
% 
%     metrics_stf_x.perf = perf_stf_x;
%     metrics_stf_x.accuracy = accuracy_stf_x;
%     metrics_stf_x.perf_std = perf_std_stf_x;
%     metrics_stf_x.missclass = missclass_stf_x;
%     metrics_stf_x.recall_c1 = recall_c1_stf_x;
%     metrics_stf_x.recall_c2 = recall_c2_stf_x;
%     metrics_stf_x.precision_c2 = precision_c2_stf_x;
%     metrics_stf_x.precision_c1 = precision_c1_stf_x;
%     metrics_stf_x.f1 = f1_stf_x;
%     metrics_stf_x.tpr = tpr_stf_x;
%     metrics_stf_x.fpr = fpr_stf_x;

end