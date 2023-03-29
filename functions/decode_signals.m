function [valmetrics,valmetrics_stf,metrics_x,metrics_stf_x] = decode_signals(X_corID,Y_corID,X_badID,Y_badID,N,nreps,navg,artifactIdx)


%% classes
c1 = 0;
c2 = 1;


foldIdx = mod(1:length(Y_corID),N);
foldIDs = unique(foldIdx);

for k = 1:length(foldIDs)

    fprintf('Iteration %d \n',k);
    
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
             fprintf('     Repetition %d \n',r);

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
                              
                    
            X_test = X_corID(:,:,testIdx);
            Y_test = Y_corID(testIdx);

                % Without spatial filter                     
                Mdl{r} = trainClassifier(X_train,Y_train);
                predictions(1:length(testIdx),r) = testClassifier(Mdl{r}, X_test);
            
                % With spatial filter 
                stf{r} = estimateSTF(X_train, 2*(Y_train - 0.5), 0.05, 'impulse');
                
                filtDataTrain = applySTF(X_train,stf{r});

                filtDataTest = applySTF(X_test,stf{r});
                
                Mdl_stf{r} = trainClassifier(filtDataTrain,Y_train); 
                predictionsSTF(1:length(testIdx),r) = testClassifier(Mdl_stf{r},filtDataTest);                                         
        end            
    else
        
    %% with averaging the trials
            

        classSize = round(length(Y_corID(trainIdx))/2); % so that the classes will be balanced

        X_train = zeros(size(X_corID,1), size(X_corID,2), 2*classSize);
        Y_train = zeros(2*classSize,1);       
        X_test = zeros(size(X_corID,1), size(X_corID,2), length(testIdx));
        Y_test = Y_corID(testIdx);
        
        t1 = trainIdx(Y_corID(trainIdx) == c1); % c1 1st class label
        t2 = trainIdx(Y_corID(trainIdx) == c2); % c2 2nd class label
        
            for r = 1:nreps % repetitions  
                fprintf('     Repetition %d \n',r);
                
                %% Training Dataset
                                      
                    for c = 1:classSize  
                        
                        % first half classes are c1
                        ridx = randperm(length(t1));  % randomize number with the length of c1 trial for training 
                        X_train(:,:,c) = mean(X_corID(:,:,t1(ridx(1:navg))),3);  % it's taking the average of c1 trial (for training) from the randomized c1 trial index (ridx) as much as the predefined number of trial that will be averaged (navg)
                        Y_train(c) = c1;
                        
                        % second half classes are c2
                        ridx = randperm(length(t2));
                        X_train(:,:,c + classSize) = mean(X_corID(:,:,t2(ridx(1:navg))),3);  % it's taking the average of c2 trial (for training) from the randomized c2 trial index (ridx) as much as the predefined number of trial that will be averaged (navg)
                        Y_train(c + classSize) = c2;
                        
                    end
                           
                %% Validation Dataset
                Y_test = Y_corID(testIdx);
                        
                    for c = 1:length(testIdx) 
                        t3 = testIdx( Y_test == Y_test(c) ); % take only trials of current class from testIdx / testIdx will contain only trial with current classes (according to c)
                                               
                        ridx = randperm(length(t3));
                        X_test(:,:,c) = mean(X_corID(:,:,t3(ridx(1:navg))),3);

                    end
                        
                %  SVM Classification without STF
                Mdl{r} = trainClassifier(X_train,Y_train); % involves feature extraction and C calculation
                predictions(1:length(testIdx),r) = testClassifier(Mdl{r},X_test); 
                       
                %% applying STF
     
                    % estimate spatio-temporal filter
                    stf{r} = estimateSTF(X_train, 2*(Y_train-0.5), 0.05, 'impulse');
   
                    % apply spatial filter for training
                    filtDataTrain = applySTF(X_train,stf{r});
            
                    % apply spatial filter for validation
                    filtDataTest = applySTF(X_test,stf{r});
            
                    % SVM Classification applying STF
                    Mdl_stf{r} = trainClassifier(filtDataTrain,Y_train); % involves feature extraction and C calculation 
                    predictionsSTF(1:length(testIdx),r) = testClassifier(Mdl_stf{r},filtDataTest);  % predicting the test set using the SVM model
   
            end
                        
    end
    
    

    
    %% Evaluating Classification Model 
       
    validation(k) = metricsperf(predictions,Y_test);
    validation_stf(k) = metricsperf(predictionsSTF,Y_test);



     %% Testing incorrect responses
     
    bestmodel = [validation.bestmodel]; 
    pred_bad(:,k) = testClassifier(Mdl{bestmodel(k)},X_badID);
%     pred_bad(:,k) = testClassifier(Mdl,X_badID);
    
    bestmodel_stf = [validation_stf.bestmodel]; 
    filtDataTest_x = applySTF(X_badID,stf{bestmodel_stf(k)});
    pred_bad_stf(:,k) = testClassifier(Mdl_stf{bestmodel_stf(k)},filtDataTest_x);

    

    
    
    testmetrics_x(k) = metricsperf(pred_bad,Y_badID);
    testmetrics_stf_x(k) = metricsperf(pred_bad_stf,Y_badID);
    
end

valmetrics = sumperf(validation);
valmetrics_stf = sumperf(validation_stf);

metrics_x = sumperf(testmetrics_x);
metrics_stf_x = sumperf(testmetrics_stf_x);


end