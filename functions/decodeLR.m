function [predictions,predictionsSTF,alg1,alg2] = decodeLR(Xraw,Y,N,navg,nreps,meg,eog,corID)


X = double(Xraw).*1e12;
X = bandpass(1,200,meg.srate,X);

Xeog = double(eog.data(1,:,corID));
Xeog = bandpass(1,30,meg.srate,Xeog);

Xrej = rejectComponent(double(Xraw),squeeze(Xeog)); % HEOG
artFind = findArtifacts(Xrej);
Xrej = Xrej.*1e12;
Xrej = bandpass(1,200,meg.srate,Xrej);


% % time of interest: -200 - 800ms
% timeCutIdx = find(meg.time>=-0.2 & meg.time<=0.1);
% timeCut = meg.time(timeCutIdx);
% X = X(:,timeCutIdx,:);
% Xrej = Xrej(:,timeCutIdx,:);       

% lowpass filter 12.5 Hz
X = lowpass(12.5,meg.srate,X);
Xrej = lowpass(12.5,meg.srate,Xrej);
        
% DownSample 50 Hz
FsNew = 50; 
[time_X_DS,X_DS] = dwsample(size(X,1),FsNew,meg.srate,X,meg.time); %%%!! change time
[time_Xrej_DS,Xrej_DS] = dwsample(size(Xrej,1),FsNew,meg.srate,Xrej,meg.time); %%%!!

% Baseline Correction
X = baseline(X_DS,time_X_DS);
Xrej = baseline(Xrej_DS,time_Xrej_DS);

% time of interest: 0 - 800ms
timeCutIdx1 = find(time_X_DS>=0 & time_X_DS<=0.8);
timeCutIdx2 = find(time_Xrej_DS>=0 & time_Xrej_DS<=0.8);
% timeCut = meg.time(timeCutIdx);
X = X(:,timeCutIdx1,:);
Xrej = Xrej(:,timeCutIdx2,:); 

%%
foldIdx = mod(1:length(Y),N); % just one possibility to define the folds
foldIDs = unique(foldIdx);
predictions = zeros(size(Y,2),nreps);
predictionsSTF = zeros(size(Y,2),nreps);

c1 = 0;
c2 = 1;


%% 10-fold Cross-Validation

    for k = 1:length(foldIDs)
 
    %%       
%     trainIdx = find(foldIdx~=foldIDs(k)); % you could opt to balance the class sizes if quite unbalanced
    trainIdx = setdiff(find(foldIdx~=foldIDs(k)),artFind);
     
    testIdx = find(foldIdx==foldIDs(k));

    
    %%
        if navg > 1
        %% with averaging the trials
            
        classSize = round(length(Y(trainIdx))/2); % just one possibility to set class sizes
        Xtrain = zeros(size(Xrej,1), size(Xrej,2), 2*classSize);
        Ytrain = zeros(2*classSize,1);
        Xtest = zeros(size(X,1), size(X,2), length(testIdx));
        t1 = trainIdx(Y(trainIdx) == c1); % c1 1st class label
        t2 = trainIdx(Y(trainIdx) == c2); % c2 2nd class label
        
        for r = 1:nreps % repetitions            
            %% Training Dataset
            for c = 1:classSize  
                
                ridx = randperm(length(t1));
                Xtrain(:,:,c) = mean(Xrej(:,:,t1(ridx(1:navg))),3);
                Ytrain(c) = c1;
                
                ridx = randperm(length(t2));
                Xtrain(:,:,c + classSize) = mean(Xrej(:,:,t2(ridx(1:navg))),3);
                Ytrain(c + classSize) = c2;
            end
            
                
            %% Validation Dataset
            for c = 1:length(testIdx)
                
                t3 = testIdx(Y(testIdx)==Y(testIdx(c))); % take only trials of current class label
                ridx = randperm(length(t3));
                Xtest(:,:,c) = mean(X(:,:,t3(ridx(1:navg))),3);               
            end
            


                    %%  SVM Classification without STF
                    alg1 = trainClassifier(Xtrain,Ytrain); % involves feature extraction and C calculation
                    predictions(testIdx,r) = testClassifier(alg1,Xtest); 
                       
                    %% applying STF
     
                    % estimate spatio-temporal filter
                    stf = estimateSTF(Xtrain, 2*(Ytrain-0.5), 0.05, 'impulse');
   
                    % apply spatial filter for training
                    filtDataTrain = applySTF(Xtrain,stf);
            
                    % apply spatial filter for validation
                    filtDataTest = applySTF(Xtest,stf);
            
                    %% SVM Classification applying STF
                    alg2 = trainClassifier(filtDataTrain,Ytrain); % involves feature extraction and C calculation 
                    predictionsSTF(testIdx,r) = testClassifier(alg2,filtDataTest);  % predicting the test set using the SVM model
            
        end
        
    else
        %% Single-trial: Training Dataset
      
        Xtrain = Xrej(:,:,trainIdx);   
        Ytrain = Y(trainIdx);        

        %% Single-trial: Validation Dataset
        
        Xtest = X(:,:,testIdx);
     
                %% SVM Classification without STF
                alg1 = trainClassifier(Xtrain,Ytrain); % involves feature extraction and C calculation
                predictions(testIdx) = testClassifier(alg1, Xtest);
            
                %% applying spatial filter
                stf = estimateSTF(Xtrain, 2*(Ytrain-0.5), 0.05, 'impulse');
                filtDataTrain = applySTF(Xtrain,stf);
                filtDataTest = applySTF(Xtest,stf);
            
                %% SVM Classification applying STF
                alg2 = trainClassifier(filtDataTrain,Ytrain); % involves feature extraction and C calculation
                predictionsSTF(testIdx) = testClassifier(alg2,filtDataTest); 
        
        end  
    
    end

end