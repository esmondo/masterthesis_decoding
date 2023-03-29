function [predictions,predictionsSTF,alg1,alg2] = decode_meg(Xmeg,Y,Z,N,navg,nreps,meg)
% X is the predictor (meg.data)
% Y is the label (target position)
% Z is the component hat needs to be rejected (eog.data)

disp(['Begin MEG decoding algorithm of navg ', num2str(navg)])

% heog to be rejected (if necessary)
heog = double(Z(1,:,:));
heog = bandpass(1,30,meg.srate,heog);

% X to be trained
Xrej = double(Xmeg);
Xrej = bandpass(1,200,meg.srate,Xrej);
% Xrej = rejectComponent(Xrej,squeeze(heog)); % HEOG
artFind = findArti_meg(Xrej);
Xrej = Xrej.*1e12;

% X to be tested
Xori = double(Xmeg);
Xori = bandpass(1,200,meg.srate,Xori);
Xori = Xori.*1e12;

% lowpass filter 12.5 Hz
Xrej = lowpass(12.5,meg.srate,Xrej);
Xori = lowpass(12.5,meg.srate,Xori);

% DownSample 50 Hz
FsNew = 50; 
[ds_Xrej,t_Xrej] = dwsample(size(Xrej,1),FsNew,meg.srate,Xrej,meg.time); 
[ds_Xori,t_Xori] = dwsample(size(Xori,1),FsNew,meg.srate,Xori,meg.time); 

% Baseline Correction
Xrej = baseline(ds_Xrej,t_Xrej);
Xori = baseline(ds_Xori,t_Xori);

% time of interest: 0 - 800ms
timeCutIdx = find(t_Xori>=0 & t_Xori<=0.8);

% timeCut = meg.time(timeCutIdx);
Xrej = Xrej(:,timeCutIdx,:); 
Xori = Xori(:,timeCutIdx,:);


%%
foldIdx = mod(1:length(Y),N); % just one possibility to define the folds
foldIDs = unique(foldIdx);
predictions = zeros(size(Y,2),nreps);
predictionsSTF = zeros(size(Y,2),nreps);

% classes
c1 = 0;
c2 = 1;


%% 10-fold Cross-Validation

    for k = 1:length(foldIDs)
        disp(['Begin 10-fold cross validation: Fold ', num2str(k)])
 
    %%       
%     trainIdx = find(foldIdx~=foldIDs(k)); % you could opt to balance the class sizes if quite unbalanced
    trainIdx = setdiff(find(foldIdx~=foldIDs(k)),artFind);
     
    testIdx = find(foldIdx==foldIDs(k));

    
    %%
        if navg > 1
        %% with averaging the trials
            
        classSize = round(length(Y(trainIdx))/2); % just one possibility to set class sizes
                                                  % the length of label for training divided by 2
        
        Xtrain = zeros(size(Xrej,1), size(Xrej,2), 2*classSize);
        Ytrain = zeros(2*classSize,1);
        Xtest = zeros(size(Xori,1), size(Xori,2), length(testIdx));
        
        t1 = trainIdx(Y(trainIdx) == c1); % c1 1st class label
        t2 = trainIdx(Y(trainIdx) == c2); % c2 2nd class label
        
            for r = 1:nreps % repetitions            
                %% Training Dataset
                
                    for c = 1:classSize 
                
                        % randomize number with the length of c1 trial for training
                        ridx = randperm(length(t1)); 
                
                        % it's taking the average of c1 trial (for training) 
                        % from the randomized c1 trial index (ridx) 
                        % as much as the predefined number of trial that will be averaged (navg)
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
                        Xtest(:,:,c) = mean(Xori(:,:,t3(ridx(1:navg))),3);               
                    end
            
                        %  SVM Classification without STF
                        alg1 = trainClassifier(Xtrain,Ytrain); % involves feature extraction and C calculation
                        predictions(testIdx,r) = testClassifier(alg1,Xtest); 
                       
                %% applying STF
     
                    % estimate spatio-temporal filter
                    stf = estimateSTF(Xtrain, 2*(Ytrain-0.5), 0.05, 'impulse');
   
                    % apply spatial filter for training
                    filtDataTrain = applySTF(Xtrain,stf);
            
                    % apply spatial filter for validation
                    filtDataTest = applySTF(Xtest,stf);
            
                    % SVM Classification applying STF
                    alg2 = trainClassifier(filtDataTrain,Ytrain); % involves feature extraction and C calculation 
                    predictionsSTF(testIdx,r) = testClassifier(alg2,filtDataTest);  % predicting the test set using the SVM model
            
            end
        
        else                
               %% Single-trial: Training Dataset
               
               for r = 1:nreps % repetitions
      
                    Xtrain = Xrej(:,:,trainIdx);   
                    Ytrain = Y(trainIdx);        

                    %% Single-trial: Validation Dataset
        
                    Xtest = Xori(:,:,testIdx);
     
                    % SVM Classification without STF
                    alg1 = trainClassifier(Xtrain,Ytrain); % involves feature extraction and C calculation
                    predictions(testIdx,r) = testClassifier(alg1, Xtest);
            
                    %% applying spatial filter
                    stf = estimateSTF(Xtrain, 2*(Ytrain-0.5), 0.05, 'impulse');
                    filtDataTrain = applySTF(Xtrain,stf);
                    filtDataTest = applySTF(Xtest,stf);
            
                    % SVM Classification applying STF
                    alg2 = trainClassifier(filtDataTrain,Ytrain); % involves feature extraction and C calculation
                    predictionsSTF(testIdx,r) = testClassifier(alg2,filtDataTest); 
               end
        
        end  
    
    end

end