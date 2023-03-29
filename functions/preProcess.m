function [dataProcessed,labelProcessed] = preProcess(meg,EOGtrain,Xtrain,Ytrain)

% Xtrain = double(meg.data);
% dataEOG = double(eog.data);
time = meg.time;
Fs = meg.srate;
targetSide = meg.side;
subjectRes = meg.response;
targetPos = meg.position;
focusRate = meg.fokus;


%% Size of necessary data

[idxChan,~,~] = size(Xtrain);

%% Bandpass Filtering
   Xtrain = bandpass(1,200,Fs,Xtrain);
   EOGtrain = bandpass(1,30,Fs,EOGtrain);
      
%% Artifact Rejection 

  % Eye Movement Rejection
  Xtrain = rejectComponent(Xtrain,squeeze(EOGtrain(1,:,:))); % HEOG

  % Artficat Detection
  artFind = findArtifacts(Xtrain);
        
  % Artifact Rejection
  Xtrain(:,:,artFind) = [];
  targetSide(artFind) = []; % removing artifacts in side of targets's tilt
  subjectRes(artFind) = []; % removing artifacts in subject's response
  targetPos(artFind) = []; % removing artifacts in position where the target appeared
  focusRate(artFind) = []; % removing artifacts in participant's focus ratings
  Ytrain(artFind) = [];
  
  %% Downsample
  Fs_new = Fs;
  [timeDS,dataDS] = dwsample(idxChan,Fs_new,Fs,Xtrain,time);

%% baseline correction

Xtrain = baseline(dataDS,timeDS);

%% time range -200 to 800ms

% timeIdx = find(timeDS>=-0.2 & timeDS<=0.8);
% timeNew = timeDS(timeIdx);


%%
% [corrPos,~,~,~,~] = hemisphaereTeilen(data,corrTargPos);

%% preprocessed data

dataProcessed = Xtrain;
labelProcessed = Ytrain;

end