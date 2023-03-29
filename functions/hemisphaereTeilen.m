function [Y,idxTargetLVF, idxTargetRVF] = hemisphaereTeilen(targetPosition_before)


Y(targetPosition_before <= 9) = 0; % LVF
Y(targetPosition_before > 9) = 1; % RVF


idxTargetLVF = find(Y == 0); % indices of correct LVF response 
idxTargetRVF = find(Y == 1); % indices of correct RVF response 



end