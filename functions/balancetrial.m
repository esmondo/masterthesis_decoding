function [Y,idxTargetLVF,idxTargetRVF,idxCorr] = balancetrial(Y)


idxTargetLVF = find(Y == 0); % indices of correct LVF response 
idxTargetRVF = find(Y == 1); % indices of correct RVF response 

numTargetLVF = sum(Y == 0);
numTargetRVF = sum(Y == 1);

if numTargetLVF >= numTargetRVF
    numBalance = numTargetRVF;
    numDiff = numTargetLVF - numTargetRVF;
    idxrandom = randperm(numBalance,numDiff);
    idxTargetLVF(idxrandom) =[];
%     idxTargetLVF(end-(numBalance-numDiff)) = [];
    idxTargetRVF = idxTargetRVF;
elseif numTargetLVF <= numTargetRVF
    numBalance = numTargetLVF;
    numDiff = numTargetRVF - numTargetLVF;
    idxrandom = randperm(numBalance,numDiff);
    idxTargetRVF(idxrandom) = [];
%     idxTargetRVF(end-(numBalance-numDiff)) = [];
    idxTargetLVF = idxTargetLVF;
end

idxCorr = unique([idxTargetLVF idxTargetRVF]);

Y = Y(idxCorr);

end