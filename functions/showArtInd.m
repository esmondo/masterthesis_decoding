function idxArt = showArtInd(thresholdNum,signal)

% artifact detection (threshold 3 pT)

% threshold value
thrUpper = thresholdNum * 10^-12; % tesla to pico-tesla
thrLower = -(thresholdNum) * 10^-12;

% setting upper limit and lower limit for the threshold
% resulting logical values based on the upper and lower limit
artUpper = signal >= thrUpper;
artLower = signal <= thrLower;

% combine them, so we have logical values of detected artifacts
% that need to be removed
artAll = artUpper + artLower;

% indexing the artifacts:
% (1) summing the sample points (matrix row-wise), so we are able to 
% mark the artifacts in each channel
sumPoints = sum(artAll,2);
sumPoints(sumPoints > 0) = 1;

% (2) since we have to remove any trials consisting artifacts, we need to
% sum the channels and mark them
sumChan = sum(sumPoints);
sumChan(sumChan > 0) = 1;

% (3) using the information above, the artifacts indices are obtained:
idxArt = find(sumChan == 1).';

end