function artFind = artReject(signal,threshold)

% peak-to-peak difference (max-min) across signals in every channel
artPeak = peak2peak(signal,2);

% logical values of artifacts greater than 6pT
artShow = artPeak > (threshold .* 10^-12);

% sum data row-wise to see which trial that has artifacts
artSum = sum(artShow);

% find indices of trials that are not zero (artifact detection)
artFind = find(artSum ~= 0).';



end