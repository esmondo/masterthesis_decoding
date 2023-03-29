function eegxcl = findArtifacts_EEG(data_eeg)

[x1, x2, x3] = size(data_eeg);

crit = 2e-4;

for ex = 1:x3
    tmp = find(data_eeg(:,:,ex) > crit, 1);
    
    if isempty( tmp ) == 0
        eegxcl( ex) = 1;
    else
        eegxcl( ex ) = 0;
    end
    
end

%%
eegxcl = logical( eegxcl );

eegabsxcl = find(eegxcl == 1); % use the logical array to reject the unwanted trials / extract only clean trials

data_eeg(:,:,eegxcl) = 0;

% v = squeeze(var(data_eeg,0,[1 2]))';
v = squeeze(mean(var(data_eeg,0,2)))';

crit_threshold = mean(v) + std(v)*4;

eegxcl = [find(v > crit_threshold) eegabsxcl]; 

end