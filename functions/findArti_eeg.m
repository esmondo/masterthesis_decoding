function eegxcl = findArti_eeg(data_eeg)

[x1, x2, x3] = size(data_eeg);

crit = 2e-4;
% crit = 1e-4;

for ex = 1:x3        
    tmp = find(abs(data_eeg(:,:,ex)) > crit, 1);
    
    if isempty( tmp ) == 0
        eegxcl( ex) = 1;
    else
        eegxcl( ex ) = 0;
    end
    
end

%%
eegxcl = logical( eegxcl );

eegabsxcl = find(eegxcl == 1); % use the logical array to reject the unwanted trials / extract only clean trials

%%
data_eeg(:,:,eegxcl) = 0;

% v = squeeze(var(data_eeg,0,[1 2]))';
v = squeeze(mean(var(data_eeg,0,2)))';

crit_threshold = mean(v) + std(v)*4;

eegxcl = [find(v > crit_threshold) eegabsxcl]; 

end

%         for ch = 1:x1
%             avg(ch) = mean(mean(squeeze(data_eeg(ch,:,:)),2));
%             sd(ch) = mean(std(squeeze(data_eeg(ch,:,:)),0,2));
%             
%             upthr = avg(ch) + std(ch)*4;
%             lowthr = avg(ch) - std(ch)*4;
%             
%             tmp1 = find(data_eeg(ch,:,ex) > upthr, 1);
%             tmp2 = find(data_eeg(ch,:,ex) < lowthr, 1);
%     
%             if isempty(tmp1) == 0
%                 eegxcl(ex) = 1;
%             elseif isempty(tmp2) == 0
%                 eegxcl(ex) = 1;
%             else
%                 eegxcl(ex) = 0;
%             end
%             
%         end