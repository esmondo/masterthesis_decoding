function megxcl = findArti_meg(data_meg)

[x1, x2, x3] = size(data_meg);

crit = 3e-12;

for ex = 1:x3
    tmp = find(abs(data_meg(:,:,ex)) > crit, 1);
    
    if isempty( tmp ) == 0
        megxcl( ex) = 1;
    else
        megxcl( ex ) = 0;
    end
    
end

%%
megxcl = logical( megxcl );

megabsxcl = find(megxcl == 1); % use the logical array to reject the unwanted trials / extract only clean trials

%%
data_meg(:,:,megxcl) = 0;

% v = squeeze(var(data_meg,0,[1 2]))';
v = squeeze(mean(var(data_meg,0,2)))';

crit_threshold = mean(v) + std(v)*4;

megxcl = [find(v > crit_threshold) megabsxcl]; 

end