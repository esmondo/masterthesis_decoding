function [X_sorted,y_sorted] = getRandTrials(X,y,categories)

% total number of trials
nTrials = size(X,1);   

rng('default')
randOrder = randperm(nTrials)';

%
randPred = X(randOrder,:);
randResp = y(randOrder,:);

X_sorted = [];
y_sorted = [];


    for ni = 1:size(categories,1)
    
        cat = categories(ni);
    
        catPred = randPred((randResp == cat),:);
        catResp = randResp((randResp == cat),:);
    
        X_sorted = [X_sorted; catPred];
        y_sorted = [y_sorted; catResp];

    end

end