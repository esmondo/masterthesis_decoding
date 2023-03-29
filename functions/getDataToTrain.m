function [X_train, y_train, X_test, y_test] = ...
    getDataToTrain(X_sorted,y_sorted,categories,vecsPerCat,foldSizes,iterationNum)

% pre-allocate data    
X_train = [];
y_train = [];
X_test = [];
y_test = [];


% Get the number of folds from the foldSizes matrix.
numFolds = size(foldSizes,2); % 2x10

catStart = 1;

% For each category (LVF & RVF)
for (catIndex = 1 : size(categories,1))

    % Get the list of fold sizes for this category as a column vector.
    catFoldSizes = foldSizes(catIndex, :)';
    
    % Set the starting index of the first fold for this category.
    foldStart = catStart;
    
    % For each fold (in each category 0 & 1)
    for (foldIndex = 1 : numFolds)
        
        % Compute the index of the last vector in this fold.
        foldEnd = foldStart + catFoldSizes(foldIndex) - 1;
        
        % Select all of the vectors in this fold.
        foldVectors = X_sorted(foldStart : foldEnd, :);
        foldCats = y_sorted(foldStart : foldEnd, :);
        
        % If this fold is to be used for validation in this round...
        if (foldIndex == iterationNum) %%%%%%
            % Append the vectors to the validation set.
            X_test = [X_test; foldVectors];
            y_test = [y_test; foldCats];
        % Otherwise, use the fold for training.
        else
            % Append the vectors to the training set.
            X_train = [X_train; foldVectors];
            y_train = [y_train; foldCats];
        end
        
        % Update the starting index of the next fold.
        foldStart = foldEnd + 1;
    end
    
    % Set the starting index of the next category.
    catStart = catStart + vecsPerCat(catIndex);
    
end

end
