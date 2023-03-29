function [Mdl,X_train] = trainClassifier(X_train,y_train)


%% Features

% X_train = permute(X_train,[3 2 1]);
X_train = reshape(permute(X_train,[3,2,1]),...
         [size(X_train,3),size(X_train,1)*size(X_train,2)]);
y_train = y_train';

%% C parameter
if sum(y_train)/length(y_train) == 0
    C = 1;
elseif sum(y_train)/length(y_train) == 1
    C = 1;
else
    C = 1/mean(diag(X_train*X_train'));
end

%% Train SVM Classifier
Mdl = fitcsvm(X_train,y_train,...
    'KernelFunction', 'linear',...
    'KernelScale', 1,...
    'BoxConstraint',C);

end