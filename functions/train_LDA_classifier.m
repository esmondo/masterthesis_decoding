function [Mdl,X_train] = train_LDA_classifier(X_train,y_train)


%% Features

% X_train = permute(X_train,[3 2 1]);
X_train = reshape(permute(X_train,[3,2,1]),...
         [size(X_train,3),size(X_train,1)*size(X_train,2)]);
y_train = y_train';


%% Train LDA Classifier

Mdl = fitcdiscr(X_train,y_train);

end