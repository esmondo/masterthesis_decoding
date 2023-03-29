function predictions = testClassifier(Mdl,X_test)

X_test = reshape(permute(X_test,[3,2,1]),...
         [size(X_test,3),size(X_test,1)*size(X_test,2)]);
     
predictions = predict(Mdl,X_test);

end