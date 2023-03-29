function accuracyM = acc(Y,predictions,nreps)

accuracy = sum(predictions==repmat(Y',1,nreps))./length(Y)*100;
accuracyM = mean(accuracy,2);

end