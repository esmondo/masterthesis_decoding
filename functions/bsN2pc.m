function bcdat = bsN2pc(data,time)

% indexing prestimulus onset from -200 ms
prestimVolt = find(time>=-0.2 & time<=0);
% prestimVal = time(prestimVolt);

% mean of prestimulus voltage
prestimMean = mean(data(prestimVolt),2);

% baseline correction
bcdat = bsxfun(@minus, data, prestimMean);


end