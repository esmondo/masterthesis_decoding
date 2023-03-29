function [maxPleftChan,maxNleftChan,maxPrightChan,maxNrightChan] = ...
    maxPmaxN(tValues,chansLeft,chansRight,timeN2pcIdx)

%% Left Hemisphere %

% occipito-temporal sensors LVF
tVal_left = tValues(chansLeft,timeN2pcIdx); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxPleft = tVal_left;
maxPleft(maxPleft<=0) = NaN; % highlighting only the positive flux (efflux)
maxPleft = max(maxPleft,[],2);
[maxPleftVal,maxPleftChan] = max(maxPleft(:)); 

maxNleft = tVal_left;
maxNleft(maxNleft>=0) = NaN; % highlighting only the negative flux (influx)
maxNleft = min(maxNleft,[],2);
[maxNleftVal,maxNleftChan] = min(maxNleft(:));


%% Right Hemisphere %

% occipito-temporal sensors RVF
tVal_right = tValues(chansRight,timeN2pcIdx); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxPright = tVal_right;
maxPright(maxPright<=0) = NaN;
maxPright = max(maxPright,[],2);
[maxPrightVal,maxPrightChan] = max(maxPright(:)); 


maxNright = tVal_right;
maxNright(maxNright>=0) = NaN;
maxNright = min(maxNright,[],2);
[maxNrightVal,maxNrightChan] = min(maxNright(:)); 






% gutSensorsLeft = [maxPleft; maxNleft];
% gutSensorsRight = [maxPright; maxNright];
% gutSensors = [maxPleft; maxNleft; maxPright; maxNright];

end