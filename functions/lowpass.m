function signal = lowpass(cutoff,srate,signal)


Fn = srate/2;      % nyquist frequency

%%% butterworth filter design
[z,p,k] = butter(4, cutoff/Fn ,'low');
% Transfer function properties:
% z = zeros
% p = poles
% k = gain

%%% convert zero-pole-gain filter to second order sections form
[sos,g] = zp2sos(z,p,k);
% sos = second-order sections
% g = gain equivalent to k

%% Filtering 
for tr = 1:size(signal,3)

    signal(:,:,tr) = filtfilt(sos, g, signal(:,:,tr)')';

end


end