function [dataDS,timeDS] = dwsample(idxChan,Fs_new,srate,data,time)

%% Downsample data

[p,q] = rat(Fs_new/srate);
t = ceil(size(data,2) / q);

timeDS = linspace(time(1), time(end), t);

dataDS = zeros(idxChan,ceil(length(timeDS)),size(data,3));

    for ds = 1:size(data,3)
        dataDS(:,:,ds) = resample(data(:,:,ds)',p,q)';
    end

end