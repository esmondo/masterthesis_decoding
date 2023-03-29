if ~exist('topoplot','file')
    addpath /export/data/reichert/toolbox/ft4topoplot/
end

cfgfile = '/export/data/reichert/capLayout/Neuromag_helmet.mat';

cfg = topoprepare(cfgfile);

% type help topoplot to learn about the optional fields in the cfg
% structure


tidx = timeNew>0.175 & timeNew<0.185; % average range around 180 ms
megframe = double (mean(tValues(:,tidx),2)); % convert to double

figure;
toposhow(cfg, megframe);