function ds = fif2mat(dataFile,varargin)
% ds = neuromag2mat(dataFile,...)
% reads Elekta  NeuroMag data from dataFile and returns a datastructure ds
% optional arguments: 'eegchans' - eeg channel numbers in scan setup
%                     'eogchans' - eog channel numbers in scan setup
%                     'emgchans' - emg channel numbers in scan setup
%                                  if no channels defined, use empty matrix
%                     'bidx' - begin with sampling point bidx
%                     'eidx' - end at sampling point eidx
% requires fieldtrip fileio

% set default values
eegChans = []; % [1:29]
eogChans = [1:2]; % [30,31]
emgChans = []; % [63:64]
bIdx=1;

% check input arguments
for k=1:2:length(varargin),
    if strcmpi(varargin{k},'eegchans'),
        eegChans = varargin{k+1};
    elseif strcmpi(varargin{k},'eogchans'),
        eogChans = varargin{k+1};
    elseif strcmpi(varargin{k},'emgchans'),
        emgChans = varargin{k+1};
    elseif strcmpi(varargin{k},'bidx'),
        bIdx = varargin{k+1};
    elseif strcmpi(varargin{k},'eidx'),
        eIdx = varargin{k+1};
    end
end

% read the data
hdr=ft_read_header(dataFile);
rawdata=ft_read_data(dataFile);

if ~exist('eIdx','var'),
    eIdx=size(rawdata,2);
end

sensitivity=1;  % assume T as default unit
srate=hdr.Fs;
megIdx=zeros(306,1);
eegIdx=zeros(size(eegChans));
eogIdx=zeros(size(eogChans));
emgIdx=zeros(size(emgChans));
mcount = 0;
ecount = 0;

% sort the data
for k=1:length(hdr.label)
    if strcmp(hdr.label{k},'STI101')
        trigIdx=k;
    elseif strcmp(hdr.label{k},'STI016') % photo diode plugged in bit 16
        auxIdx=k;
    elseif strcmp(hdr.label{k},'STI102')
        respIdx=k;
    elseif ~isempty(strfind(hdr.label{k},'MEG')==1)
        mcount = mcount + 1;
        megIdx(mcount)=k;
    elseif ~isempty(strfind(hdr.label{k},'EEG'))
        ecount = ecount + 1;
        eegIdx(ecount)=k;
    elseif ~isempty(strfind(hdr.label{k},'BIO')) || ~isempty(strfind(hdr.label{k},'EOG'))
        tmp  = str2double(hdr.label{k}(4:end));
        if any(tmp==eogChans)          
            eogIdx(eogChans==tmp)=k; 
        elseif any(tmp==emgChans)
            emgIdx(emgChans==tmp)=k; 
        end
    end
end

% put the data in a structure
ds=struct;
ds.srate = srate;
ds.sensitivity = sensitivity;
ds.time = linspace(0,(eIdx-bIdx+1)/srate,eIdx-bIdx+1);
ds.meg = single(rawdata(megIdx,bIdx:eIdx,:));
%channels
ds.eeg = single(rawdata(eegIdx,bIdx:eIdx,:));
ds.eog = single(rawdata(eogIdx,bIdx:eIdx,:));
ds.emg = single(rawdata(emgIdx,bIdx:eIdx,:));
ds.response = int16(rawdata(respIdx,bIdx:eIdx,:));
ds.trigger = int16(rawdata(trigIdx,bIdx:eIdx,:));
% ds.aux = int16(rawdata(auxIdx,bIdx:eIdx,:));
ds.megLabel=cell(length(megIdx),1);[ds.megLabel{:}]=deal(hdr.label{megIdx});
ds.eegLabel=cell(length(eegIdx),1);[ds.eegLabel{:}]=deal(hdr.label{eegIdx});
ds.eogLabel=cell(length(eogIdx),1);[ds.eogLabel{:}]=deal(hdr.label{eogIdx});
ds.emgLabel=cell(length(emgIdx),1);[ds.emgLabel{:}]=deal(hdr.label{emgIdx});
