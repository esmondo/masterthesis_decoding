%% Christoph adapted this script initially created by christian

addpath(genpath('/export/data/esmondo'))
% datapath = '/data/reichert/MEG/4Esmondo/data/rawdata/';
% logpath = '/data/reichert/MEG/4Esmondo/data/logfiles/';
datapath = '/export/data/esmondo/';
logpath = '/export/data/esmondo/';


if ~exist('ft_defaults','file')
    addpath /data/reichert/SVN/fieldtrip.git/trunk
    ft_defaults
end

% ID and Date
subj = { 'mk48', 'kn19', 'hl64', 'pd46', 'hs85', 'cf15', 'ma72',...
    'xp38', 'kq70', 'ip21', 'lv24', 'ka46', 'xd74', 'sl74', 'os39', 'mf73'};


% for s = 3; %1:length(subj)

s=1;

%% hier laden wir die MEG/EEG Daten
meg = struct;

flist = dir(sprintf([datapath '*n2pc_%02d*.fif'],s));
filename = [datapath flist.name];

meg.header  = ft_read_header( filename );
meg.srate = meg.header.Fs;
eog.srate = meg.header.Fs;
grad.srate = meg.header.Fs;
chanidx = strmatch('MEG',meg.header.label); % hier suchen wir nach channels, die mit MEG anfangen
meg.channels = chanidx( 1:3:end ); % beginnend mit dem ersten nehmen wir immer nur jeden dritten, da das die Magnetometer sind
grad.channels = setdiff(chanidx, meg.channels);
eog.channels = strmatch( 'EOG',meg.header.label ); % das sind die EOG channels

X = single( ft_read_data( filename ));

meg.data = X( meg.channels,: );
grad.data =  X( grad.channels,: );
eog.data = X( eog.channels,: );
trigchan = X( end-2,: );

clear X


% Triggercodes.m

%% Results-File laden
flist = dir([logpath subj{s} '_coloredgabors.mat']);
load([logpath flist.name])

if s==16
    results(501:1200,:)=[];
end
%% hier werden nur die trigger onsets genommen
trigchan( find( diff( trigchan ) < 1 )+1) = 0;

%% settings fuer epochierung
trigs_used          = [0 10 20 30 40 50 60];
trigchan(~ismember(trigchan, trigs_used))=0;
trialstart          = find( trigchan == 20 );

fprintf('subj %i: %i trials found\n',s, length(trialstart));

trialdur            = round( -1*meg.srate:2*meg.srate );
trials              = bsxfun( @plus,trialstart',trialdur );
[x2 x3]             = size(trials);


x1                  = size( eog.data,1 );
eog.data            = ( reshape( eog.data( :,trials' ),[x1 x3 x2]));
eog.time            = trialdur/eog.srate;

x1                  = size( meg.data,1 );
meg.data            = ( reshape( meg.data( :,trials' ),[x1 x3 x2]));
meg.time            = trialdur/meg.srate;

x1                  = size( grad.data,1 );
grad.data            = ( reshape( grad.data( :,trials' ),[x1 x3 x2]));
grad.time            = trialdur/grad.srate;

%% Remove excessive trials
if s == 11
    meg.data(:,:,101:143)       = [];
    eog.data(:,:,101:143)       = [];
    grad.data(:,:,101:143)       = [];
%     meg.timepoints(101:143)     = [];
end
if s == 12
    meg.data(:,:,1:34)          = [];
    eog.data(:,:,1:34)          = [];
    grad.data(:,:,1:34)          = [];
%     meg.timepoints(1:34)        = [];
end
if s == 13
    meg.data(:,:,1:11)          = [];
    eog.data(:,:,1:11)          = [];
    grad.data(:,:,1:11)          = [];
%     meg.timepoints(1:11)        = [];
end

if s == 16
    meg.data(:,:,1:700)          = [];
    eog.data(:,:,1:700)          = [];
    grad.data(:,:,1:700)          = [];
%     meg.timepoints(1:700)        = [];
end
%% Results-File zuf??gen
keycode_side = [55 56];
results( find( results( :,1 ) == keycode_side( 1 )),1 ) = 0;        % links geantwortet
results( find( results( :,1 ) == keycode_side( 2 )),1 ) = 1;        % 2  % rechts geantwortet

keycode_fokus =[49 50 51 52 53];
results( find( results( :,5 ) == keycode_fokus( 1 )),5 ) = 1;     % completely off-task
results( find( results( :,5 ) == keycode_fokus( 2 )),5 ) = 2;     
results( find( results( :,5 ) == keycode_fokus( 3 )),5 ) = 3;     
results( find( results( :,5 ) == keycode_fokus( 4 )),5 ) = 4;    
results( find( results( :,5 ) == keycode_fokus( 5 )),5 ) = 5;     % completely on-task

side=[];
side( find( results( :,3 )<0 )) = 0;    % tilt nach links
side( find( results( :,3 )>0 )) = 1; %2 % tilt nach rechts
results = [results side'];

meg.response        = results(:,1)';
meg.RT              = results(:,2)';
meg.winkel          = abs(results(:,3))';
meg.position        = results(:,4)';
meg.fokus        	= results(:,5)';    
meg.side            = results(:,7)';


%% Speichern
filename = sprintf('%sMW%02dmegSSS.mat',datapath,s);
filenamegrad = sprintf('%sMW%02dgradSSS.mat',datapath,s);

save( filename,'-v7.3','meg','eog' );
save( filenamegrad,'-v7.3','grad' );

% end


