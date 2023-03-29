clear all
close all
clc


% ID and Date
subj = { 'mk48', 'kn19', 'hl64', 'pd46', 'hs85', 'cf15', 'ma72',...
    'xp38', 'kq70', 'ip21', 'lv24', 'ka46', 'xd74', 'sl74', 'os39', 'mf73'};
date = { '190412', '190520', '190521', '190524', '190529', '190531', '190614',...
    '190617', '190710', '190711', '190716', '190718', '190812', '190820', '190822', ...
    '190823'};

L4 = [2, 7, 8, 10, 14, 15 ]; 
L5 = [1, 3, 4, 5, 6, 9, 12, 13, 16]; 
L6 = 11; 

for s = 1:16

if ismember( s,L4 )
    appendix = {'','-1','-2','-3','-4'};   
elseif ismember( s,L5 )
    appendix = {'','-1','-2','-3','-4','-5'}; % fuer den Fall, dass wir mehr Daten haben
elseif ismember( s,L6 )
    appendix = {'','-1','-2','-3','-4','-5', '-6'}; 
end

% 1, 3, 4, 5, 6, 9, 12, 13, 16  appendix 5
% 2, 7, 8, 10, 14, 15  appendix 4
% 11 appendix 6

%% hier laden wir die MEG/EEG Daten
addpath( genpath( '/export/data/reichert/toolbox/fieldtripSVN' )); % hier laden wir die Funktionen aus fieldtrip

meg.data = []; % hier legen wir eine structure an, in die die Daten geschrieben werden
eog.data = []; % hier legen wir eine structure an, in die die Daten geschrieben werden
trigchan = [];

for ap = 1:length( appendix ) % das ist die Schleife ueber die appendizes
    if s<10;
        filename = (['/export/data/wienke/data/rawdata/' date{ s } '/mindwander_n2pc_0', num2str(s), '_', subj{ s } '' appendix{ ap } '.fif']);
    else
        filename = (['/export/data/wienke/data/rawdata/' date{ s } '/mindwander_n2pc_' ,num2str( s ),'_' subj{ s } '' appendix{ ap } '.fif']);
    end;

    if ap == 1;
        meg.header  = ft_read_header( filename );
        meg.srate = meg.header.Fs;
        eog.srate = meg.header.Fs;
        meg.channels = strmatch('MEG',meg.header.label); % hier suchen wir nach channels, die mit MEG anfangen
        meg.channels = meg.channels( 1:3:end ); % beginnend mit dem ersten nehmen wir immer nur jeden dritten, da das die Magnetometer sind
        eog.channels = strmatch( 'EOG',meg.header.label ); % das sind die EOG channels
    end
    X = single( ft_read_data( filename ));
%     X( :,2:2:end,: ) = []; 
    
%     if ap  == 1
%         meg.timepoints = 1:length(X);
%     else
%         meg.timepoints = [meg.timepoints meg.timepoints(end)+1:meg.timepoints(end)+length(X)];
%     end
    
    meg.data = cat( 2,meg.data,X( meg.channels,: ));
    eog.data = cat( 2,eog.data,X( eog.channels,: ));
    trigchan = cat( 2,trigchan,X( end-2,: ));
    
    clear X
end

% Triggercodes.m

% eeg.srate = 1e3; 
% meg.srate = 1e3;

%% Results-File laden
load(['/export/data/wienke/data/rawdata/' subj{s} '_coloredgabors.mat'])

if s==16
    results(501:1200,:)=[];
end
%% hier werden nur die trigger onsets genommen
trigchan( find( diff( trigchan ) < 1 )+1) = 0;

rmpath( genpath( '/export/data/reichert/toolbox/fieldtripSVN' )); % ab hier werden wir nicht mehr auf die fieldtrip skripte zurueckgreifen

%% vorverarbeitung EEG
addpath(genpath('/export/data/duerschm/allscripts'));

% Notch filter details
% nff             = 50:50:300; % notch filter frequencies
% F               = [nff-2; nff+2]';
% 

% band pass filter eog
eog.data            = double(eog.data);
eog                 = ieegBandpassFilter( eog,[1 30]);
eog.data = single( eog.BP - repmat( mean( eog.BP,2 ),[1 length( eog.BP )]));
eog = rmfield( eog, 'BP' );

% meg.data        = double( meg.data );


%     meg             = ieegNotchFilter( meg,F ); %hier sagt er out of memory
% for c = 1:size(meg.data,1 ); % hier filtern wir pro Kanal, da die Daten sonst f�r Matlab zu gro� sind
%     disp( c );
%     meg2.data = double( meg.data( c,: ));
%     meg2.srate = meg.srate;
%     meg2 = ieegNotchFilter( meg2,F );
%     meg.data( c,: ) = meg2.data;
% end;

%% Bandpass Filter
for c = 1:size( meg.data,1 );
    disp( c )
    meg2.data = double( meg.data( c,: ));
    meg2.srate = meg.srate;
    meg2 = ieegBandpassFilter( meg2,[.1 30]);
    meg.data( c,: ) = single( meg2.BP - repmat( mean( meg2.BP,2 ),[1 length( meg2.BP )]));
end
clear meg2
%% Bandpass Filter Hilbert
% for c = 1:102;
%     disp( c )
%     meg2.data = meg.data( c,: );
%     meg2.srate = meg.srate;
%     meg2 = ieegBandpassFilterAmplitude( meg2,[1 200]);
%     meg.data( c,: ) = single( meg2.BP - repmat( mean( meg2.BP,2 ),[1 length( meg2.BP )]));
% end

%% settings fuer epochierung
trigs_used          = [0 10 20 30 40 50 60];
trigchan(~ismember(trigchan, trigs_used))=0;
trialstart          = find( trigchan == 20 );
trialdur            = round( -1*meg.srate:2*meg.srate );
trials              = bsxfun( @plus,trialstart',trialdur );
[x2 x3]             = size( trials );

% meg.timepoints      = meg.timepoints/meg.srate;
% meg.timepoints      = meg.timepoints(trialstart);


x1                  = size( eog.data,1 );
eog.data            = ( reshape( eog.data( :,trials' ),[x1 x3 x2]));
eog.time            = trialdur/eog.srate;

x1                  = size( meg.data,1 );
meg.data            = ( reshape( meg.data( :,trials' ),[x1 x3 x2]));
meg.time            = trialdur/meg.srate;
%     meg                 = rmfield(meg,'BP' )

% baseline correction
bsl                 = find( meg.time>-.5&meg.time<0 );
meg.data            = single( meg.data - repmat( mean( meg.data( :,bsl,: ),2 ),[1 length( meg.data )]));

%% �berz�hlige Trials entfernen
if s == 11
    meg.data(:,:,101:143)       = [];
    eog.data(:,:,101:143)       = [];
%     meg.timepoints(101:143)     = [];
end
if s == 12
    meg.data(:,:,1:34)          = [];
    eog.data(:,:,1:34)          = [];
%     meg.timepoints(1:34)        = [];
end
if s == 13
    meg.data(:,:,1:11)          = [];
    eog.data(:,:,1:11)          = [];
%     meg.timepoints(1:11)        = [];
end

if s == 16
    meg.data(:,:,1:700)          = [];
    eog.data(:,:,1:700)          = [];
%     meg.timepoints(1:700)        = [];
end
%% Results-File zuf�gen
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

% filename = '/export/data/wienke/data/rawdata/euclidDistance.mat';
% load( filename )
% meg.distance=[];
% for k = 1:length(meg.position)
% meg.distance(k) = d2c(meg.position(k));
% end

% x=[7, 10];
% meg.distance2( ismember( meg.position, x )) = 1;
% va( 1 ) = mean( d2c( x ));
% 
% x=[4, 8, 11, 13];
% meg.distance2( ismember( meg.position, x )) = 2;
% va( 2 ) = mean( d2c( x ));
% 
% x=[5, 14];
% meg.distance2( ismember( meg.position, x )) = 3;
% va( 3 ) = mean( d2c( x ));
% 
% x=[1, 9, 12, 16];
% meg.distance2( ismember( meg.position, x )) = 4;
% va( 4 ) = mean( d2c( x ));
% 
% x=[2, 6, 15, 17];
% meg.distance2( ismember( meg.position, x )) = 5;
% va( 5 ) = mean( d2c( x ));
% 
% x=[3, 18];
% meg.distance2( ismember( meg.position, x )) = 6;
% va( 6 ) = mean( d2c( x ));
% 
% for k = 1:6
%     meg.distance2( meg.distance2 == k ) = va( k );
% end

%% Artifact Rejection
%% reject eye movements
addpath(genpath('/export/data/duerschm/allscripts'))
% meg.cleandata = rejectComponent( meg.data,squeeze( eog.data( 1,:,: )));
% meg.cleandata = rejectComponent( meg.cleandata,squeeze( eog.data( 2,:,: )));

meg.data = single( rejectComponent( meg.data,squeeze( eog.data( 1,:,: ))));
meg.data = single( rejectComponent( meg.data,squeeze( eog.data( 2,:,: ))));
rmpath(genpath('/export/data/duerschm/allscripts'))

%% Common average referencing
% to subtract the average over all electrodes from each electrodes for each time 
% point. This distributes the "responsibility" over all electrodes, rather than 
% assigning it to only one of them
% http://imaging.mrc-cbu.cam.ac.uk/meg/IntroEEGMEG

% mean(meg.data gibt mittlere Amplitude aller Elektroden f�r jeden
% Zeitpunkt und trial; also 1xNtimepointsxNtrials
% repmat dann um das ganze f�r 102 Kan�le zu erweitern und 1 weil schon NtimexNtrial 
meg.data = meg.data-repmat( mean( meg.data),[size( meg.data,1 ) 1] );



%% Doppelte Tastendr�cke
% excl = find( meg.response == 999);
% 
% meg.data( :,:,excl )      = [];
% eog.data( :,:,excl )      = [];
% meg.position( :,excl )    = [];
% meg.RT( :,excl )          = [];
% meg.fokus( :,excl )       = [];
% meg.winkel( :,excl )      = [];
% meg.distance( :,excl )    = [];
% meg.distance2( :,excl )   = [];
% meg.response( :,excl )    = [];
% meg.side( :,excl )        = [];
% meg.timepoints( :,excl )  = [];
% 
% clear excl 

%% Absoluter Wert 3e-12
% [x1, x2, x3] = size(meg.data);
% crit = 3e-12;
% 
% for ex = 1:x3
%     tmp = find( meg.data( :,:,ex) > crit);
%     if isempty( tmp ) == 0
%         excl( ex) = 1;
%     else
%         excl( ex ) = 0;
%     end
% end
% excl = logical( excl );
% 
% meg.data( :,:,excl )      = [];
% eog.data( :,:,excl )      = [];
% meg.position( :,excl )    = [];
% meg.RT( :,excl )          = [];
% meg.fokus( :,excl )       = [];
% meg.winkel( :,excl )      = [];
% % meg.distance( :,excl )    = [];
% % meg.distance2( :,excl )   = [];
% meg.response( :,excl )    = [];
% meg.side( :,excl )        = [];
% % meg.timepoints( :,excl )  = [];
% 
% clear excl tmp ex crit

%% Blinks
[x1, x2, x3] = size(eog.data);

for tr = 1:x3
EM_ch1(tr) = max(abs(eog.data(1,find( meg.time > -.5 & meg.time < .5 ),tr)))';
EM_ch2(tr) = max(abs(eog.data(2,find( meg.time > -.5 & meg.time < .5 ),tr)))'; % channel 2 vetikale Augenbewegungen?
end

% excl=find(EM_ch2 > 1e-4);
excl=find(EM_ch2 > mean(EM_ch2)+std(EM_ch2));
% excl=find(EM_ch2 > mean(EM_ch2)+std(EM_ch2)*2);

meg.data( :,:,excl )      = [];
eog.data( :,:,excl )      = [];
meg.position( :,excl )    = [];
meg.RT( :,excl )          = [];
meg.fokus( :,excl )       = [];
meg.winkel( :,excl )      = [];
% meg.distance( :,excl )    = [];
% meg.distance2( :,excl )   = [];
meg.response( :,excl )    = [];
meg.side( :,excl )        = [];
% meg.timepoints( :,excl )  = [];

clear excl tr EM_ch1 EM_ch2


% Varianz der MEG-Daten
% [x1, x2, x3]  = size( meg.data );

% X           = reshape( meg.data,[x1*x2 x3]);
% v           = var( X );
% crit        = mean( v )+std( v )*5;
% crit        = 1.5e-24;
% crit        = 3e-24;
% excl        = find( v>crit );

% meg.data( :,:,excl )      = [];
% eog.data( :,:,excl )      = [];
% meg.position( :,excl )    = [];
% meg.RT( :,excl )          = [];
% meg.fokus( :,excl )       = [];
% meg.winkel( :,excl )      = [];
% meg.distance( :,excl )    = [];
% meg.distance2( :,excl )   = [];
% meg.response( :,excl )    = [];
% meg.side( :,excl )        = [];
% meg.timepoints( :,excl )  = [];

% clear excl crit X v 

%% Speichern
if s<10;
    filename = (['/export/data/wienke/data/Esmondo/MW0' num2str( s ) 'meg.mat']);
else
    filename = (['/export/data/wienke/data/Esmondo/MW' num2str( s ) 'meg.mat']);
end
save( filename,'-v7.3','meg','eog' );

clearvars -except subj date s L4 L5 L6

end
% 
% X = mean(meg.data,3)*1e+14;
% Z = X +repmat([1:5:102*5]',[1 6001]);
% plot( meg.time, Z )