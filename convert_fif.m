clear all
close all
clc



subj  = {'pi51',... 
         'ew16',...
         'pd46',... % kein logfile vorhanden
         'bf98',...
         'wl41',...
         'qp14',...
         'hs85',...
         'mf73',...
         'mu66',...
         'mk48',... 
         'bh27',...
         'ip21',...
         'gj66',... % cans use gj666, 1 logfile for 2 sessions (wrong naming)
         'pi51(2)',... % cnat use this, second measurment of same particpants
         'gj66(2)',... % cant use gj66, 1 logfile for 22 sessions
         'pz96',...
         'ap38',...
         'sc17'...
         };





date = { '200707', '200714', '200717', '200729', '200731', '200806', '200813',...
    '200814', '200818', '200908', '200910', '201001', '201002', '201006','201015','201019', '201028', '201028'};

L5 = [2 3 4 5 6 8 9 11 12  16 17 18];
L6 = [1 7 10];

subs = [1:2, 4:12, 16:18];

sidx = 1;


for s = subs;
    
%     s = subs(i);
    
    addpath( genpath( '/export/data/reichert/toolbox/fieldtripSVN' )); % hier laden wir die Funktionen aus fieldtrip
    
    if s < 10
    %if subs( s ) < 10
        filename = (['/export/data/database/MEG/mindwandering_sl_sss/mindwander_sl_0' num2str(s) '_' subj{s} '-rall_mci_sss_ds4.fif']);
    elseif s >=10
        filename = (['/export/data/database/MEG/mindwandering_sl_sss/mindwander_sl_' num2str(s) '_' subj{s} '-rall_mci_sss_ds4.fif']);
    end
    
    meg.header          = ft_read_header( filename );
    meg.srate           = meg.header.Fs;
    eog.srate           = meg.header.Fs;
    meg.channels        = strmatch('MEG',meg.header.label); % hier suchen wir nach channels, die mit MEG anfangen
    meg.channels        = meg.channels( 1:3:end ); % beginnend mit dem ersten nehmen wir immer nur jeden dritten, da das die Magnetometer sind
    eog.channels        = strmatch( 'EOG',meg.header.label ); % das sind die EOG channels
    eeg.channels        = strmatch( 'EEG',meg.header.label ); % das sind die EEG channels
    
    X = ft_read_data( filename );
    rmpath( genpath( '/export/data/reichert/toolbox/fieldtripSVN' )); % hier laden wir die Funktionen aus fieldtrip
    
    eeg.srate = meg.srate;
    
    %% sorting the data
    meg.data            = X( meg.channels,: );
    eog.data            = X( eog.channels,: );
    eeg.data            = X( eeg.channels,: );
    trigchan            = X( end-2,: );
    
    %% hier werden nur die trigger onsets genommen
    trigchan(find(diff(trigchan) < 1) +1) = 0;
    trigchan(find(trigchan < 0)) = [];
    
    
    %% Results-File laden
    load(['/export/data/Vogelgesang/betamotor/data/chris3/logfiles/' subj{ s } '_coloredgaborsEEG.mat'])
    
    
    rmpath( genpath( '/export/data/reichert/toolbox/fieldtripSVN' )); % ab hier werden wir nicht mehr auf die fieldtrip skripte zurueckgreifen
    
    %% vorverarbeitung EEG
    addpath(genpath('/export/data/duerschm/allscripts'));
    
    % Notch filter details
    nff             = 50:50:200; % notch filter frequencies
    F               = [nff-2; nff+2]';
    
    
    % band pass filter eog
    eog.data            = double(eog.data);
    eog                 = ieegBandpassFilter( eog,[1 30]);
    eog.data = single( eog.BP - repmat( mean( eog.BP,2 ),[1 length( eog.BP )]));
    
    
    meg.data        = double( meg.data );
    
    %% rereferencing of EEG
    eeg.data = eeg.data( 1:29,:,: )-repmat( eeg.data( 30,:,: )./2,[29 1 1] );
    
    %% Notch filter
    %     meg             = ieegNotchFilter( meg,F ); %hier sagt er out of memory
    fprintf('notch filtering \n' )
    for c = 1:102; % hier filtern wir pro Kanal, da die Daten sonst f�r Matlab zu gro� sind
        meg2.data = meg.data( c,: );
        meg2.srate = meg.srate;
        meg2 = ieegNotchFilter( meg2,F );
        meg.data( c,: ) = meg2.data;
    end;
    
    %% Bandpass Filter
    for c = 1:102;
        disp( c )
        meg2.data       = meg.data( c,: );
        meg2.srate      = meg.srate;
        meg2            = ieegBandpassFilter( meg2,[1 200]);
        meg.data( c,: ) = single( meg2.BP - repmat( mean( meg2.BP,2 ),[1 length( meg2.BP )]));
    end
    
        %% notch filter eeg
        for c = 1:29; % hier filtern wir pro Kanal, da die Daten sonst f�r Matlab zu gro� sind
            eeg2.data       = eeg.data( c,: );
            eeg2.srate      = meg.srate;
            eeg2            = ieegNotchFilter( eeg2,F );
            eeg.data( c,: ) = eeg2.data;
        end;
    
        %% Bandpass Filter
        for c = 1:29;
            disp( c )
            eeg2.data       = eeg.data( c,: );
            eeg2.srate      = meg.srate;
            eeg2            = ieegBandpassFilter( eeg2,[1 200]);
            eeg.data( c,: ) = single( eeg2.BP - repmat( mean( eeg2.BP,2 ),[1 length( eeg2.BP )]));
        end
    
    %% settings fuer epochierung
    trigs_used          = [0 10 20 30 40 50 60];
    trigchan(~ismember(trigchan, trigs_used))=0;
    trialstart          = find( trigchan == 20 );
    trialdur            = round( -1*meg.srate:3*meg.srate );
    trials              = bsxfun( @plus,trialstart',trialdur );
    [x2 x3]             = size( trials );
    
    
    %% zuerst schauen wir nach onsets des visual displays (vpo) und ob tatsächlich vor jedem VD auch ein Fixationskreuz präsentiert wurde
    vpo             = find( trigchan== 20 );
    idx = 1;
    for o = 1:length( vpo );
        waitingasterisks( o,length( vpo ));
        fco                 = find( trigchan( 1:vpo( o )) == 10,1,'last' );
        vdlatency           = vpo( o ) - fco;
        if vdlatency < meg.srate*2;
            trialnr( idx )  = o;
            t2vd( idx )     = vdlatency; % time to visual display onset
            idx = idx+1;
        end
    end
    
    %% jetzt müssen wir fehlende trials aus der Matrix nehmen, in der wir Verhaltensergebnisse haben
    twof                = setdiff( 1:length( results ),trialnr ); % trials without fixation
    results( twof,: )   = [];
    
    
    %     meg.timepoints      = meg.timepoints/meg.srate;
    %     meg.timepoints      = meg.timepoints(trialstart);
    
    
    x1                  = size( eog.data,1 );
    eog.data            = ( reshape( eog.data( :,trials' ),[x1 x3 x2]));
    eog.time            = trialdur/eog.srate;
    
    x1                  = size( meg.data,1 );
    meg.data            = ( reshape( meg.data( :,trials' ),[x1 x3 x2]));
    meg.time            = trialdur/meg.srate;
    %     meg                 = rmfield(meg,'BP' )
    
    
    meg.t2vd            = t2vd;
    meg.trialnr         = trialnr;
    
        x1                  = size( eeg.data,1 );
        eeg.data            = ( reshape( eeg.data( :,trials' ),[x1 x3 x2]));
        eeg.time            = trialdur/eeg.srate;
    
    
    % baseline correction
    %     bsl                 = find( meg.time>-.5&meg.time<0 );
    %     meg.data            = single( meg.data - repmat( mean( meg.data( :,bsl,: ),2 ),[1 length( meg.data )]));
    %     eeg.data            = single( eeg.data - repmat( mean( eeg.data( :,bsl,: ),2 ),[1 length( eeg.data )]));
    %
    
    %% variance
    vco                 = 5;
    
    [x1 x2 x3]          = size( meg.data );
    X                   = reshape( meg.data,[x1*x2 x3]);
    megabsxcl           = find( max( X )>.5e-11 );
    X( :,megabsxcl )    = 0;
    v                   = var( X );
    crit                = mean( v )+std( v )*vco;
    %     crit                = 1.5e-24;
    megxcl              = [find( v>crit ) megabsxcl];
    
    
    [x1 x2 x3]          = size( eog.data );
    X                   = reshape( eog.data,[x1*x2 x3]);
    v                   = var( X );
    crit                = mean( v )+std( v )*vco;
    eogxcl              = find( v>crit );
    
    
        [x1 x2 x3]          = size( eeg.data );
        X                   = reshape( eeg.data,[x1*x2 x3]);
    
        eegabsxcl           = find( max( X )>2e-4 );
        % if for some reason the general threshold is too low we find an
        % individual threshold
        if length( eegabsxcl )>200;
            x                   = max( X );
            p                   = cdf( 'Normal',x,mean( x ),std( x ));
            crit                = min( x( find( p>.95 )));
            eeg.cutoff          = crit;
            eegabsxcl           = find( x>crit );
        end
        
        X( :,eegabsxcl )    = 0;
        v                   = var( X );
        crit                = mean( v )+std( v )*vco;
        eegxcl              = [find( v>crit ) eegabsxcl];
    
    
%     xcl                 = unique( [eogxcl megxcl] );
        xcl                 = unique( [eogxcl megxcl eegxcl] );
    
    
    
    %% trial exclusion
    eeg.data( :,:,xcl ) = [];
    meg.data( :,:,xcl ) = [];
    eog.data( :,:,xcl ) = [];
    results( xcl,: )    = [];
    meg.t2vd( xcl )     = [];
    meg.trialnr( xcl )  = [];
    
    eeg.trialsxcl       = xcl;
    meg.trialsxcl       = xcl;
    eog.trialsxcl       = xcl;
    
    
    
    %% reject eye movements
    addpath(genpath('/export/data/duerschm/allscripts'))
    
    %     eeg.data = single( rejectComponent( eeg.data,squeeze( eog.data( 1,:,: ))));
    %     eeg.data = single( rejectComponent( eeg.data,squeeze( eog.data( 2,:,: ))));
    %     eeg.data = eeg.data - repmat( mean( eeg.data ),[29 1 1]);
    %
    %     meg.data = single( rejectComponent( meg.data,squeeze( eog.data( 1,:,: ))));
    %     meg.data = single( rejectComponent( meg.data,squeeze( eog.data( 2,:,: ))));
    %     meg.data = meg.data - repmat( mean( meg.data ),[102 1 1]);
    
    %% downsampling
    %     eeg.data        = eeg.data( :,1:2:end,: );
    %     eeg.time        = eeg.time( 1:2:end );
    %     eeg.srate       = 500;
    %
    %     eog.data        = eog.data( :,1:2:end,: );
    %     eog.time        = eog.time( 1:2:end );
    %     eog.srate       = 500;
    %
    %     meg.data        = meg.data( :,1:2:end,: );
    %     meg.time        = meg.time( 1:2:end );
    %     meg.srate       = 500;
    
    
    
    
    %% Speichern
    if s < 10;
        filename = (['/export/data/Vogelgesang/CFC/Data/preprocessed_data/MW0' num2str( s ) 'meg.mat']);
    else
        filename = (['/export/data/Vogelgesang/CFC/Data/preprocessed_data/MW' num2str( s ) 'meg.mat']);
    end
    save( filename,'-v7.3','meg','eeg','eog', 'results');
    
    %     if s<10;
    %         filename = (['/export/data/Vogelgesang/betamotor/data/chris3/MW0' num2str( s ) 'eeg.mat']);
    %     else
    %         filename = (['/export/data/Vogelgesang/betamotor/data/chris3/MW' num2str( s ) 'eeg.mat']);
    %     end
    %     save( filename,'-v7.3','eeg','eog','results');
    
    %     clearvars -except subj s L5 L6 date subs
    
end