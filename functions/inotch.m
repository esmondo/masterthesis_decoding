function ieeg = inotch( varargin )

ieeg = varargin{ 1 };
freq = varargin{ 2 };

[szx,szy,szz] = size( ieeg.data );
nbd = ieeg.data; %notch band data
nbd                       = reshape( permute( nbd,[1 3 2]),[szz*szx szy]);
% fprintf('asterisks denote frequency\n')
for f = 1:size( freq,1 )
    [a,b]                   = butter( 2,[freq( f,1 )/( ieeg.srate/2 ), freq( f,2 )/( ieeg.srate/2 )],'stop');
    nbd                     = permute( filtfilt( a,b,permute( nbd,[2 1]) ),[2 1]);
end
ieeg.data      = permute( reshape( nbd,[szx szz szy]),[1 3 2]);