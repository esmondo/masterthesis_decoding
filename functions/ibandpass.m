function ieeg = ibandpass( varargin ) % varargin> variable Anzahl an input argumenten

ieeg = varargin{ 1 };
freq = varargin{ 2 };

[szx,szy,szz] = size( ieeg.data );
ieeg.BP = zeros( szx,szy,szz,size( freq,1 ));
% fprintf('asterisks denote frequency\n')
for f = 1:size( freq,1 )
%     if ismember( f/50,1:10000 ) == 1;
%         fprintf( '*\n' )
%     else
%         fprintf( '*' )
%     end
    
    X                       = reshape( permute( ieeg.data,[1 3 2]),[szz*szx szy]);
    [a,b]                   = butter( 2,[freq( f,1 )/( ieeg.srate/2 ), freq( f,2 )/( ieeg.srate/2 )]);
    ffX                     = permute( filtfilt( a,b,permute( X,[2 1]) ),[2 1]);
%         ffX                     = permute( filter( a,b,permute( X,[2 1]) ),[2 1]);
    ieeg.BP( :,:,:,f )      = permute( reshape( ffX,[szx szz szy]),[1 3 2]);
    
end