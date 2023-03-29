function idxArtifacts = findArti(data_meg) 

 [x1 x2 x3]          = size( data_meg );
    X                   = reshape( data_meg,[x1*x2 x3]);
    megabsxcl           = find( max( X )> 3e-12 );
    X( :,megabsxcl )    = 0;
    v                   = var( X );
    crit                = mean( v )+std( v )*4;
    
    idxArtifacts              = [find( v>crit ) megabsxcl];

end