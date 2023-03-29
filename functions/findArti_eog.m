function eogxcl = findArti_eog(data_eog)

[x1, x2, x3] = size(data_eog);

crit = 




    X                   = reshape( eog.data,[x1*x2 x3]);
    v                   = var( X );
    crit                = mean( v )+std( v )*vco;
    eogxcl              = find( v>crit );

end