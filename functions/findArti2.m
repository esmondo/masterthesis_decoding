function megxcl = findArti2(data_meg)


[x1, x2, x3] = size(data_meg);
crit = 3e-12;

for ex = 1:x3
    tmp = find( data_meg( :,:,ex) > crit);
    if isempty( tmp ) == 0
        megxcl( ex) = 1;
    else
        megxcl( ex ) = 0;
    end
end

xcl = logical( megxcl );

megxcl = find(xcl == 1); % use the logical array to reject the unwanted trials / extract only clean trials

end