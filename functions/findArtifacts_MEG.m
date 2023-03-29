function artFind = findArtifacts_MEG(data)

crit =  3e-12;
artShow = abs(data) > crit;
artSum = sum(artShow);

artFind = permute(artSum,[3 2 1]);
artFind = sum(artFind,2);
artFind = find(artFind ~= 0).';

end