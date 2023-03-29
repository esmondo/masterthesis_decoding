function artFind = findArtifacts(data)

artShow = abs(data) > 3e-12;
artSum = sum(artShow);

artFind = permute(artSum,[3 2 1]);
artFind = sum(artFind,2);
artFind = find(artFind ~= 0).';

end