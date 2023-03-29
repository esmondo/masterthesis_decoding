function [vecsPerCat] = getVecsPerCat(X, y, categories)


% Get the number of categories present.
numCats = length(categories);

% 'vecsPerCat' will store the number of input vectors belonging to each
% category of label --> 0 for LVF ; 1 for RVF
vecsPerCat = zeros(numCats, 1);

% For each category...
for (i = 1 : numCats)
    
    % Get the ith category; store the category value in column 1.
    cat = categories(i);
    
    % Count the number of input vectors with that category.
    vecsPerCat(i, 1) = sum(y == cat);    
    
end

end