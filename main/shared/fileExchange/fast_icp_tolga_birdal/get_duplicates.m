
% Returns the duplicates of a vector
% 
% Author: Tolga Birdal

function duplicates = get_duplicates(X)
%[uniqueX i j] = unique(X,'first');
% duplicates = 1:length(X);
% duplicates(i) = [];
%duplicates = find(not(ismember(1:numel(X),i)));
uniqueX = unique(X);
countOfX = hist(X,uniqueX);
index = (countOfX~=1 & countOfX~=0);
duplicates = uniqueX(index);
end