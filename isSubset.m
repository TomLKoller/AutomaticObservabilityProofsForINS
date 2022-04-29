function isSub = isSubset(Subset,SuperSet)
isSub=all(ismember(Subset,SuperSet));
end

