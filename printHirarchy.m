function  printHirarchy(filtered_obs_from, check_names)
greater=containers.Map;
for outer_index=1:size(filtered_obs_from,1)
    for index=1:numel(filtered_obs_from{outer_index})
        if ~isKey(greater,getKey(filtered_obs_from{outer_index}{index}))
            greater(getKey(filtered_obs_from{outer_index}{index}))=[outer_index];
        else
        greater(getKey(filtered_obs_from{outer_index}{index}))=[
        greater(getKey(filtered_obs_from{outer_index}{index})), outer_index];
        end
    end
end

skips=zeros(numel(greater.keys),1);

    keys=greater.keys();
for index=1:numel(greater.keys)
    if skips(index)
        continue
    end
    skips(index)=1;
    key=(keys{index});
   set=greater(key);
   if (numel(getArray(key))==1) && (numel(set) ==1)
       continue
   end
   equals=cell(0,1);
   equals{1}=getArray(key);
   for inner_index=1:numel(greater.keys)
       if skips(inner_index)
           continue
       end
        other_key=keys{inner_index};
        if isSubset(getArray(other_key),set) && isSubset(getArray(key), greater(other_key))
            equals{end+1}=getArray(other_key);
            skips(inner_index)=1;
        end
        
   end
   %cell2mat(equals)
   %equals 
   for el=1:numel(equals)
      fprintf("[%s]", join(check_names(equals{el}),","));
      if el~=numel(equals)
          fprintf("=");
      end
      set=setdiff(set,equals{el},'stable');
   end
   if numel(set) >0
    fprintf("> [%s]\n",join(check_names(set),",")); 
   else
    fprintf("\n");
   end
end
end




function key=getKey(indices)
key=vec2str(indices,[],[],0);
end

function array=getArray(key)
    ar=split((key),',');
    array=[cellfun(@str2num,ar)];
end
