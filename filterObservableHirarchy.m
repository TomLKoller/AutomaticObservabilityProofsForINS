function filtered_obs_from = filterObservableHirarchy(obs_from)
filtered_obs_from=cell(numel(obs_from),1);
for outer_index=1:size(obs_from,1)
    new_list=cell(0,1);
    for index=1:numel(obs_from{outer_index})
        is_child=0;
        for inner_index=1:numel(obs_from{outer_index})
            %check if a subset exists
            if (inner_index ~=index) && isSubset(obs_from{outer_index}{inner_index}, obs_from{outer_index}{index})
                is_child=1;
                break
           end
        end
        %Append if no parent has been found
        if ~is_child
            new_list{end+1}=obs_from{outer_index}{index};
        end
    end
    filtered_obs_from{outer_index}=new_list;
end
end


