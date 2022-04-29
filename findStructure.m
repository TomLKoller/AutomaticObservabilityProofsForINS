meas_functions;

num_functions=numel(meas_functions_);
global obs_from;
obs_from=cell(num_functions,1);
for f_index=1:num_functions
    obs_from{f_index}=cell(0,1);
end
for f_index=1:num_functions
   traverseStructure([f_index],zeros(num_functions,1),meas_functions_,0,1);
end
  %filtered_obs_from= filterObservableHirarchy(obs_from);
   
function traverseStructure(function_indices, pred_result,meas_functions_,start_Omega,depth)
% if depth ==3 
%     return
% end
disp(function_indices);
if depth > 1
[check_result, new_Omega]=obs_proof(meas_functions_(function_indices(end)),start_Omega);
else 
[check_result, new_Omega]=obs_proof(meas_functions_(function_indices));
end
global obs_from;
for c_index=1:numel(check_result)
    if check_result(c_index) > pred_result(c_index)
        obs_from{c_index,1}{end+1}=function_indices;
    end
end
for c_index=function_indices(1):numel(check_result)
    if ~check_result(c_index)
        traverseStructure([function_indices, c_index],check_result,meas_functions_,new_Omega,depth+1);
    end
end
end