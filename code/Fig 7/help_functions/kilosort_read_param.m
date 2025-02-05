function param = kilosort_read_param(params_file, param_str, param_type)
% kilosort_load_params(params_file, param_str)
% load param from params.py
param = [];
fid = fopen(params_file,'r'); 
while ~feof(fid)
    line = fgetl(fid); % read line by line
    param = sscanf(line, [param_str ' = ' param_type]);
    if ~isempty(param)
        break;
    end
end
fclose(fid);

end
