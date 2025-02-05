function [unit_id, unit_quality, is_noise] = kilosort_read_cluster_group(filename)
% [unit_id, unit_quality, is_noise] = kilosort_read_cluster_groups(filename)
% read cluster quality information

fid = fopen(filename);
C = textscan(fid, '%s%s');
fclose(fid);

unit_id = cellfun(@str2double, C{1}(2:end));
unit_quality = C{2}(2:end);
is_noise = cellfun(@(x)strcmp(x,'noise'),unit_quality);
