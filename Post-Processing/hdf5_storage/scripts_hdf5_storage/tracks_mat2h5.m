function tracks_mat2h5(filepath)

% tracks_mat2h5(filepath)
%
% 03/2020 - Thomas Basset
%
% Convert track files as "tracks_sample.mat" to .h5 files
% ___________________________________________
% INPUT
% filepath : path to .mat file (without .mat)
% ___________________________________________

%load track structure
tracks=load([filepath '.mat']);
field=fieldnames(tracks);
tracks=tracks.(field{1});

%save in .h5
fields=fieldnames(tracks);
for k=1:length(fields) %for each field
    fieldk=fields{k};
    tracks_temp=vertcat(tracks.(fieldk));
    h5create([filepath '.h5'],['/' fieldk],[length(tracks_temp) 1],'ChunkSize',[1000 1],'Deflate',9)
    h5write([filepath '.h5'],['/' fieldk],tracks_temp)
end

end