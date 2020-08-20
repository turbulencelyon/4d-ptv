%% Convert "tracks_sample.mat" to "tracks_sample.h5"
filepath='C:\Users\Thomas Basset\Desktop\scripts\hdf5_storage\tracks_sample';
tracks_mat2h5(filepath);

%% Read "tracks_sample.h5" and return the same structure as "tracks_sample.mat"
tracks=tracks_h52mat(filepath);

%% To read given fields
fields={'x','y','z','t'};
tracks=tracks_h52mat(filepath,fields);
