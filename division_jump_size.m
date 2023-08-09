%division_jump_size: estimates the mean division jump size using the data
%provided in the directory fetch_dir below. Also gives a vector of jump
%sizes for each dividing cell across all videos.

%Created by: Shahzeb Raja Noureen
%Date created: 11/12/2020
%Last modified: 11/12/2020

clear
close all

%directory address where data folders are located
fetch_dir=pwd;

%total number of data to read and use to calculate estimate
R=8;

%coordiantes of a cell after division
jump_coords=[];

%initialise the size jump size to 0
jump_size=0;

%preallocate an empty array for recording the jump sizes
jump_size_vec=[];

%tolerate for comapring flaots
tol=1e-5;

%total number of daughter tracks
num_tracks_all=0;

%for each video 1:R
for r=1:R
%%data prep    
%data import
fetch_folder_d=sprintf('/Data/Video %d/Daughter/Results.csv',r);

%full file name od daughter cell
full_file_name_d=fullfile(pwd, fetch_folder_d);

df_d=readtable(full_file_name_d,'Range','C:I');

%converting pixel coordianates to microns by multiplying by scale factor
%0.619
X_full_d=df_d.X*0.619;
Y_full_d=df_d.Y*0.619;
d_tracks=df_d.Track;

%convert the cell into chars (char)
tracks_str=string(d_tracks);

%daughter cell track id
track_id_d=unique(d_tracks);

%daughter cell track id in string format
track_id_d=string(track_id_d);

%number of daughter cell tracks
num_tracks=length(track_id_d);

%total number of daughter cell tracks
num_tracks_all=num_tracks_all+num_tracks;

%adding a 0 to tracks single digit track ids (1a,1b,etc) so that they can
%be sorted in the ascending order (1a,1b,2a,2b)
for i=1:num_tracks
    if size(char(track_id_d(i)),2)==2
        track_id_d(i)="0"+track_id_d(i);
    end
end
%sort the track_id into ascending order (1a,1b,2a,2b,...)
track_id_d=sortrows(track_id_d);


%do the above with track_str also
len_tracks_str=size(tracks_str,1);

for i=1:len_tracks_str
    if size(char(tracks_str(i)),2)==2
        tracks_str(i)="0"+tracks_str(i);
    end
end


%%mother cells data
fetch_folder_m=sprintf('/Data/Video %d/Mother/Results.csv',r);

%full file directory for the mother cell data
full_file_name_m=fullfile(fetch_dir, fetch_folder_m);

%reading columns C:I of mother data
df_m=readtable(full_file_name_m,'Range','C:I');

%converting pixel coordianates to microns by multiplying by scale factor
%0.619
X_full_m=df_m.X*0.619;
Y_full_m=df_m.Y*0.619;
tracks_m=df_m.Track;

%%analysis of the data

%clear the data frames df_d and df_m to save memory
clear df_d  df_m;


%for each track 1:num_tracks
for i=1:num_tracks
    
    %daughter tracks
    j=track_id_d(i);
    X_d=X_full_d(tracks_str==j);
    Y_d=Y_full_d(tracks_str==j);
    
    
    %length of daughter track
    lenX=length(X_d);
    
    %covert j (string) to char
    jchar=char(j);
    
    %extract the track num and store it as k
    k=str2double(jchar(1:end-1));
    
    %fetch data for track k
    X_m=X_full_m(tracks_m==k);
    Y_m=Y_full_m(tracks_m==k);
    
    
    %to deal with data inconsistency
    %if frameof mitosis is repeated in daughter data and mother data then
    if sqrt((X_m(end)-X_d(1))^2+(Y_m(end)-Y_d(1))^2)<tol
        
        pos_mitosis=[X_d(1) Y_d(1)];
        
         %jump size in in components dx and dy
        jump_coords=[(X_d(1)-X_d(2)),(Y_d(1)-Y_d(2))];

    else        
        pos_mitosis=[X_m(end) Y_m(end)];
        %jump size in in components dx and dy
        jump_coords=[(X_m(end)-X_d(1)),(Y_m(end)-Y_d(1))];
        
    end
    
    %calculate the mitotic jump size
    jump_size=sqrt(jump_coords(1)^2+jump_coords(2)^2);
    
    %check if the jump_size is 0 and find out where (which data, mother
    %and which daughter it is). Print the video number, mother id and the
    %daughter id of the cell with 0 jump size
%     if jump_size==0
%         fprintf('%d,%d,%s\n',r,k,j);
%     end
    
    %allocate the jump size to the jump size vector
    jump_size_vec=[jump_size_vec,jump_size];
    
end
end

%calcualte the total number of divisions (number of daughter cells/2)
num_divisions_all=num_tracks_all/2;

%mean jump size
mean_jump_size=mean(jump_size_vec)