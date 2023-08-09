function [theta_vec,mu,var_angles,d_mean]=persistence_angles_func(len_pers_track)
% Analyses persistence angles.
%outputs: theta_vec: vector of persistence angles
% mu: mean direction of movement
% var_angles: variance in the direction of movement
% d_mean: mean distance moved per cell per persistent jump

%Created by: Shahzeb Raja Noureen
%Created on: 30/11/2020
%Last modified: 30/11/2020

fetch_dir=pwd;

%number of video files to use
R=8;

%to contain angle differences for all tracks as a flattened
%array
theta_vec=[];

%displacement vector
d_vec=[];


%for each video 1:R
for r=1:R
    
%%data import
fetch_folder_d=sprintf('Data/Video %d/Daughter/Results.csv',r);

full_file_name_d=fullfile(fetch_dir, fetch_folder_d);

df_d=readtable(full_file_name_d,'Range','C:I');

%converting pixel coordianates to microns%converting pixel coordianates to microns by multiplying by scale factor
%0.619
X_full_d=df_d.X*0.619;
Y_full_d=df_d.Y*0.619;
d_tracks=df_d.Track;

%convert the cell into chars (char)
tracks_str=string(d_tracks);

%daughter cell unique track id
track_id_d=unique(d_tracks);

%daughter cell track id in string format
track_id_d=string(track_id_d);

%number of daughter cell tracks
num_tracks=length(track_id_d);

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
fetch_folder_m=sprintf('Data/Video %d/Mother/Results.csv',r);

%full file directory for the mother cell data
full_file_name_m=fullfile(fetch_dir, fetch_folder_m);

%reading columns C:I of mother data
df_m=readtable(full_file_name_m,'Range','C:I');

%converting pixel coordianates to microns by multiplying by scale factor
%0.619
X_full_m=df_m.X*0.619;
Y_full_m=df_m.Y*0.619;
tracks_m=df_m.Track;



% num_tracks=2;


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
    if X_m(end)==X_d(1) && Y_m(end)==Y_d(1)
        
        %if track length >= 19 frames (3rs)
        if lenX>=len_pers_track+1
            %take persistent track from 2:19
            pers_X=X_d(2:len_pers_track+1);
            pers_Y=Y_d(2:len_pers_track+1);
    
        %otherise if the track is smaller than 18 frame (the cells may have
        %left the field of vision before the 18th frame or entered a follicle)
        else

            pers_X=X_d(2:end);
            pers_Y=Y_d(2:end);
        end
        
        
         %jump size in in components dx and dy
        jump_coords=[(X_d(1)-X_d(2)),(Y_d(1)-Y_d(2))];

    else
        if lenX>=len_pers_track
            
            %set persistence track as the first 18 frame (corresponding to
            %3 hours of motion)
            pers_X=X_d(1:len_pers_track);
            pers_Y=Y_d(1:len_pers_track);
    
        %otherise if the track is smaller than 18 frame (the cells may have
        %left the field of vision before the 18th frame or entered a follicle)
        else
            
            %persistent track
            pers_X=X_d;
            pers_Y=Y_d;
        end
        
        %jump size in in components dx and dy
        jump_coords=[(X_m(end)-X_d(1)),(Y_m(end)-Y_d(1))];
        
    end
    
    %division angle
    angle_mitosis=atan2(jump_coords(2),jump_coords(1));    
    
    %displacment vectors for persistent tracks
    pers_dx=(pers_X(1:end-1)-pers_X(2:end));
    pers_dy=(pers_Y(1:end-1)-pers_Y(2:end));
    
    %angle of displacement
    t=atan2(pers_dy,pers_dx);
    
    %displacement
    d=sqrt(pers_dx.^2+pers_dy.^2);
    
    d_vec=[d_vec;d];
    
    %angle of displacmenet for non-zero displacment only
    t(d==0)=[];
    
    %angle of derivaion from division angle
    theta=t-angle_mitosis;
    
    theta_vec=[theta_vec;theta];
    
end

end

%map theta_nz over the range [-pi,pi]
theta_vec(theta_vec>pi)=mod(theta_vec(theta_vec>pi),-2*pi);
theta_vec(theta_vec<-pi)=mod(theta_vec(theta_vec<-pi),2*pi);

%mean components of the direction of movement
Sbar=mean(sin(theta_vec));
Cbar=mean(cos(theta_vec));

%mean direction
mu=atan(Sbar/Cbar);

Rbar=sqrt(Sbar.^2+Cbar.^2);

%variance
var_angles=1-Rbar;


%mean displacement between frames during persistence (use non-zero
%displacements only to calculate the mean displacement)
d_mean=mean(d_vec(d_vec~=0));