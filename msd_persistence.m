%msd_persistence.m: Calculates mean square displacement for cells after
%division to check for persistence. Track length up to 18 length
%corresponding to 180 minutes (or 3 hrs) is used as is the suggested
%duration of persistence according to data.

%Coded by: Shahzeb Raja Noureen
%Date created: 11/2020
%Last modified: 24/01/2021

close all
clear

M=8; %8 datasets in total

%msd ensemble vector
msd_ens=[];

%time of time for persisent movement = 18*10=180 minutes (10 is the time
%between successive observations of position)
len_max_track=18;

%occurance counts
c=zeros(1,len_max_track-1);

%number of all tracks
num_tracks_all=0;

%loop through each data set (video)
for m=1:M
    
    full_file_name_d=pwd+"/Data/Video "+num2str(m)+"/Daughter/Results.csv";

    %read daughter cell data coloumn C:I
    df_d=readtable(full_file_name_d,'Range','C:I');

    %converting pixel coordianates to microns%converting pixel coordianates
    %to microns by multiplying by scale factor 0.619
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
    %full file directory for the mother cell data
    full_file_name_m=pwd+"/Data/Video "+num2str(m)+"/Mother/Results.csv";

    %reading columns C:I of mother data
    df_m=readtable(full_file_name_m,'Range','C:I');

    %converting pixel coordianates to microns by multiplying by scale factor
    %0.619
    X_full_m=df_m.X*0.619;
    Y_full_m=df_m.Y*0.619;
    tracks_m=df_m.Track;

    %time between frames
    dt=10;
    
    %length of track_vec
    num_tracks=length(track_id_d);
    
    %update the total number of tracks
    num_tracks_all=num_tracks_all+num_tracks;
    
    %for foe each track, calculate the msd
    for j=1:num_tracks
        
        %daughter tracks
        i=track_id_d(j);
        X_d=X_full_d(tracks_str==i);
        Y_d=Y_full_d(tracks_str==i);
        
        %covert j (string) to char
        jchar=char(i);

        %extract the track num and store it as k
        k=str2double(jchar(1:end-1));

        %fetch data for track k
        X_m=X_full_m(tracks_m==k);
        Y_m=Y_full_m(tracks_m==k);
        
        %length of daughter track
        lenX=length(X_d);
        
        
        len_pers_track=len_max_track;
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

        end
        
        %set X and Y persistent tracks
        X=pers_X';
        Y=pers_Y';
        
        %number of point in trajectory
        N=length(X);
        
        %maximum lag
        max_lag=N-1;
        
        msd_time_ave=zeros(1,len_max_track-1);
        
        
        % for each value of time lag
        for lag=1:max_lag
            %temporarily msd variable
            msd_temp_sq=(X(lag+1:end)-X(1:end-lag)).^2+(Y(lag+1:end)-Y(1:end-lag)).^2;
            
            %time average
            msd_time_ave(lag)=sum(msd_temp_sq(1:N-lag));
            
            %occurance counts
            c(lag)=c(lag)+N-lag;
            
        end
        
        msd_ens=[msd_ens; msd_time_ave];
        
    end
end
    
%calculate the msd time-ensemble average
msd_time_ens_ave=sum(msd_ens,1)./c;

%lag vector
lag_vec=1:(len_max_track-1);

%lagged time
lagged_dt=lag_vec*dt;

%fit a line throuhg logged time and msd data
p=polyfit(log(lagged_dt),log(msd_time_ens_ave),1);

%the slope of the fitted line
alpha=p(1);

%record the hurst's coefficient
hurst_coeff=alpha/2;

mat_file_name="msd_persistence_"+m;
% save(mat_file_name+".mat",'msd_time_ens_ave','lagged_dt');

loglog(lagged_dt,msd_time_ens_ave,'o-');
hold on
loglog(lagged_dt,exp(p(2))*lagged_dt.^alpha,'x');


figure;

ax=gca;
green=[0, 166/255, 81/255];

plot(lagged_dt,msd_time_ens_ave,'.','MarkerSize',20);
hold on
plot(lagged_dt,exp(p(2))*lagged_dt.^alpha,'-','Color',green,'LineWidth',2);
legend('MSD data','fitted values','Location','northwest');
xlabel('time (mins)');
ylabel('MSD (\mu m^2)');

ax.FontSize=14;

fprintf("Hurst coeffceint = %f \n",hurst_coeff);
%save the hurst coefficients
% save("hurst_coefficients.mat",'hurst_coeff');
