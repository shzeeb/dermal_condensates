%msd_non_dividers: Estimate the diffusion coefficient by calculating the
%msd for all the tracks in each data set.
%variables/objects: 'lagged_dt','msd_time_ens_ave',
%'mean_D' (in part 2),'CI_D','CI_lower_msd','CI_upper_msd','fitted_vals','model' are
%saved once during each looping and the variable 'D' is saved once at the
%end.

%Coded by: Shahzeb Raja Noureen
%Created on: 24/01/2020
%Last modified: 24/01/2020

%this code uses amalgamates all the track to calcualte the overall msd. Then
%fit a line of best fit to this data and estimate the overall (mean)
%diffusion coefficient

close all
clear

M=8; %8 videos in total

%initialise amalgamated track vectors
X_rec_full=[];
Y_rec_full=[];

%initialse track id vector
track_id=[];

%temporary value of the old track number
track_num_temp=0;

%pull all the data in one matrix
for m=1:M
    %full path to load cell tracking data
    full_file_name=pwd+"/Data/Video "+num2str(m)+"/Non-dividing/Results.csv";
    
    %read the data as table
    df=readtable(full_file_name);
    
    %converting pixel coordianates to microns%converting pixel coordianates
    %to microns by multiplying by scale factor 0.619
    X_full=df.X*0.619;
    Y_full=df.Y*0.619;

    track_id=[track_id;df.Track+track_num_temp];
    track_num_temp=max(track_id);
    
    X_rec_full=[X_rec_full;X_full];
    Y_rec_full=[Y_rec_full;Y_full];
end

X_full=X_rec_full;
Y_full=Y_rec_full;

%time between frames
dt=10;

%find the unique tracks
track_vec=unique(track_id);

%length of track_vec
num_tracks=length(track_vec);

%finding the largest track size

%length of the largest track, initialised to 0
len_max_track=0;

%size of each tracks: for determning what acceptable max lag should be.
len_tracks=zeros(1,num_tracks);

%find the length of the largest track
for i=1:num_tracks
    j=track_vec(i);
    len_tracks(i)=sum(track_id==j);
end

len_max_track=max(len_tracks);

%msd ensemble vector
msd_ens=[];

%occurance count vector
c=zeros(1,len_max_track-1);

%for each track calculate the msd
for j=1:num_tracks

    %track number
    i=track_vec(j);

    %track coordinates for track i
    track_X=X_full(track_id==i);
    track_Y=Y_full(track_id==i);

    %length of each track j for each video m
%     len_track(j)=length(track_X);
    
    X=track_X';
    Y=track_Y';

    %number of point in trajectory
    N=length(X);

    %maximum lag
    max_lag=N-1;

    msd_time_ave=zeros(1,len_max_track-1);


    % for each value of time lag
    for lag=1:max_lag
        msd_temp_sq=(X(lag+1:end)-X(1:end-lag)).^2+(Y(lag+1:end)-Y(1:end-lag)).^2;

        %time average
        msd_time_ave(lag)=sum(msd_temp_sq(1:N-lag));
        c(lag)=c(lag)+N-lag;
    end
    
    %msd ensemble
    msd_ens=[msd_ens; msd_time_ave];
end

%msd time-ensemble average
msd_time_ens_ave=sum(msd_ens,1)./c;

lag_vec=1:(len_max_track-1);
lagged_dt=lag_vec*dt;

%fit a linear regression model
model=fitlm(lagged_dt,msd_time_ens_ave,'Intercept',false);

%extract coefficients of the regression model
model.coefCI;

%%calculate the variance of the estimated msd
% mean of lagged time
mean_lagged_dt=mean(lagged_dt);
x_minus_xbar_sq=(lagged_dt-mean_lagged_dt).^2;
var_lagged_dt=var(lagged_dt);
mse=model.MSE;

var_mean_msd=mse*(1/N+x_minus_xbar_sq/((N-1)*var_lagged_dt));

%fitted values
fitted_vals=model.Fitted;


%95% CI for the response (msd)
CI_lower_msd=fitted_vals'-2*sqrt(var_mean_msd);
CI_upper_msd=fitted_vals'+2*sqrt(var_mean_msd);

model_coeff=model.Coefficients;

%mean slope (msd)
mean_slope=table2array(model_coeff(1,1));

%coefficients of the model
coeff_CI=model.coefCI;
%CI for the slope
CI_slope=coeff_CI;

% D
D=mean_slope/4;

%CI D
CI_D=CI_slope/4;

ax=gca;

green=[0, 166/255, 81/255];

plot(lagged_dt,msd_time_ens_ave,'.');
hold on
plot(lagged_dt,model.Fitted,'-','Color',green,'LineWidth',2);

xlim([0 max(lagged_dt)]);
xlabel('time (mins)');
ylabel('MSD (\mu m^2)');

legend('MSD data','Fitted values');
ax.FontSize=14;



disp("Diffusion coefficient, D = "+D);

% save("msd_non-dividers_"+m+".mat",'lagged_dt','msd_time_ens_ave','CI_D','CI_lower_msd','CI_upper_msd','fitted_vals','model')
% save("diffusion_coefficient_non-dividers.mat",'D');
