%Division angle data analysis. Imports division angle data. Shifts data and
%mod(data,180) to map the angles to [0,pi]. Plots the histograms. Also fits
%a VM(mu,kappa) distribution to the amalgamated angle data using least
%squares.

%Created by: Shahzeb Raja Noureen
%Date created: 30/09/2020
%Last modified: 05/12/2020


clear
close all

%directory for reading the data

fetch_folder=pwd+"/Data/Division_angle_data";


%number of data files to read (max 7)
num_files = 7;

%angle vector to contain amalgamated angles
angles=[];


%vector for mean direction
mu=zeros(1,num_files);

%unadjusted (raw) division angles
raw_data=[];

%for each data file
for i = 1:num_files
    %case file name
    base_file_name = sprintf('video_%d_angle_length.csv',i);
    
    %full file name (directory/base_file_name)
    full_file_name = fullfile(fetch_folder, base_file_name);
    
    %read data from the file
    df = readtable(full_file_name);
    
    raw_data=[raw_data;df.Angle];
    
    %extract angles from data frame
    data=mod(df.Angle+180,180);
    
    
    %fit a kernel to the data density
    [kernel,theta]=ksdensity(data,linspace(0,180,1000));
    
    %find the maximum of the kernel curve
    max_kernel=max(kernel);
    
    %find angle (theta) at the maximum
    theta_max_kernel=theta(kernel==max_kernel);
    
    %shift the data by 90-theta_max_kernel degrees
    data=mod(data+90-theta_max_kernel,180);

    %append the data to the angles vector
    angles = [angles; data];
    
    
    %find Sbar and Cbar (see Mardia (2000), Directional statistics, p. 15)
    Sbar=mean(sin(data*pi/180));
    Cbar=mean(cos(data*pi/180));

    %find the mean direction (see Mardia (2000), Directional statistics, p. 15)
    if Cbar>=0
        mu(i)=atan(Sbar/Cbar);
    else
            mu(i)=atan(Sbar/Cbar)+pi;
    end
end

%% plot a histogram of the raw data
raw_data_rad=deg2rad(raw_data);
figure;
h=histogram(raw_data_rad,20,'Normalization','pdf');
h.FaceColor=[190/255, 191/255, 193/255];
xlabel('Division angle (rad)');
ylabel('Probability density');
ax=gca;
ax.FontSize=16;

set(gcf,'renderer','Painters')
saveas(ax,"raw_div_angles_hist"+".eps",'epsc');

%find the peaks if the raw angles distribution

[density,xi]=ksdensity(raw_data_rad);
peaks=findpeaks(density);

len_peaks=length(peaks);
mode_vec=zeros(len_peaks,1);

for i=1:len_peaks
    mode_vec(i)=xi(peaks(i)==density);
end


%% fitting von-Mises distribution minimising sum of squared residuals
%convert angles into radians
angles=angles*pi/180;

%find Sbar and Cbar
Sbar=mean(sin(angles));
Cbar=mean(cos(angles));

Rbar=sqrt(Sbar.^2+Cbar.^2);
var_angles=1-Rbar;

%find the overall mean
if Cbar>=0
    mu_all=atan(Sbar/Cbar);
else
    mu_all=atan(Sbar/Cbar)+pi;
end


%grid of parameter kappa
kappa=0:0.01:3;

%length of kappa
lenk=length(kappa);

%%von-mises pdf
%Besel function of first kind and order zero
I0 = @(k,r) 1/(factorial(r)^2)*(k/2)^(2*r);

%von-mises pdf
g = @(x,mu,kappa,sum_I0) 1/(2*pi*sum_I0)*exp(kappa*cos(x-mu));

%number of bins to use for histogram
nbins=50;

%support of the pdf (mesh x between -pi and pi)
x=linspace(0,pi,nbins);
%sum of squared residuals
ssr=zeros(1,lenk);

%find bin counts
bin_counts=histcounts(angles,nbins,'Normalization','pdf');


%for each value of kappa mesh
for j=1:lenk

    %evaluate the Besel function
    sum_I0=0;
    for r=0:1000
        sum_I0=sum_I0+I0(kappa(j),r);
    end
    
    %evaluate the VM pdf
    y=g(x,mu_all,kappa(j),sum_I0);

    %calculate the SSR
   ssr(j)=sum((y-bin_counts).^2);

end

%find kappa which minimises the SSR
kappa_star=kappa(ssr==min(ssr));

figure;
h=histogram(angles,20,'Normalization','pdf','FaceColor','b'); hold on
h.FaceColor=[190/255, 191/255, 193/255];

%evaluate the Besel function for kappa_star
sum_I0=0;
for r=0:1000
    sum_I0=sum_I0+I0(kappa_star,r);
end


%plots VM distribution curve
plot(x,g(x,mu_all,kappa_star,sum_I0),'-','Color',[0, 166/255, 81/255],'LineWidth',2);

% darkBackground(h);
%title the plot and label the axes
% title(sprintf('von-Mises with $\\mu=%1.3f$ and $\\kappa=%1.3f$',mu_all,kappa_star),'Interpreter','latex');

xlabel('Division angle (rad)');
ylabel('Probability density');
legend('Data','VM fit');
ax=gca;
ax.FontSize=16; 

ax.XColor=[77/255, 77/255, 79/255];
ax.YColor=[77/255, 77/255, 79/255];
% set(gcf,'renderer','Painters')

% print('divisions_histogram','-depsc');
% print('div_angles_adjusted_2','-depsc');