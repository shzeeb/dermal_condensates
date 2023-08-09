%this scrript finds an optimal parameter kappa_star of von-Mises which
%minimises the sum of squared residuals between the empirical distribution
%and the fitted curve.

%Created by: Shahzeb Raja Noureen
%Created on: 30/11/2020
%Last modified: 30/10/2020

clear;

%grid of parameter kappa
kappa=0:0.01:3;

%length of kappa
len_kappa=length(kappa);

%%von-mises pdf
%besel function of first kind and order zero
I0 = @(k,r) 1/(factorial(r)^2)*(k/2)^(2*r);

%von-mises pdf
g = @(x,mu,kappa,sum_I0) 1/(2*pi*sum_I0)*exp(kappa*cos(x-mu));

%number of bins to use for histogram
nbins=50;

%support of the pdf (mesh x between -pi and pi)
x=linspace(-pi,pi,nbins);

%max track length we'd like to include for our estiamte (18 max)
len_track=18;


%sum of squared residuals
ssr=zeros(1,len_kappa);

%find theta_vec and mu (mean direction), variance of the angles and the mean displacement
%using the function persistence_func
[theta_vec,mu,var_angles,d]=persistence_angles_func(len_track);

%for testing the fit on simulated data
% mu=0.5;
% kappa_test=1.5;
% theta_vec=circ_vmrnd(mu,kappa_test,100000);


%do a histogram of theta_vec using nbins number of bins
h=histogram(theta_vec,nbins,'Normalization','pdf');


%for each value of kappa mesh
for j=1:len_kappa

    %evaluate the Besel function
    sum_I0=0;
    for r=0:1000
        sum_I0=sum_I0+I0(kappa(j),r);
    end

    %evaluate the VM pdf
    y=g(x,mu,kappa(j),sum_I0);

    %calculate the SSR
    ssr(j)=sum((y-h.Values).^2);

end

%find kappa which minimises the SSR
kappa_star=kappa(ssr==min(ssr));

T=table(mu,var_angles,kappa_star,d);
T.Properties.VariableNames={'mean direction (\mu)','variance of angles (\sigma^2)','optimal kappa (\kappa^\star)','mean displacement (d)'}


%% plotting histogram and optimal pdf


%do a histogram of theta_vec using nbins number of bins
h=histogram(theta_vec,20,'Normalization','pdf');
h.FaceColor=[190/255, 191/255, 193/255];
hold on
    
    
sum_I0=0;

for r=0:1000
    sum_I0=sum_I0+I0(kappa_star,r);
end

plot(x,g(x,mu,kappa_star,sum_I0),'-','Color',[0, 166/255, 81/255],'LineWidth',2);
xlabel('Persistence (deviation) angle (rad)');
ylabel('Probability density');

%set label and ticks font size
ax=gca;
ax.FontSize=16;
ax.XColor=[77/255, 77/255, 79/255];
ax.YColor=[77/255, 77/255, 79/255];

legend('Data','VM fit');
% set(gcf,'renderer','Painters')
% 
% print("persistence_histogram",'-depsc');
