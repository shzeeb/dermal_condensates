%plot mean proportions and 95% confidence interval for 8 randomly chosen
%repeats

clear
close all

%to plot experimental proportions or simulated ones
simulated=1;

if simulated==0
    
    
    %proportions experimental data
    data=[16 27.27272727 16.66666667 18.42105263 27.08333333 16.12903226 5.454545455 27.90697674;...
        26.66666667 25 4.651162791 5.128205128 22.72727273 4.444444444 13.04347826 12.76595745];
    
    %mean proportions
    props=[mean(data(1,:)),mean(data(2,:))];
    
    len_data=size(data,2);
    
    %standard error in proportions
    std_err=[sqrt(var(data(1,:))/len_data) sqrt(var(data(2,:))/len_data)];
    
    %upper and lower limit of CI
    CI_upper=std_err; %times by 2.306 for 95% CI (value given for N=8 and with N-1=7 DOF use the t distribution)
    CI_lower=std_err;
    
else 

    %load data: you must load the .mat file containing the proportions data
    %from the simulations. The names chosen are: 
    % (1) dermal_condensates_no_progeny
    % (2) dermal_condenate_v9_no_persistence
    % (3) dermal_condensates_with_persistence
    %if you change the name of the saved files then you will have to call
    %it by the new name here.
    
    file_name="dermal_condensates_with_persistence";

    data_temp_1=[];
    data_temp_2=[];

    num_data=4;

    for i=1:num_data
        df=load(file_name+"_"+num2str(i)+".mat");
        
        len_props=size(df.prop_div_cells,2);
        
        idx=randperm(len_props,8);
        
        data_temp_1=[data_temp_1 df.prop_div_cells(idx)];
        data_temp_2=[data_temp_2 df.prop_nondiv_cells(idx)];
    end

    data=[data_temp_1;data_temp_2]*100;


    props=mean(data,2);

    len_data=size(data,2);


    std_err=[sqrt(var(data(1,:))/len_data) sqrt(var(data(2,:))/len_data)];


    CI_upper=2.306*std_err;
    CI_lower=2.306*std_err;
end


% plot proportions


figure;
%plot bar chart
b1=bar(1,props(1));
hold on
b2=bar(2,props(2));
hold on

nbars=2;
ngroups=2;

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar([1 2], props, CI_upper, CI_lower, '.r','LineWidth',1);
end


%set edge colour to none
b1.EdgeColor='none';
b2.EdgeColor='none';

%face colours
light_grey=[190/255, 191/255, 193/255];
dark_grey=[77/255, 77/255, 79/255];

%set set colour
b1.FaceColor=light_grey;
b2.FaceColor=dark_grey;

%x axis ticks
names = {'Dividing cells';'Non-dividing cells'};

%set x axis ticks
set(gca,'xtick',1:2,'xticklabel',names)

%y label
ylabel('Proportion %');

%show the proportions as text on top of each bar
xtips1 = b1.XEndPoints;
ytips1 = b1.YEndPoints;
labels1 = string(b1.YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',16)
xtips2 = b2.XEndPoints;
ytips2 = b2.YEndPoints;
labels2 = string(b2.YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',16)

ax=gca;
ax.FontSize=16;
ylim([0 20]);

% figure; histogram(df.prop_div_cells,'Normalization','pdf')
% figure;
% histogram(df.prop_nondiv_cells,'Normalization','pdf')

set(gcf,'renderer','Painters')

%name and save figure as an eps graphic

if simulated==0
    fig_name="observed_proportions";
    
else
fig_name=file_name;
end

error_bars_div=[props(1)-CI_lower(1) props(1)+CI_upper(1)];
error_bars_nondiv=[props(2)-CI_lower(2) props(2)+CI_upper(2)];
error_bars=[error_bars_div;error_bars_nondiv];

%save figure
% print(fig_name,'-depsc');

%save data as csv
% writematrix(error_bars,"error_bars_"+fig_name+".csv");
