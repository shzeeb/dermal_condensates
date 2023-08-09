%plots the distribution of initial distance of a cells (dividing or non-dividing
%from the nearest follicle given the cells get adsorbed eventually.

% close all;
clear;

type=["ND" "D"];

figure;

% t = tiledlayout(1,1);
for i=1:length(type)
    %non-dividing cells
    file_name=sprintf('pdf_dist_v2_type_%s.mat',type(i));
    
    %import data and variables
    data=open(file_name);

    r_vec=data.r_vec;
    norm_bins=data.norm_bins;
    
    A=trapz(r_vec(2:end),norm_bins);
    norm_dens=norm_bins/A;

    %plot the density
    

    if i==1
%         yyaxis left
        plot(r_vec(2:end),norm_dens,'-','LineWidth',2,'Color',[24/255, 163/255, 63/255],'DisplayName','non-dividing cells');
        hold on

        xlabel('Initial distance from follicle (microns)');
        ylabel('Probability density','Color','k');

%         df=[r_vec(2:end)',norm_dens'];
%         writematrix(df,'non_dividing_sim.csv','Delimiter','comma');


    elseif i==2
%         yyaxis right
        plot(r_vec(2:end),norm_dens,'-','LineWidth',2,'Color',[62/255, 79/255, 188/255],'DisplayName','dividing cells');
        
        hold off
%             set(gca,'YColor',[24/255, 163/255, 63/255])
%         df=[r_vec(2:end)',norm_dens'];
%         writematrix(df,'dividing_sim.csv','Delimiter','comma');

    end
%     set(gca,'Box','off');
    set(gca,'FontSize',16);


    xlim([35,100]);
    xticks(35:20:100);
    ylim([0,0.1]);


%     set(gcf,'renderer','Painters')
%     %name and save figure as an eps graphic
%     fig_name="distance_density_"+type(i);
%     print(fig_name,'-depsc');
end

 

legend