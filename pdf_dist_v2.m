%analysis of the distances between cells and follicles for divided and
%removed cells and for cells which did not divide and got adsorbed. Uses
%the technique for finding the minimum distance rather than considering
%disjoint rings separating each follicle region from another.
%End result: a probability distirbution of starting a distance r away from a
%follicle conditional on the fact that the cell divided (or did not divide)
%and got adsorbed later on.

%Created by: Shahzeb Raja Noureen
%Date created: 14/12/2020
%Last modified: 15/01/2020

%clear all variables and close all open figures
clear
close all

%import data
data1=open("dermal_condensates_with_persistence_1.mat");
data2=open("dermal_condensates_with_persistence_2.mat");
data3=open("dermal_condensates_with_persistence_3.mat");
data4=open("dermal_condensates_with_persistence_4.mat");


X_init_mat=[data1.X_init_mat data2.X_init_mat data3.X_init_mat data4.X_init_mat];
Y_init_mat=[data1.Y_init_mat data2.Y_init_mat data3.Y_init_mat data4.Y_init_mat];



%dividers or non-dividers. Insert 'D' for dividers and 'ND' for
%non-dividers
type='D';

%number of follicles on lattice
num_fol=data1.num_fol;

%size of domain
Lx=data1.Lx;
Ly=data1.Ly;

%number of cells
% N=data.N;

N=size(data1.removed_mat,1);

M_vec=[data1.M data2.M data3.M data4.M];
M=sum(M_vec);

%matrix containing indicators whether a cells adsorbed or not for each
%repeat
removed_mat=[data1.removed_mat data2.removed_mat data3.removed_mat data4.removed_mat];

%matrix containing indicators whether a cells was iniitially adsorbed or
%not for each repeat
removed_init_mat=[data1.removed_init_mat data2.removed_init_mat data3.removed_init_mat data4.removed_init_mat];
% removed_init_mat=[data.removed_init_mat;zeros(N-data.N_init,M)];

%martrix containing indicators whether a cell proliferated or not for each
%repeat
prolif_mat=[data1.prolif_mat data2.prolif_mat data3.prolif_mat data4.prolif_mat];

%follicle equations
fplus=data1.fplus;
fminus=data1.fminus;

%radius of follicle
r_0=data1.r_0;

%separation distance between follicles
d_sep=data1.d_sep;

%logical. set to 1 if need ghost follicles. Otherwise set to 0
ghost=0;

%generate centre coordinates of follicle centres
[xfc,yfc]=follicles(Lx,Ly,num_fol,r_0,d_sep,ghost);


%if ghost follicles enabled then find the new number of follicles including the
%ghost follicles
if ghost==1
    num_fol=size(xfc,1);
end

%prealloc empty distance distace vector
dist_vec=[];

%prealloc memory for a distance matrix
dist_mat=zeros(N,num_fol,M);


%for each repeat
for m=1:M
    
    %fetch initial positions of all cells for m-th repeat
    X_init=X_init_mat(:,m);
    Y_init=Y_init_mat(:,m);
    
    %if type==D then execute code for dividers
    switch type
        case 'D'   
            for k=1:num_fol   
                %calculate teh distance between each follicle and each cell
                %in turn. Only save values for dividers which proliferated,
                %were not initially adsorbed and which got adsorbed at a
                %later stage during simulation
                dist_=sqrt((X_init-xfc(k)).^2+(Y_init-yfc(k)).^2).*(prolif_mat(:,m)~=0).*(removed_init_mat(:,m)==0).*(removed_mat(:,m)==1);
                
                %do the above but for the distance matrix
                dist_mat(:,k,m)=dist_;%.*(prolif_mat(:,m)~=0).*(removed_init_mat(:,m)==0).*(removed_mat(:,m)==1);
                
                %save results of dist_ inside the dist_vec array
%                 dist_vec=[dist_vec;dist_];
            end
        
        %otherwise if type==ND then execute code for non-dividers
        case 'ND'
            for k=1:num_fol
                %calculate teh distance between each follicle and each cell
                %in turn. Only save values for cell which did not proliferate,
                %were not initially adsorbed and which got adsorbed at a
                %later stage during simulation
                dist_=sqrt((X_init-xfc(k)).^2+(Y_init-yfc(k)).^2).*(prolif_mat(:,m)==0).*(removed_init_mat(:,m)==0).*(removed_mat(:,m)==1);
                
                %do the above but for the distance matrix
                dist_mat(:,k,m)=dist_.*(prolif_mat(:,m)==0).*(removed_init_mat(:,m)==0).*(removed_mat(:,m)==1);
                
                %save results of dist_ inside the dist_vec array
%                 dist_vec=[dist_vec;dist_];
            end
    end 
end



%prealloc memory to contain minimum distances between each cell and a
%follicle
min_dist=[];

%follicle index for corresponding to the minimum distance
idx=[];

%for each repeat, search for the minimum distance between each cell and
%each follicle. Save the minimum distances and the follicle index for that
%minimum distance.
for m=1:M
    dist_vec=dist_mat(:,:,m);
    [min_dist_temp,idx_temp]=min(dist_vec,[],2);
    min_dist=[min_dist; min_dist_temp];
    idx=[idx;idx_temp];
end

%keep non-zero minimum dist and the corresponding follicle index
min_dist_nz=min_dist(min_dist~=0);
idx_nz=idx(min_dist~=0);

%binning the data
%width of each annulus
dr=1;

%number of annuli (bins)
P=ceil((max(min_dist_nz)-r_0)/dr);

%radii of each annulus relative to the radius of its follicle
r_vec=r_0:dr:P*dr+r_0;

% figure;
% histogram(min_dist_nz,P,'BinLimits',[35 max(min_dist_nz)])


%divding data into the the different annuli and normalising it by the area
%of the annulus
%prealloc empty matrix for binned data
binned_data=[];

%max size of binned data
max_size=0;

%prealloc vector for counts inside each bin
bin_counts=zeros(1,P);

%for each annulus (p)
for p=1:P
    
    %fetch distance data to go inside bin p
    bin_data=min_dist_nz(min_dist_nz>=r_vec(p) & min_dist_nz<=r_vec(p)+dr);
    
    %length of bin data (num rows)
    len_bin_data=size(bin_data,1);
    
    
    %if size of new data smaller than the current max of data
    if max_size>len_bin_data
        
        bin_data=[bin_data;zeros(max_size-len_bin_data,1)];
    
    %otherwise if new binning data is larger than current maximum
    elseif max_size<len_bin_data
        num_zeros_r=len_bin_data-max_size;
        num_zeros_c=size(binned_data,2);
        
        binned_data=[binned_data;zeros(num_zeros_r,num_zeros_c)];
    end

    %record binning data in the binned_matrix 
    binned_data(:,p)=bin_data;
    
    %calculate the bin counts
    bin_counts(p)=nnz(binned_data(:,p));
    
    %update the current maximum data size
    max_size=max(len_bin_data,max_size);
    
end

%areas of each annulus
area_rings=pi*(r_vec(2:end).^2-r_vec(1:end-1).^2);
% circumf=2*pi*r_vec;

% area_2=2*pi*r_vec*dr;

%normalised bin counts for distances inside each annulus ring
norm_bins=bin_counts./area_rings;

%save the data as a .mat file for plotting later
file_name=sprintf('pdf_dist_v2_type_%s.mat',type);
save(file_name);



%% time to absorptions

t_removed=[data1.t_removed_mat data2.t_removed_mat data3.t_removed_mat ...
    data4.t_removed_mat];
t_rem_div=t_removed.*(prolif_mat~=0);
t_rem_nondiv=t_removed.*(prolif_mat==0);


% figure;
% histogram(t_rem_nondiv(t_rem_nondiv~=0));
% 
% title('Non-dividing cells cells','FontSize',18,'Interpreter','latex');
% xlabel('Adsorption time','FontSize',20,'Interpreter','latex');
% ylabel('Frequency','FontSize',20,'Interpreter','latex');
% 
% figure;
% histogram(t_rem_div(t_rem_div~=0));
% 
% title('Dividing cells cells','FontSize',18,'Interpreter','latex');
% xlabel('Adsorption time','FontSize',20,'Interpreter','latex');
% ylabel('Frequency','FontSize',20,'Interpreter','latex');

figure;
h=histogram(t_removed(t_removed~=0),'Normalization','pdf');
h.FaceColor=[190/255, 191/255, 193/255];

% title('All cells','FontSize',18,'Interpreter','latex');
xlabel('adsorption time');
ylabel('probability density');
ax=gca;
ax.FontSize=16;

set(gcf,'renderer','Painters')

%name and save figure as an eps graphic
fig_name="tta_all_cells";
% print(fig_name,'-depsc');
