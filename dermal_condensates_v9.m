%Name: dermal_condensates_v9
%Description: cells diffused and divide and get adsorbed when hit with a
%hair follicle. At division, cells perform a jump over 10 mins. No
%persistent motion here.


clear;

% parameters not required to reallocate before each repeat
% initalise simulation parameters here:
%domain size
ax = 0;
Lx = 634.15;
ay = 0;
Ly = 634.15;
D = 0.078036;%0.1504;%0.15836;%0.2643;%0.413304;%0.2196;%0.5465;%0.01207175; % diffusion coefficient. WT D = 0.657
dt = 0.1; %dt needs to be sufficiently small such that the position jumps do not take a partile outside of the domain
radius_fol = 35; %radius of hair follicle (fixed for all follicles)
mean_dist_fol = 170; %mean distance between the centre of one follicle to the centre of an adjacent one
jump_mitosis = 0.64271; %jump size per agent at mitosis. Units: microns per minute
mu_mit=1.562;
kappa_mit=2.85;
test=0;
ghost=0;

if test==1
    tfin = 200; %
    M=1; %number of repeats
    num_iter = ceil(tfin/dt);
    N_init = 5; %initial number of agents on the domain. WT n0 = 1936
    num_fol = 0; %number of hair follicles
    num_prolif_total=4; %from Mingliang's thesis: 276 events from E13.5-E15
    num_non_div_cells=1;
    video=1;
    Save=0;
    
elseif test==0
    tfin = 1800; %
    M=500; %number of repeats
    num_iter = ceil(tfin/dt);
    N_init = 1936; %initial number of agents on the domain. WT n0 = 1936
    num_fol =14; %number of hair follicles
    num_prolif_total=623; %from Mingliang's thesis: 276 events from E13.5-E15
    num_non_div_cells=num_prolif_total;
    r_0=radius_fol;
    d_sep=mean_dist_fol;
    video=0;
    Save=1;
    
else
    error("'test' must be either 0 or 1.");
end

rec_step=10;

%times to record state of system at regular intervals
rec_times = [0:rec_step:tfin];
% rec_times=[rec_times,inf];

% rec_times=linspace(1,1800,4);
% rec_times=[600 1200 1800];

%number of points to plot the follcles (npts only used to plot follicles).
num_pts = 1000;

% x points to evaluate fplus and fminus at.
xpts = zeros(num_pts,num_fol);

%vectors to e populated by the y cooridinates of each hair follicle.
circ_plus = zeros(num_pts,num_fol);
circ_minus = zeros(num_pts,num_fol);

%equation of follicles (essentially the equation of a circle with centre
%with centre at (centre_fol_x,centre_fol_y) and radius radius_fol,
%rearranged
fplus = @(xpts,centre_fol_x,centre_fol_y,radius_fol) sqrt(radius_fol^2 - (xpts-centre_fol_x).^2) + centre_fol_y;
fminus = @(xpts,centre_fol_x,centre_fol_y,radius_fol) -sqrt(radius_fol^2 - (xpts-centre_fol_x).^2) + centre_fol_y;


if test==1
    if num_fol==1
        xfc(1)=Lx/2;
        yfc(1)=Ly/2;
    end
    
elseif test==0
    [xfc,yfc]=follicles(Lx,Ly,num_fol,r_0,d_sep,ghost);
    
else
    error("'test' must be either 0 or 1.");
end


% xfc(1)=radius_fol;
% xfc(2)=xfc(1)+mean_dist_fol;
% xfc(3)=xfc(1)+2*mean_dist_fol;
% xfc(4)=xfc(1)+3*mean_dist_fol;
% 
% xfc(5)=radius_fol+mean_dist_fol/2;
% xfc(6)=xfc(5)+mean_dist_fol;
% xfc(7)=xfc(5)+2*mean_dist_fol;
% xfc(8)=xfc(5)+3*mean_dist_fol;
% 
% xfc(9)=radius_fol;
% xfc(10)=xfc(1)+mean_dist_fol;
% xfc(11)=xfc(1)+2*mean_dist_fol;
% xfc(12)=xfc(1)+3*mean_dist_fol;
% 
% xfc(13)=radius_fol+mean_dist_fol/2;
% xfc(14)=xfc(5)+mean_dist_fol;
% xfc(15)=xfc(5)+2*mean_dist_fol;
% xfc(16)=xfc(5)+3*mean_dist_fol;
% 
% 
% 
% yfc(1)=radius_fol;
% yfc(2)=radius_fol;
% yfc(3)=radius_fol;
% yfc(4)=radius_fol;
% 
% yfc(5)=yfc(1)+mean_dist_fol;
% yfc(6)=yfc(1)+mean_dist_fol;
% yfc(7)=yfc(1)+mean_dist_fol;
% yfc(8)=yfc(1)+mean_dist_fol;
% 
% yfc(9)=yfc(1)+2*mean_dist_fol;
% yfc(10)=yfc(1)+2*mean_dist_fol;
% yfc(11)=yfc(1)+2*mean_dist_fol;
% yfc(12)=yfc(1)+2*mean_dist_fol;
% 
% yfc(13)=yfc(1)+3*mean_dist_fol;
% yfc(14)=yfc(1)+3*mean_dist_fol;
% yfc(15)=yfc(1)+3*mean_dist_fol;
% yfc(16)=yfc(1)+3*mean_dist_fol;


% xfc(1) = 170.75;
% xfc(2) = xfc(1) + mean_dist_fol;
% xfc(3) = xfc(2) + mean_dist_fol;
% xfc(4) = 170.75+84.56;
% xfc(5) = xfc(4) + mean_dist_fol;
% xfc(6) = xfc(5) + mean_dist_fol;
% xfc(7) = xfc(1);
% xfc(8) = xfc(2);
% xfc(9) = xfc(3);
% 
% 
% yfc(1) = 170.95;
% yfc(2) = 170.95;
% yfc(3) = 170.95;
% yfc(4) = 317.08;
% yfc(5) = 317.08;
% yfc(6) = 317.08;
% yfc(7) = 317.08+146.15;
% yfc(8) = 317.08+146.15;
% yfc(9) = 317.08+146.15;


% centre_fol_x(1) = 158.5;
% centre_fol_x(2) = 158.5;
% centre_fol_x(3) = 475.5;
% centre_fol_x(4) = 475.5;
% 
% centre_fol_y(1) = 158.5;
% centre_fol_y(2) = 475.5;
% centre_fol_y(3) = 475.5;
% centre_fol_y(4) = 158.5;



% follicles=[1:7,9:num_fol-1];
%plotting hair follicles:
%for each of the num_fol follicles,
for i=1:num_fol
%     i=follicles(j);
    %generate a vector time points to evaluate the circle at. These points
    %should have range [centre_fol_x - radius_fol, centre-fol_x +
    %radius_fol
    xpts(:,i) = linspace(xfc(i)-radius_fol,xfc(i)+radius_fol,num_pts);
    
    %evaluate the y coordiantes of the follicle at xpts and plot the follicle
    circ_plus(:,i) = fplus(xpts(:,i)',xfc(i),yfc(i),radius_fol);
%     plot(xpts,circ_plus,'-b');
%     hold on

    circ_minus(:,i) = fminus(xpts(:,i)',xfc(i),yfc(i),radius_fol);
%     plot(xpts,circ_minus,'-b');
%     hold on
end
% % 
%set axes limits and aspect ratio of the figure.
% xlim([0 Lx]);
% ylim([0 Ly]);
% pbaspect([1 1 1]);

%initialise vectors to hold the initial and total number of removed agents
%for each repeat
num_removed_init=zeros(1,M);
num_removed_total=zeros(1,M);


%% the crux of the algorithm

% to write video

%%to write video
if video==1
    vid_name=sprintf('hair_follicle_no_persistence_test.avi');
    v = VideoWriter(vid_name);
    v.FrameRate = 20;
    v.Quality = 100;
    open(v);
    set(gcf,'position',[0,0,1080,1920]);
end

prop_div_cells=zeros(1,M);
prop_nondiv_cells=zeros(1,M);

num_prolif_rem_vec=zeros(1,M);
num_nondiv_rem_vec=zeros(1,M);
sum_outin_vec=zeros(1,M);
sum_inin_vec=zeros(1,M);
num_prolif_tot_vec=zeros(1,M);


%for each of M repeats
parfor m=1:M
    
    rng('shuffle');
    
    % parameters to be reinitialised before each repeat
    dist_mat = []; %preallocate matrix to hold distance of each agent from nearest follicle
    
    %uniformly distributed random times for prolifeataion events to take place.
    %the (t-10) ensures the last prolif event ends before tfin
    tprol = rand(num_prolif_total,1)*(tfin-10);
%     tprol=10:12:tfin;
%     tprol=tprol';
%     tprol=linspace(10,tfin-10,num_prolif_total)';
    tprol = sort(tprol); %tprol sorted in ascending order
    tprol = [tprol; inf]; %adding inf to the tprol to calibrate the size of the vector
    tprol_current = tprol(1); %current proliferaiton time: this is the time at
    %next proliferation event would take place.
    
    proliferating = zeros(N_init,1); %cells which are currently proliferating
    proliferated_tot = zeros(N_init,1); %cells which have prolifearted upto time t
    prolif_agent_index = []; %vector holding the indices of currently proliferating agents
    angle_mitosis_vec = []; %vector to hold angle of mitosis of currently proliferating agents
    time_prolif_stop = []; %vector to hold time when agents come out of their proliferating phase
    num_prolif_current=0; %preallocate the number of currently prolferating agents
    num_prolif_vec = []; %to hold the number of proliferating agents with time
    
    %initalise the recording index, rec_index = 1 means the first frame
    %recorded will be rec_times(1)
    rec_index = 1;
    
    
    %preallocate the number of removed cells (cells entering follicles)
    sum_removed = zeros(1,num_iter);
    
    %initial cell positions
    X_init = ax*ones(N_init,1) + Lx*rand(N_init,1);
    Y_init = ay*ones(N_init,1) + Ly*rand(N_init,1);
    
%     X_init=600;
%     Y_init=600;
    
%     X_init=319.3813;
%     Y_init=253.1323;
    
    
    
    N = N_init;
    
    %Don't need to keep track of all the positions at every step
    X = X_init;
    Y = Y_init;
    
    if num_fol>0
        %determine distance of each cell from each follicle
        for k=1:num_fol
            dist_mat(:,k) = sqrt((X-xfc(k)).^2+(Y-yfc(k)).^2);
        end

        %Find the minimum distance to any follicle
        min_dist_from_fol = min(dist_mat,[],2);

        %Mark all cells within a given distance as removed
        removed_init = min_dist_from_fol <= radius_fol;

        %initially removed cells
        num_removed_init(m)=sum(removed_init);
        
    else
        removed_init=zeros(N_init,1);
    end
    
    
    %counter for number of proliferation events (must be > 0 for indexing
    %proliferation and stoppping times
    i = 1;
    
    X_old = [];
    Y_old = [];
    X_mitosis=[];
    Y_mitosis=[];
    
    t=0; %initalise current time to 0
    
    
    
    %6hours prior copies of X and Y
%     Xprev=X;
%     Yprev=Y;

    removed=removed_init;
    removed_mat=removed;
    
    %select cells which are non already inside a follicle to be sampled as
    %on dividing cells
    non_div_cells_out=(1:N)'.*(1-removed_init);
    
    %delete zero values
    non_div_cells_out(non_div_cells_out==0)=[];
    
    %sample non dividng cells from the cells available from above
    non_div_sample=datasample(non_div_cells_out,num_non_div_cells,'Replace',false);
    
    %randomly sample the same num of non dividng cells as dividng cells to
    %tracks
    non_div_vec=zeros(N,1);
    
    %create a logical vector where a 1 denotes a non-dividng cell. We use
    %this vector later as we do not want these agents to divide
    non_div_vec(non_div_sample)=1;
    
    
    for j = 1:num_iter
        
        if video==1
            %for the movie -- ignore
            if t == rec_times(rec_index)
                clf('reset');
                if num_non_div_cells>0
                    plot(X(non_div_vec==1),Y(non_div_vec==1),'o','Color',[0, 166/255, 81/255],'MarkerFaceColor',[0, 166/255, 81/255]);
                    hold on
                end

                plot(X(non_div_vec==0 & removed_init==0),Y(non_div_vec==0 & removed_init==0),'or','MarkerFaceColor','r');
                hold on

                plot(X(removed_init==1),Y(removed_init==1),'oc','MarkerFaceColor','c');
                hold on

                for p = 1:num_fol
                    plot(xpts(:,p),circ_plus(:,p),'-k');
                    hold on
                    plot(xpts(:,p),circ_minus(:,p),'-k');
                    hold on
                end


                hold off

                time_in_hours = t/60;
                title_ = sprintf('T = %f',t);
                title(title_);
                xlim([ax, Lx]);
                ylim([ay, Ly]);
                axis square

                frame = getframe(gcf);
                writeVideo(v,frame);
                rec_index = rec_index + 1;
    % %             pause;
            end
        end
        
        %current time
        t=j*dt;
        
        
        %update cell position
        X=X+sqrt(2*D*dt)*randn(N,1).*(1-removed).*(1-proliferating);
        Y=Y+sqrt(2*D*dt)*randn(N,1).*(1-removed).*(1-proliferating);
        
        
        
        %implementing periodic boundaries
        
        %if an agents new position X is outside of the left boundary then,
        %the agent's new position is given by Lx+X
        X(X < 0) = Lx + X(X < 0);
        
        %otherwise, if an agents new position X is outside of the left
        %boundary then,
        
        
        %the agent's new position is given by X-Ly
        X(X > Lx) = X(X > Lx)-Lx;
        
        %if an agents new position Y is outside of the bottom boundary (y = 0) then,
        
        %the agent's new position is given by Ly+Y
        Y(Y < 0) = Ly + Y(Y < 0);
        
        %otherwise, if an agents new position Y is outiside of the top
        %boundary L_y then,
        
        %the agent's new position is given by Y-Ly
        Y(Y > Ly) = Y(Y > Ly)-Ly;
        
        
%         %update position based on the new position of a cell being outside of the domain or inside. We
%         %reflect the position is the new cell poisiton is outside of domain
%         X = (ax-X).*(X<ax)+(2*Lx-X).*(X>Lx)+X.*(X>ax & X<Lx);
%         Y = (ay-Y).*(Y<ay)+(2*Ly-Y).*(Y>Ly)+Y.*(Y>ay & Y<Ly);
        
        num_prolif_vec = [num_prolif_vec num_prolif_current];
        
        %     Xtrack(:,j) = X;
        %     Ytrack(:,j) = Y;
        
        if sum(proliferating)~=N
            if t >= tprol_current
                
                %agents which are not in an follicle (removed) and
                %are not proliferating already
                ind_agents_avail_prolif = (1:N)'.*(1-removed).*(1-proliferating).*(1-non_div_vec);
                ind_agents_avail_prolif(ind_agents_avail_prolif==0) = [];
                
                %sample an agent from the pool of non-removed and non-proliferating
                %agents
                new_prolif_agent_index = datasample(ind_agents_avail_prolif,1);
                
                
                %position coordinates of mitosis event (the old position of the
                %dividing cells)
                X_mitosis = [X_mitosis X(new_prolif_agent_index)];
                Y_mitosis = [Y_mitosis Y(new_prolif_agent_index)];
                
                %add new entry to X and Y for the position of new cell
                X=[X;X(new_prolif_agent_index)];
                Y=[Y;Y(new_prolif_agent_index)];
                
                X_init=[X_init;X(end)];
                Y_init=[Y_init;Y(end)];
                
                %agents currently proliferating
                proliferating(new_prolif_agent_index) = 1;
                proliferating = [proliferating;1];
                
                %all the agents proliferated so far
                proliferated_tot(new_prolif_agent_index) = 1;
                proliferated_tot = [proliferated_tot;1];
                
                %sum(proliferated_tot) might not always add up to
                %num_prolif_tot as some cells may divide more than once.
                
                
                %sample the direction in which the new cells will be propelled
%                 angle_mitosis = randn*0.5367+1.638;%1.696;
                angle_mitosis = circ_vmrnd(mu_mit,kappa_,mit,1);
                
                %         angle_mitosis = pi*rand;
                
                %index of the currently proliferating agents
                prolif_agent_index=[prolif_agent_index,new_prolif_agent_index,N+1];
                
                %angle of mitosis of the currently proliferating agents
                angle_mitosis_vec=[angle_mitosis_vec,angle_mitosis, angle_mitosis+pi];
                
                %increase the number of currently proliferatin agents by 2 (since two
                %agents are unger groing prolifeartion at a time (the parent and its
                %progeny)
                num_prolif_current=num_prolif_current+2;

                
                % Time to stop proliferating for each pair of agents (parent
                % and daughter)
                time_prolif_stop=[time_prolif_stop,t+10,t+10];
                
                %increment number of cells
                N=N+1;
                
                %increase size of non_div_vec to account for the new cell
                non_div_vec=[non_div_vec;0];
                
                %increase size of removed_init to account for the new cell
                removed_init=[removed_init;0];
                
                removed=[removed;0];
                    
                %increment the number of prolferation events happened til time
                %t
                i=i+1;
                tprol_current = tprol(i);
            end
        end
        
        %if there are currently proliferating agents then,
        if num_prolif_current > 0
            
            %if curent time >= min(time_prolif_stop) (equivalent to
            %time_prolif_stop(1) as time_prolif_stop is an non-decreasing
            %sequence
            if t>time_prolif_stop(1)
                X_old=[X_old; [X(prolif_agent_index(1)) X(prolif_agent_index(2))]];
                Y_old=[Y_old; [Y(prolif_agent_index(1)) Y(prolif_agent_index(2))]];
                
                %discount for the two agents which have just finished
                %proliferating
                num_prolif_current=num_prolif_current-2;
                
                %reset the corresponding entries in proliferating to 0
                proliferating(prolif_agent_index(1))=0;
                proliferating(prolif_agent_index(2))=0;
                
                %reset time_prolif_stop to nothing
                time_prolif_stop(1:2)=[];
                
                %reset prolif_agent_index to nothing
                prolif_agent_index(1:2)=[];
                
                %reset angle_mitosis_vec to nothing
                angle_mitosis_vec(1:2)=[];
                
            end
            
            %move each proliferating agent a distance jump_mitosis
            for kk=1:num_prolif_current
                
                %find new position of kkth cell
                xi_x = X(prolif_agent_index(kk)) + dt*jump_mitosis*cos(angle_mitosis_vec(kk));
                xi_y = Y(prolif_agent_index(kk)) + dt*jump_mitosis*sin(angle_mitosis_vec(kk));
                
                
                %account for periodic boundaries
                
                if xi_x < 0
%                     xi_x=ax-xi_x;
                    
                    xi_x = Lx + xi_x;
                    
                elseif xi_x > Lx
                    
%                     xi_x=2*Lx-xi_x;
                    xi_x = xi_x-Lx;
                
                
                elseif xi_y < 0
                    
%                     xi_y=ay-xi_y;
                    xi_y = Ly + xi_y;
                    
                elseif xi_y > Ly
%                     xi_y=2*Ly-xi_y;
                    xi_y= xi_y-Ly;
                    
                end
                
                %update the X and Y vectors to the new position
                X(prolif_agent_index(kk)) = xi_x;
                Y(prolif_agent_index(kk)) = xi_y;
            end
        end
        
        if num_fol>0
        %distance matrix
            dist_mat = [];
            for k=1:num_fol
                dist_mat = [dist_mat, sqrt((X-xfc(k)).^2+(Y-yfc(k)).^2)];
            end




            %Find the minimum distance to any follicle
            min_dist_from_fol = min(dist_mat,[],2);

            %Mark all cells within a given distance as removed
            removed = min_dist_from_fol <= radius_fol;

            %Count how many cells are now removed
            sum_removed(j) = sum(removed);
        end
        
    end
    
    num_removed_total(m)=sum(removed);
    
    inin=[removed_init];
    outin=removed-inin;
    prolif_and_removed=proliferated_tot.*outin;
    
    non_div_and_removed=non_div_vec.*outin;
    
    num_prolif_rem_vec(m)=sum(prolif_and_removed);
    num_nondiv_rem_vec(m)=sum(non_div_and_removed);
    sum_outin_vec(m)=sum(outin);
    sum_inin_vec(m)=sum(inin);
    
    sum_prolif_tot=sum(proliferated_tot);
    
    num_prolif_tot_vec(m)=sum_prolif_tot;
    
    prop_div_cells(m)=sum(prolif_and_removed)/num_prolif_total;

%     noprolif_and_removed=(1-proliferated_tot).*outin;
%     prop_nondiv_cells(m)=sum(noprolif_and_removed)/(N-num_prolif_total);
    
    prop_nondiv_cells(m)=sum(non_div_and_removed)/num_prolif_total;
    
    if video==1
        close(v)
    end
end
%%
% sum_removed=sum(removed);
% sum_inin=sum(removed_init);
% sum_outin=sum(outin);
% sum_prolif_removed=sum(prolif_and_removed);
% sum_noprolif_removed=sum(non_div_and_removed);
% 
% fprintf("sum_removed=%d\nsum_inin=%d\nsum_outin=%d\nsum_prolif_removed=%d\nsum_noprolif_removed=%d\n",sum_removed,sum_inin,...
%     sum_outin,sum_prolif_removed,sum_noprolif_removed);
% 
% save(vid_name+".mat");
%% find distance from cell to condensate

% dist_to_condensates_mat = zeros(num_prolif_total,num_fol);
% 
% for i=1:num_fol
%     dist_to_condensates_mat(:,i) = sqrt((X_mitosis-centre_fol_x(i)).^2+(Y_mitosis-centre_fol_y(i)).^2);
% end


% num_removed_diff=abs(num_removed_total-num_removed_init);
% mean_removed_total=mean(num_removed_total);


%%
% mean_removed_diff=mean(num_removed_diff);


%% plot position of mitoses

% figure;
% plot(X_mitosis,Y_mitosis,'.r');
% xlabel('$x\;\mu$m','interpreter','latex');
% ylabel('$y\;\mu$m','interpreter','latex');
% title('Locations of mitosic events','interpreter','latex');
% xlim([0,Lx]);
% ylim([0,Ly]);
% pbaspect([1 1 1]);

%% cells which proliferated and eventually ended up in a follicle

% diff_angle=abs(angle2fol-angle_of_agent);
% diff_angle(isnan(diff_angle))=0;

% for i=1:length(rec_times)    
%     x_cond=eucl_dist(:,i).*cos(diff_angle(:,i)).*(1-removed_mat(:,i));
%     x_cond_copy=x_cond;
%     y_cond=eucl_dist(:,i).*sin(diff_angle(:,i)).*(1-removed_mat(:,i));
%     x_cond(x_cond>400)=[];
%     % x(isnan(x))=[];
%     y_cond(x_cond_copy>400)=[];
%     % y(isnan(xcopy))=[];
% 
%     compass(x_cond,y_cond);
%     hold on
% end

% X_removed_and_prolif = X.*removed.*proliferated_tot;
% Y_removed_and_prolif = Y.*removed.*proliferated_tot;
% 
% %prob of proliferating and ending up in a follicle (possibly no correct,
% %check with Kit)
% prob = nnz(X_removed_and_prolif)/nnz(removed);

%% plot direction of mitoses

% len_X_old= size(X_old,1);
% dist_mitosis=sqrt((X_old(:,2)-X_old(:,1)).^2+(Y_old(:,2)-Y_old(:,1)).^2);
% 
% figure;
% for i = 1:len_X_old
%     if sqrt((X_old(i,2)-X_old(i,1))^2+(Y_old(i,2)-Y_old(i,1))^2)<=2*min(dist_mitosis)+eps
%         plot([X_old(i,2),X_old(i,1)],[Y_old(i,2),Y_old(i,1)],'-b');
%         hold on
%     end
% end
% 
% xlabel('$x\;\mu$m','interpreter','latex');
% ylabel('$y\;\mu$m','interpreter','latex');
% title('Angle of mitosis and mitotic jumps','interpreter','latex');
% xlim([0,Lx]);
% ylim([0,Ly]);
% pbaspect([1 1 1]);

%% angle of deviation of a follicle




% 
% for i=1:N
%     plot([dist2folx(i,3)],[dist2foly(i,3)],'o');
%     hold on
% end

%%

% abs_diff_angle=abs(diff_angle);
% 
% condensate_cells=abs_diff_angle.*removed;
% condensate_cells(condensate_cells==0)=[];
% condensate_cells(isnan(condensate_cells))=[];
% 
% 
% intercond_cells=abs_diff_angle.*(1-removed);
% intercond_cells(intercond_cells==0)=[];
% intercond_cells(isnan(intercond_cells))=[];
% 
% figure;
% polarhistogram(condensate_cells,100);
% figure;
% polarhistogram(intercond_cells,100);
% 
% legend('condensate cells','intercondensate cells');

%% out out and out in analysis

mean_prop_div_cells=mean(prop_div_cells);
mean_prop_nondiv_cells=mean(prop_nondiv_cells);

var_prop_div_cells=var(prop_div_cells);
var_prop_nondiv_cells=var(prop_nondiv_cells);

%plotting

% figure;
% 
% b=bar([prop_div_cells', prop_nondiv_cells']);
% 
% b(2).FaceColor='g';
% b(2).EdgeColor='none';
% b(1).EdgeColor='none';
% 
% 
% hold on
% plot(mean_prop_div_cells*ones(1,M),'r-','LineWidth',2)
% hold on
% plot(mean_prop_nondiv_cells*ones(1,M),'k-','LineWidth',2)
% legend('prop div','prop non-div','mean div','mean non-div')

%% bar chart for the mean

% figure;
% b2=bar(18.45);
% hold on
% b3=bar(2, 12.75);
% 
% b2.EdgeColor='none';
% b3.EdgeColor='none';
% b3.FaceColor='r';
% 
% names = {'Dividing cells';'Non-dividing cells'};
% 
% set(gca,'xtick',[1:2],'xticklabel',names)
% 
% ylabel('Proportion %');

if Save==1
    save("dermal_condenate_v9_no_persistence.mat",'-v7.3');
end


% inin=[removed_init; zeros(num_prolif_total,1)];
% outin=removed-inin;
% prolif_and_removed=proliferated_tot.*outin;
% 
% prop_div_cells=sum(prolif_and_removed)/num_prolif_total;
% 
% noprolif_and_removed=(1-proliferated_tot).*outin;
% prop_nondiv_cells=sum(noprolif_and_removed)/(N-num_prolif_total);

%% scrap code

% Plot at t = 0
% subplot(1,2,1);
% 
% %plot the cells
% scatter(x0,y0,10,'ok','Filled');
% hold on
% 
% npts = 100;
% 
% %plot follicles
% for i = 1:nf
%     xpts = linspace(xfc(i)-rf,xfc(i)+rf,npts);
%     ypts = linspace(yfc(i)-rf,yfc(i)+rf,npts);
%     
%     circ_plus = fplus(xpts,ypts,xfc(i),yfc(i),rf,npts);
%     plot(xpts,circ_plus,'-b');
%     hold on
%     
%     circ_minus = fminus(xpts,ypts,xfc(i),yfc(i),rf,npts);
%     plot(xpts,circ_minus,'-b');
%     hold on
% end
% 
% %other plotting parameters
% xlim([ax, Lx]);
% ylim([ay, Ly]);
% pbaspect([1 1 1]);
% 
% title_init = sprintf("t = %d, %d agents, D = %f",0,n0,D);
% title(title_init)
% 
% 
% % plot at t = tfin
% subplot(1,2,2);
% 
% Xend = X;
% Yend = Y;
% 
% %plot cells
% scatter(Xend,Yend,10,'ok','Filled');
% hold on
% 
% 
% %plot follicles
% for i = 1:nf
%     xpts = linspace(xfc(i)-rf,xfc(i)+rf,npts);
%     ypts = linspace(yfc(i)-rf,yfc(i)+rf,npts);
%     
%     circ_plus = fplus(xpts,ypts,xfc(i),yfc(i),rf,npts);
%     plot(xpts,circ_plus,'-b');
%     hold on
%     
%     circ_minus = fminus(xpts,ypts,xfc(i),yfc(i),rf,npts);
%     plot(xpts,circ_minus,'-b');
%     hold on
% end
% 
% 
% %other plotting parameters
% xlim([ax, Lx]);
% ylim([ay, Ly]);
% pbaspect([1 1 1]);
% title_end = sprintf("t = %d, %d agents, D = %f",tfin,n0,D);
% title(title_end)
% 
% %% cell tracking
% 
% dx = zeros(n0,jmax-1);
% dy = zeros(n0,jmax-1);
% 
% 
% for j = 2:jmax
%     dx(:,j-1) = Xtrack(:,j)-Xtrack(:,j-1);
%     dy(:,j-1) = Ytrack(:,j)-Ytrack(:,j-1);
% end
% 
% 
% 
% dz = sqrt(dx.^2 + dy.^2);
% 

