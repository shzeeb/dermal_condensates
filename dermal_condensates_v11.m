%Name: dermal_condensates_v11
%Description: Cell proliferation occurs but the parent cell does not
%produce a daughter cell. The parent cell perform the solo mitotic jump
%followed by persistent motion for 180 mins (3 hrs).


clear;

% parameters not required to reallocate before each repeat
% initalise simulation parameters here:
%domain size
ax = 0;
Lx = 634.15;
ay = 0;
Ly = 634.15;
D = 0.078036;%0.0928;%0.1504;%0.15836;%0.2643;%0.413304;%0.2196;%0.5465;%0.01207175; % diffusion coefficient. WT D = 0.657
dt = 0.1; %dt needs to be sufficiently small such that the position jumps do not take a partile outside of the domain
tfin = 1800; %
M=1; %number of repeats
num_iter = ceil(tfin/dt);
N_init=1936; %initial number of agents on the domain. WT n0 = 1936
num_fol=14; %number of hair follicles
r_0 = 35; %radius of hair follicle (fixed for all follicles)
d_sep = 170; %mean distance between the centre of one follicle to the centre of an adjacent one
jump_mitosis = 0.64271;%0.6325; %jump size per agent at mitosis. Units: microns per minute
num_prolif_total=623; %averaged over all videos
num_non_div_cells=num_prolif_total;
mu_mit=1.562;
kappa_mit=2.85;
mu_pers=-0.0782;
kappa_pers=1.07;
ghost=0; %indicator variable to introduce ghost nodes (1) or not (0). Always set to 0 in this simulatuion code.
persist=1;

rec_step=10;

%times to record state of system at regular intervals
rec_times = 0:rec_step:tfin;

% rec_times=[rec_times,inf];

% rec_times=linspace(1,1800,4);
% rec_times=[600 1200 1800];

%number of points to plot the follcles (npts only used to plot follicles).
num_pts = 1000;

% x points to evaluate fplus and fminus at.
xpts = zeros(num_pts,1);

%vectors to e populated by the y cooridinates of each hair follicle.
circ_plus = zeros(num_pts,1);
circ_minus = zeros(num_pts,1);

fplus = @(xpts,xfc,yfc,r_0) sqrt(r_0^2 - (xpts-xfc).^2) + yfc;
fminus = @(xpts,xfc,yfc,r_0) -sqrt(r_0^2 - (xpts-xfc).^2) + yfc;


[xfc,yfc]=follicles(Lx,Ly,num_fol,r_0,d_sep,ghost);


% list_fol=[1:7,9:15];

%plotting hair follicles:
%for each of the num_fol follicles,

% for i = 1:num_fol
% %     i=list_fol(j);
%     %generate a vector time points to evaluate the circle at. These points
%     %should have range [xfc - r_0, xfc +
%     %r_0
%     xpts = linspace(xfc(i)-r_0,xfc(i)+r_0,num_pts);
%     
%     %evaluate the y coordiantes of the follicle at xpts and plot the follicle
%     circ_plus=fplus(xpts,xfc(i),yfc(i),r_0);
%     plot(xpts,real(circ_plus),'-b');
%     hold on
% 
%     circ_minus=fminus(xpts,xfc(i),yfc(i),r_0);
%     plot(xpts,real(circ_minus),'-b');
%     hold on
% end
% 
% % set axes limits and aspect ratio of the figure.
% xlim([0 Lx]);
% ylim([0 Ly]);
% pbaspect([1 1 1]);

%initialise vectors to hold the initial and total number of removed agents
%for each repeat
num_removed_init=zeros(1,M);
num_removed_total=zeros(1,M);


%%the crux of the algorithm

% to write video

vid_name="test_case_3_persistence";
set(gca,'nextplot','replacechildren');
v = VideoWriter(vid_name+".avi");
open(v);

%%recording variables. Recoded once at the end of a each repeat

%proportion of dividing cells
prop_div_cells=zeros(1,M);

%proportion of non-dividing cells
prop_nondiv_cells=zeros(1,M);

%number of dividing cells which eventually got adsorbed
num_prolif_rem_vec=zeros(1,M);

%%number of non-dividing cells which eventually got adsorbed
num_nondiv_rem_vec=zeros(1,M);

%number of out-in cells (cells which start outside and end up inside
%a follicle after some arbitrary time)
sum_outin_vec=zeros(1,M);

%number of in-in cells (initially removed cells)
sum_inin_vec=zeros(1,M);

%number of proliferations
num_prolif_tot_vec=zeros(1,M);

%division positions x coordinate
X_mitosis_mat=zeros(num_prolif_total,M);

%division positions y coordinate
Y_mitosis_mat=zeros(num_prolif_total,M);

%initial position of all cells x coordiantes
X_init_mat=zeros(N_init,M);

%initial position of all cells y coordiantes
Y_init_mat=zeros(N_init,M);

%removed cells per repeat. Indicator matrix. 0=not removed. 1=removed
removed_mat=zeros(N_init,M);

%X position of each cells at different recording times
X_mat=zeros(N_init,M);

%Y position of each cells at different recording times
Y_mat=zeros(N_init,M);

%dividing cells per repeat. Indicator matrix. 0= cells didnt divide. 1=cell
%divided once, >1 = cell divided more than once
prolif_mat=zeros(N_init,M);

%initially removed cells. Indicator matrix. 0=not removed. 1=removed
removed_init_mat=zeros(N_init,M);

%time of adsorption of each removed cell
t_removed_mat=zeros(N_init,M);


%for each of M repeats
for m=1:M
    % parameters to be reinitialised before each repeat
    dist_mat = []; %preallocate matrix to hold distance of each agent from nearest follicle
    
    %uniformly distributed random times for prolifeataion events to take place.
    %the (t-10) ensures the last prolif event ends before tfin
    tprol = rand(num_prolif_total,1)*(tfin-10);
%     tprol=[10;20;25];
%     tprol=10:12:tfin;
%     tprol=tprol';
%     tprol=linspace(10,tfin-10,num_prolif_total)';
    tprol = sort(tprol); %tprol sorted in ascending order
%     tprol=10;
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
    
    %time of adsorption of each removed cell for a single repeat
    t_removed=zeros(N_init,1);


    %initalise the recording index, rec_index = 1 means the first frame
    %recorded will be rec_times(1)
    rec_index = 2;
    
    
    
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
    dist2folx=[];
    dist2foly=[];
    
    %determine distance of each cell from each follicle
    for k=1:num_fol
        dist_mat(:,k) = sqrt((X-xfc(k)).^2+(Y-yfc(k)).^2);
    end
    
    %Find the minimum distance to any follicle
    min_dist_from_fol = min(dist_mat,[],2);
    
    %Mark all cells within a given distance as removed
    removed_init = min_dist_from_fol <= r_0;
    
    %initially removed cells
    num_removed_init(m)=sum(removed_init);
    
    
    %counter for number of proliferation events (must be > 0 for indexing
    %proliferation and stoppping times
    i = 1;
    
    X_old = [];
    Y_old = [];
    X_mitosis=[];
    Y_mitosis=[];
    
    t=0; %initalise current time to 0
    
    removed=removed_init;
    angle2fol=[];
    angle_of_agent=[];
    diff_angle=[];
    diffX=[];
    diffY=[];
    dist2folx_mat=[];%zeros(N,length(rec_times));
    dist2foly_mat=[];%zeros(N,length(rec_times));
    diffX_mat=[];
    diffY_mat=[];
    eucl_dist=[];
    removed_mat=[];
    
    %6hours prior copies of X and Y
    Xprev=X;
    Yprev=Y;
    removed_mat=removed_init;
    
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
    
    %persistence vector. column 1 is either 0,1or 2. A 1 is for daughter A
    %and a 2 is for daughter B. Could 2 contains time when the cells start
    %moving persistently
    persistence=zeros(N,2);
    
    %number of currently persisting agents
    num_persist_current=0;
    
    persistent_agents=[];
    
    iter_pers=[];
    
    angle_pers=zeros(N,1);
    
    %vector to contain differnt stages for persistence (18 stages)
%     persist_stage=zeros(N,2);
    
    %angle pers old
    angle_pers_old=0;
    
    %set of times when pairs of cells should stop persistent motion
    t_pers_stop=[];
    
    %indices of persistent agents
    pers_agent_index=[];
    
    %mitosis angles of persistent cells
    angle_mit_pers=zeros(N,1);
%     angle_mit_pers_b=zeros(N,1);
    
    X_mat=[];
    Y_mat=[];
    
    for j = 1:num_iter
        
        %current time
        t=j*dt;
        
        %distace vector
        
        
        
        %update cell position
        X=X+sqrt(2*D*dt)*randn(N,1).*(1-removed).*(1-proliferating).*(persistence(:,1)==0);
        Y=Y+sqrt(2*D*dt)*randn(N,1).*(1-removed).*(1-proliferating).*(persistence(:,1)==0);
        
        
        
        %implementing periodic boundaries
        
        %if an agents new position X is outside of the left boundary then,
        %the agent's new position is given by Lx+X
        X(X < 0) = Lx + X(X < 0);
        
        %otherwise, if an agents new position X is outside of the left
        %boundary then,
        
        
        %the agent's new position is given by X-Lx
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
                ind_agents_avail_prolif = (1:N)'.*(1-removed).*(1-proliferating).*(persistence(:,1)==0).*(1-non_div_vec);%.*(proliferated_tot~=0); %uncomment (1-proliferated_tot) if dont want cells to divide twice
                ind_agents_avail_prolif(ind_agents_avail_prolif==0) = [];
                
                %sample an agent from the pool of non-removed and non-proliferating
                %agents
                
                if ~isempty(ind_agents_avail_prolif)
                    new_prolif_agent_index = datasample(ind_agents_avail_prolif,1);
                end
                
                %position coordinates of mitosis event (the old position of the
                %dividing cells)
                X_mitosis = [X_mitosis X(new_prolif_agent_index)];
                Y_mitosis = [Y_mitosis Y(new_prolif_agent_index)];
                
                %add new entry to X and Y for the position of new cell
%                 X=[X;X(new_prolif_agent_index)];
%                 Y=[Y;Y(new_prolif_agent_index)];
                
%                 X_init=[X_init;X(end)];
%                 Y_init=[Y_init;Y(end)];
                
                %agents currently proliferating
                proliferating(new_prolif_agent_index) = 1;
%                 proliferating = [proliferating;1];
                
                %all the agents proliferated so far
                proliferated_tot(new_prolif_agent_index) = proliferated_tot(new_prolif_agent_index)+1;
%                 proliferated_tot = [proliferated_tot;1];
                
                %sum(proliferated_tot) might not always add up to
                %num_prolif_tot as some cells may divide more than once.
                
                
                %sample the direction in which the new cells will be propelled
%                 angle_mitosis = randn*0.5367+1.4511;%1.696;
                angle_mitosis=circ_vmrnd(mu_mit,kappa_mit,1);
                
                
                %index of the currently proliferating agents
                prolif_agent_index=[prolif_agent_index,new_prolif_agent_index];
                
                %angle of mitosis of the currently proliferating agents
                angle_mitosis_vec=[angle_mitosis_vec,angle_mitosis];
                
                %increase the number of currently proliferatin agents by 2 (since two
                %agents are unger groing prolifeartion at a time (the parent and its
                %progeny)
                num_prolif_current=num_prolif_current+1;

                
                % Time to stop proliferating for each pair of agents (parent
                % and daughter)
                time_prolif_stop=[time_prolif_stop,t+10];
                
                %increment number of cells
%                 N=N+1;
                
                %increase size of non_div_vec to account for the new cell
%                 non_div_vec=[non_div_vec;0];
                
                %increase size of persistence matrix to account for the new
                %cell
%                 persistence=[persistence;[0 0]];
                
%                 persist_stage=[persist_stage;[0 0]];
                
                %increase size of angle_mit_pers to account for the new
                %cell
%                 angle_mit_pers=[angle_mit_pers;0];
                
%                 angle_pers=[angle_pers;0];
                
%                 removed=[removed;0];
                
                %increment the number of prolferation events happened til time
                %t
                i=i+1;

                %select the currnet proliferation time
                tprol_current = tprol(i);
            end
        end
        
        %if there are currently proliferating agents then,
        if num_prolif_current > 0
            
            %if curent time >= min(time_prolif_stop) (equivalent to
            %time_prolif_stop(1) as time_prolif_stop is an non-decreasing
            %sequence
            if t>time_prolif_stop(1)
                X_old=[X_old; [X(prolif_agent_index(1))]];
                Y_old=[Y_old; [Y(prolif_agent_index(1))]];
                
                %discount for the two agents which have just finished
                %proliferating
                num_prolif_current=num_prolif_current-1;
                
                %reset the corresponding entries in proliferating to 0
                proliferating(prolif_agent_index(1))=0;
%                 proliferating(prolif_agent_index(2))=0;

                if persist==1
                    %indices of persistent agents
                    pers_agent_index=[pers_agent_index, prolif_agent_index(1)];

                    %set persistence of last proliferated agent
                    persistence(pers_agent_index(end),1)=1;
    %                 persistence(pers_agent_index(end),1)=2;

                    persistence(pers_agent_index(end),2)=t;
    %                 persistence(pers_agent_index(end),2)=t;

                    %persistence stages of cells (num between 1:100)
    %                 persist_stage(pers_agent_index(end-1),1)=1;
    %                 persist_stage(pers_agent_index(end),1)=1;
    %                 
    %                 persist_stage(pers_agent_index(end-1),2)=t+10;
    %                 persist_stage(pers_agent_index(end),2)=t+10;

                    iter_pers=[iter_pers, t+10];

                    %persistent agents vector (used during iteration movement
                    %update
                    persistent_agents=[persistent_agents, prolif_agent_index(1)];

                    %num of currently persisting agents
                    num_persist_current=num_persist_current+1;

                    %old angle of persistence for persistent cells
                    angle_pers(pers_agent_index(end))=circ_vmrnd(mu_pers,kappa_pers,1);

                    %angle mitosis for persistence
                    angle_mit_pers(prolif_agent_index(1))=angle_mitosis_vec(1);


                    %time when persistence motion should stop
                    t_pers_stop=[t_pers_stop, time_prolif_stop(1)+180];
                
                end
                
                %reset time_prolif_stop to nothing
                time_prolif_stop(1)=[];
                
                %reset prolif_agent_index to nothing
                prolif_agent_index(1)=[];
                
                %reset angle_mitosis_vec to nothing
                angle_mitosis_vec(1)=[];
                
            end
            
            
            
            %move each proliferating agent a distance jump_mitosis
            for kk=1:num_prolif_current
                
                %find new position of kkth cell
                xi_x = X(prolif_agent_index(kk)) + dt*jump_mitosis*cos(angle_mitosis_vec(kk))*(1-removed(prolif_agent_index(kk)));
                xi_y = Y(prolif_agent_index(kk)) + dt*jump_mitosis*sin(angle_mitosis_vec(kk))*(1-removed(prolif_agent_index(kk)));
                
                
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
        
        %distance matrix
        dist_mat = [];
        for k=1:num_fol
            dist_mat = [dist_mat, sqrt((X-xfc(k)).^2+(Y-yfc(k)).^2)];
        end
        

        %Find the minimum distance to any follicle
        min_dist_from_fol = min(dist_mat,[],2);
        
        removed_old=removed;
        %Mark all cells within a given distance as removed
        removed = min_dist_from_fol <= r_0;
        removed_new=removed-removed_old;
        
        %record the time of adsorption
        t_removed(logical(removed_new))=t;
   
        
        %Count how many cells are now removed
        sum_removed(j) = sum(removed);
        
        
        if persist==1
            %%implement persistence
            %if the numebr of persisting cells is > 0
            if num_persist_current>0
                d=0.1*dt*1.538379;

    %             angle_pers=vmrand(-0.04914,1.2,[num_persist_current,1]);
    %             for n=1:num_persist_current

                    X(pers_agent_index)=X(pers_agent_index)+d*cos(angle_pers(pers_agent_index)+angle_mit_pers(pers_agent_index)).*(1-removed(pers_agent_index));%.*(persistence(pers_agent_index(n),1)~=0);
                    Y(pers_agent_index)=Y(pers_agent_index)+d*sin(angle_pers(pers_agent_index)+angle_mit_pers(pers_agent_index)).*(1-removed(pers_agent_index));%.*(persistence(pers_agent_index(n),1)~=0);
    %             end



                %chaning orientation of persistent cells
                orient_idx=(t>iter_pers);
                nnz_orient_idx=nnz(orient_idx);

                if nnz_orient_idx>0
    %                 for b=1:nnz_orient_idx
                        iter_pers(orient_idx)=iter_pers(orient_idx)+10;

                        rand_orient=circ_vmrnd(mu_pers,kappa_pers,[nnz_orient_idx,1]);%jkappa from data 0.8053

                        angle_pers(persistent_agents(orient_idx))=rand_orient;

                        cat_vec=[iter_pers;persistent_agents];

    %                     sorted_=sortrows(cat_vec',1,"ascend");
    %                     
    %                     iter_pers=sorted_(:,1)';
    %                     persistent_agents=sorted_(:,2)';

    %                     rand_orient=vmrand(-0.04914,1,[nnz_orient_idx,1]);
    %                     for b=1:nnz_orient_idx
    %                         angle_pers(
    %                     end

    %                 end
                end


    %             if j>=iter_pers(1)
    %                 iter_pers(end+1:end+2)=iter_pers(1:2)+100;
    %                 
    %                 persistent_agents(end+1:end+2)=persistent_agents(1:2);
    %                 
    %                 
    %                 angle_pers(persistent_agents(1:2))=vmrand(-0.04914,1,[2,1]);%var=1.2
    %                 
    %                 iter_pers(1:2)=[];
    %                 persistent_agents(1:2)=[];
    %                 
    %             end

    %             for n=1:length(pers_agent_index)

                %prolivisonally set to normal (may need to change later)
    %             pers_dx=randn(-0.04062,1.77975,[2,1]);
    %             pers_dy=randn(0.01050,2.64978,[2,1]);




    %             if sum(persist_stage(:,1)==1)>=1
    %                 angle_pers(persist_stage(:,1)==1)=vmrand(-0.04914,1.2,[sum(persist_stage(:,1)==1),1]);
    %                 
    %             else
    %                 angle_pers(persist_stage(:,1)<=100 & persist_stage(:,1)>1)=angle_pers_old(persist_stage(:,1)<=100 & persist_stage(:,1)>1);
    %             end
    %             
    %             
    %             angle_pers_old=angle_pers;
    %             
    %             persist_stage(:,1)=persist_stage(:,1)+1;



                %sample deviation angles from the von Mises distribution
    %             angle_pers=vmrand(-0.04914,1.2,[N,1]);

    %             angle_pers_vec=[angle_pers, angle_pers(angle_pers~=0)+pi];

    %             d=gprnd(0.285819,0.754854,0.49,[N,1])/10;
    %             d=zeros(N,1);
    %             d(persistence==1)=gprnd(0.285819,0.754854,0.49,[nnz(persistence==1),1]);

                %displacement persistent cells per min
    %             d=1.538379;
    %             
    %             X=X+d*cos(angle_pers+angle_mit_pers).*(1-removed).*(persistence(:,1)~=0);
    %             Y=Y+d*sin(angle_pers+angle_mit_pers).*(1-removed).*(persistence(:,1)~=0);

    %             Y(persistence==1)=Y(persistence==1)+d.*sin(angle_pers);


                %persistent motion of daughter b
    %             angle_pers=zeros(N,1);

    %             d=zeros(N,1);
    %             
    %             d(persistence==2)=gprnd(0.285819,0.754854,0.49,[nnz(persistence==2),1]);

    %             X=X+0.1*d.*cos(angle_pers+angle_mit_pers).*(1-removed).*(persistence(:,1)==2);
    %             Y=Y+0.1*d.*sin(angle_pers+angle_mit_pers).*(1-removed).*(persistence(:,1)==2);





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

    %             X(persistence==2)=X(persistence==2)+d.*cos(angle_pers+pi);
    %             Y(persistence==2)=Y(persistence==2)+d.*sin(angle_pers+pi);

            %stopping persistent motion

                if t>=t_pers_stop(1)
                    persistence(pers_agent_index(1),1)=0;
                    persistence(pers_agent_index(1),2)=0;

                    angle_mit_pers(pers_agent_index(1))=0;

                    angle_pers(pers_agent_index(1))=0;

                    iter_pers(logical((cat_vec(2,:)==pers_agent_index(1))))=[];
                    persistent_agents(logical((cat_vec(2,:)==pers_agent_index(1))))=[];


    %                 iter_pers(end-1:end)=[];
    %                 
    %                 persistent_agents(end-1:end)=[];

    %                 persist_stage(pers_agent_index(1:2),1)=0;
    %                 persist_stage(pers_agent_index(1:2),1)=0;
    %                 
    %                 persist_stage(pers_agent_index(1:2),2)=0;
    %                 persist_stage(pers_agent_index(1:2),2)=0;


                    pers_agent_index(1)=[];

                    num_persist_current=num_persist_current-1;
                    t_pers_stop(1)=[];
                    
                end
            end
        end

        
%         if sum(persist_stage==100)~=0
%             persist_stage(persist_stage==100)=1;
%         end
            
        
        %for the movie -- ignore
        if t == rec_times(rec_index)
            
%             X_mat=[X_mat;zeros(N-size(X_mat,1),size(X_mat,2))];
%             Y_mat=[Y_mat;zeros(N-size(Y_mat,1),size(Y_mat,2))];
%             
%             X_mat=[X_mat,X];
%             Y_mat=[Y_mat,Y];
            
            plot(X,Y,'.k');
            hold on
            
%             plot(X,Y,'.k'); hold on
            
            for p = 1:num_fol
                plot(xpts(:,p),circ_plus(:,p),'-b');
                hold on
                plot(xpts(:,p),circ_minus(:,p),'-b');
                hold on
            end
            
            
            hold off
            
            time_in_hours = t/60;
            title_ = sprintf('T = %f',t);
            title(title_);
            xlim([ax, Lx]);
            ylim([ay, Ly]);
            pbaspect([1 1 1]);
            
            frame = getframe(gcf);
            writeVideo(v,frame);
              pause;
            rec_index = rec_index + 1;
            
        end
        
    end
    
    X_init_mat(:,m)=X_init;
    Y_init_mat(:,m)=Y_init;
    
    removed_mat(:,m)=removed;
    prolif_mat(:,m)=proliferated_tot;
    
    removed_init_mat(:,m)=removed_init;
    
    num_removed_total(m)=sum(removed);
    
    inin=removed_init;
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
    
    X_mitosis_mat(:,m)=X_mitosis;
    Y_mitosis_mat(:,m)=Y_mitosis;
    
    t_removed_mat(:,m)=t_removed';
    
    
    close(v)
end
%% dividing and removed vs non-divding and removed statistics

% sum_removed=sum(removed);
% sum_inin=sum(removed_init);
% sum_outin=sum(outin);
% sum_prolif_removed=sum(prolif_and_removed);
% sum_noprolif_removed=sum(non_div_and_removed);
% 
% fprintf("sum_removed=%d\nsum_inin=%d\nsum_outin=%d\nsum_prolif_removed=%d\nsum_noprolif_removed=%d\n",sum_removed,sum_inin,...
%     sum_outin,sum_prolif_removed,sum_noprolif_removed)
% 
% save(vid_name+".mat");

%% out out and out in analysis

% prop_div_cells=prop_div_cells.*num_prolif_tot_vec/num_prolif_total;

mean_prop_div_cells=mean(prop_div_cells);
mean_prop_nondiv_cells=mean(prop_nondiv_cells);

var_prop_div_cells=var(prop_div_cells);
var_prop_nondiv_cells=var(prop_nondiv_cells);

abs_diff_var=(var_prop_div_cells-var_prop_nondiv_cells);

rel_diff_var=abs_diff_var/var_prop_nondiv_cells*100;

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
% 
% legend('prop div','prop non-div','mean div','mean non-div')
% xlabel('proportion','interpreter','latex');
% ylabel('repeats','interpreter','latex');
% title(sprintf('mean div = %2.3f, mean non-div = %2.3f, variance div=%f, variance non-div=%f',mean_prop_div_cells,mean_prop_nondiv_cells,var_prop_div_cells,var_prop_nondiv_cells),'interpreter','latex');

%% bar chart for the mean

% figure;
% b2=bar(mean_prop_div_cells);
% hold on
% b3=bar(2, mean_prop_nondiv_cells);
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

save('dermal_condensates_no_progeny.mat');
%%

% size_Xmat=size(X_mat);
% 
% for i=1:size_Xmat(1)
%     for j=1:size_Xmat(2)-1
%         if X_mat(i,j)~=0 && X_mat(i,j+1)~=0
%             plot([X_mat(i,j) X_mat(i,j+1)],[Y_mat(i,j) Y_mat(i,j+1)],'-k');
%             hold on
%         end
%     end
% end
