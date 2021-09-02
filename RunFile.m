%% MATLAB code for hybrid model of pannus formation in a PIP joint

% Code written by FR Macfarlane (frm3@st-andrews.ac.uk)

% NOTE: User will be asked if they want to reset the initial condition
% % Details can be found in Section 5. 
% % If this is the first time running type 'Y'

%% Section 1: Defaults (optional)
clear                                               % Clear all variables
close all                                           % Close all open figures
clc                                                 % Clear command line

%% Section 2: Set up spatial domain
Cwidth=0.04;                                        % Initial width of cartilage (cm)
Bwidth=0.1;                                         % Initial height of bone in domain (cm)
Cellint=0.004;                                      % Initial x position boundary of cells (x>0, x<xmin+Cellint & x>xmax-Cellint, x<xmax)
WSpace=Cellint*2;                                   % Whitespace at sides of joint (cm)
JointSpace=0.01;                                    % Initial gap between the bones (cm)
domxmax=0.3;                                        % Maximum x position of whole domain (cm)
domxmin=0;                                          % Minimum x position of whole domain (cm)
domymax=JointSpace+(2*Cwidth)+(2*Bwidth);           % Maximum y position of whole domain (cm)
domymin=0;                                          % Minimum y position of whole domain (cm)
Dom=[domxmin domxmax domymin domymax];              % Whole domain vector (for use in functions)
dec=3;                                              % Decimal point of space step (important for link between off lattice and on lattice)
dx=1*10^-(dec);                                     % Space step in x direction
dy=1*10^-(dec);                                     % Space step in y direction
x=Dom(1):dx:Dom(2);                                 % Define x coordinates of grid
y=Dom(3):dy:Dom(4);                                 % Define y coordinates of grid
x1=round(x,dec);                                    % MATLAB formatting - need to make sure x values are divisble by dx
y1=round(y,dec);                                    % MATLAB formatting - need to make sure y values are divisble by dy
Nx=length(x);                                       % Calculate the number of grid positions in x direction
Ny=length(y);                                       % Calculate the number of grid positions in y direction

dtcells=1*10^(-2);                                  % Time-step of cellular mechanisms (days)
dtpde=1*10^(-4);                                    % Time-step of chemical mechanisms (days)
T=20/dtpde;                                         % Set final time of simulations

PlotTime=1/dtpde;                                   % Plot the results every n time-steps
Data_No_Time=1;                                     % Save the numbers of cells and global concentrations/densities every n time-steps
Data_Pos_Time=5/dtpde;                              % Save the spatial bositions of cells and concentrations/densities every n time-steps                   

%% Section 3: Initial condition of cartilage, bone and MMPs

Cart=zeros(Nx,Ny);                                  % Define matrix to store cartilage density at each position
% Set vector for the boundaries of two areas of cartilage initially:
% [xmin xmax ymin1 ymax1 ymin2 ymax2]
CDom=[WSpace Dom(2)-WSpace Dom(3) Cwidth+Bwidth Dom(4)-Cwidth-Bwidth Dom(4)];
CDom=round(CDom,dec);                               % Round to ensure values will match x coordinates
CI(1)=find(x1==CDom(1));                            % Define CI(1) to be the index of Cart that matches the domain min x coordinates
CI(2)=find(x1==CDom(2));                            % Define CI(1) to be the index of Cart that matches the domain max x coordinates
CI(3)=find(y1==CDom(3));                            % Define CI(3) to be the index of Cart that matches the domain min area 1 y coordinates
CI(4)=find(y1==CDom(4));                            % Define CI(4) to be the index of Cart that matches the domain max area 1 y coordinates
CI(5)=find(y1==CDom(5));                            % Define CI(5) to be the index of Cart that matches the domain min area 2 y coordinates
CI(6)=find(y1==CDom(6));                            % Define CI(6) to be the index of Cart that matches the domain max area 2 y coordinates
Cart(CI(1):CI(2),CI(3):CI(4))=1;                    % Set cartilage to be one inside defined area 1
Cart(CI(1):CI(2),CI(5):CI(6))=1;                    % Set cartilage to be one inside defined area 2

Bone=zeros(Nx,Ny);                                  % Define matrix to store bone density at each position
% Set vector for the boundaries of two areas of bone initially:
% [xmin xmax ymin1 ymax1 ymin2 ymax2]
BDom=[WSpace+Cwidth Dom(2)-WSpace-Cwidth Dom(3) Bwidth Dom(4)-Bwidth Dom(4)];
BDom=round(BDom,dec);                               % Round to ensure values will match x coordinates
BI(1)=find(x1==BDom(1));                            % Define BI(1) to be the index of Bone that matches the domain min x coordinates
BI(2)=find(x1==BDom(2));                            % Define BI(1) to be the index of Bone that matches the domain max x coordinates
BI(3)=find(y1==BDom(3));                            % Define BI(3) to be the index of Bone that matches the domain min area 1 y coordinates
BI(4)=find(y1==BDom(4));                            % Define BI(4) to be the index of Bone that matches the domain max area 1 y coordinates
BI(5)=find(y1==BDom(5));                            % Define BI(5) to be the index of Bone that matches the domain min area 2 y coordinates
BI(6)=find(y1==BDom(6));                            % Define BI(6) to be the index of Bone that matches the domain max area 2 y coordinates
Bone(BI(1):BI(2),BI(3):BI(4))=1;                    % Set Bone to be one inside defined area 1
Bone(BI(1):BI(2),BI(5):BI(6))=1;                    % Set Bone to be one inside defined area 2

Cart=Cart-Bone;                                     % Set cartilage to be cart-bone (we don't want bone and cartilage on same positions)

InCart=sum(sum(Cart));                              % Find total initial density of cartilage
InBone=sum(sum(Bone));                              % Find total initial density of bone

DegChem=zeros(Nx,Ny);                               % Set initial degradation MMP chemical concentration to be zero on all lattice positions

%% Section 4: Set simulation parameters

RF=0.00065;                                         % Radius of fibroblasts (cm)
RM=0.00105;                                         % Radius of macrophages (cm)

ProbMovF=(dtcells/(dx*dx))*((8.64*10^(-7)));        % Probability of fibroblast moving at each time-step
ProbMovM=(dtcells/(dx*dx))*((8.64*10^(-7)));        % Probability of macrophage moving at each time-step

ProlifF=dtcells*((0.33));                           % Probability of fibroblast dividing at each time-step
ProlifM=dtcells*((0.33));                           % Probability of macrophage dividing at each time-step

DeathF=dtcells*((0.03));                            % Probability of fibroblast dying at each time-step 
DeathM=dtcells*((0.033));                           % Probability of macrophage dying at each time-step

DegFSecProdRate=dtcells.*(13.912*10^(-5));          % Concentration of degradation chemokine produced by fibroblasts at each time-step
DegMSecProdRate=dtcells.*(5.2717*10^(-5));          % Concentration of degradation chemokine produced by macrophages at each time-step

DiffDegChem=dtpde.*(6.59*10^(-2));                  % Diffusion rate of degradation chemokine
DecayDegChem=dtpde.*0.138;                          % Decay rate of degradation chemokine

DegCartRate=dtpde.*(4.44*10^3)/(1.15);              % Degradation rate of cartilage by degradation chemokine
DegBoneRate=0.1*dtpde.*(4.44*10^3)/(1.15);          % Degradation rate of bone by degradation chemokine

%% Section 5: Set initial condition of fibroblast and macrophages

NF=200;                                             % Initial number of fibroblasts
NM=200;                                             % Initial number of macrophages

centresF=zeros(NF,2);                               % Set up matrix to store x and y coordinates of fibroblasts centres
centresM=zeros(NM,2);                               % Set up matrix to store x and y coordinates of macrophage centres

IXV=[Dom(1) Dom(1)+Cellint Dom(2)-Cellint Dom(2)];  % Boundaries for initial x positions of cells
IYV=[Dom(3) Dom(4)];                                % Boundaries for initial y positions of cells

% Check if reseting the initial condition (if you wish to use the same
% initial condition as previous cases set as No (N))
prompt=input('Do you want to reset the initial condition? Type Y for yes, N for no: ','s');

if prompt =='Y'
    'Initial condition has been reset'
    
    % Run the initial condition function for larger cell population (macrophages) first (see function at end of file)
    [centresM]=initialcondition(centresM,centresF,RM,RF,IXV,IYV);
    
    % Run the initial condition function for other cell population (fibroblasts) (see function at end of file)
    [centresF]=initialcondition(centresF,centresM,RF,RM,IXV,IYV);
    
    % Save the initial conditino for future runs
    save('Initial_Condition.mat','centresF','centresM');
else 
    'Initial condition has been loaded from saved'
    
    % Load saved initial condition
    load('Initial_Condition.mat');
end

%% Section 6: Set up of storage vectors for saving data

NFT=[];                                             % Set up vector to store number of fibroblasts at each time-point
NMT=[];                                             % Set up vector to store number of macrophages at each time-point
TSumDegChem=[];                                     % Set up vector to store global MMP concentration at each time-point
TSumCart=[];                                        % Set up vector to store global cartilage density at each time-point
TSumBone=[];                                        % Set up vector to store global bone density at each time-point
NFT(1)=NF;                                          % Set initial number of fibroblasts to the first stored position
NMT(1)=NM;                                          % Set initial number of macrophages to the first stored position
TSumDegChem(1)=sum(sum(DegChem));                   % Set global MMP concentration to the first stored position
TSumCart(1)=sum(sum(Cart));                         % Set global cartilage density to the first stored position
TSumBone(1)=sum(sum(Bone));                         % Set global bone density to the first stored position

StorecentresF={};                                   % Set up vector to store spatial positions of fibroblasts at each time-point
StorecentresM={};                                   % Set up vector to store spatial positions of fibroblasts at each time-point
StoreDegChem={};                                    % Set up vector to store spatial positions of MMPs at each time-point
StoreCart={};                                       % Set up vector to store spatial positions of cartilage at each time-point
StoreBone={};                                       % Set up vector to store spatial positions of bone at each time-point
StorecentresF{1}=centresF;                          % Set initial spatial positions of fibroblasts to be the first stored values
StorecentresM{1}=centresM;                          % Set initial spatial positions of macrophages to be the first stored values
StoreDegChem{1}=DegChem;                            % Set initial spatial positions of MMPs to be the first stored values
StoreCart{1}=Cart;                                  % Set initial spatial positions of cartilage to be the first stored values
StoreBone{1}=Bone;                                  % Set initial spatial positions of bone to be the first stored values

%% Section 7: Run simulation time-loop

CTIME=dtcells/dtpde;                                % Set the time to run cell mechanisms (as the time-step is different than the chemical mechanisms)

% Set figure defaults
figure('units','normalized','outerposition',[0 0.3 0.5 0.5])  
set(0,'defaultAxesFontSize',20)
set(0, 'defaultAxesTickLabelInterpreter','latex')

for t=1:T                                           % For all time-steps until time t = T
    
    %% Section 7a: Cell movement (see function at end of file)
    if mod(t,CTIME)==0 && NF>0                      % If the time-step is one where cell mechanisms occur and fibroblasts are in the system
        
        % Run movement function for fibroblasts
        [centresF]=movement(centresF,centresM,NF,RF,RM,Dom,ProbMovF,Cart,Bone,x1,y1,dec,CDom,InCart,InBone,dx);
    end
    
    if mod(t,CTIME)==0 && NM>0                      % If the time-step is one where cell mechanisms occur and macrophages are in the system
        
        % Run movement function for macrophages
        [centresM]=movement(centresM,centresF,NM,RM,RF,Dom,ProbMovM,Cart,Bone,x1,y1,dec,CDom,InCart,InBone,dx);
    end
    
    %% Section 7b: Cell proliferation (see function at end of file)
    if mod(t,CTIME)==0 && NF>0                      % If the time-step is one where cell mechanisms occur and fibroblasts are in the system
        
        % Run proliferation function for fibroblasts
        [centresF,NF]=proliferation(centresF,centresM,NF,RF,RM,Dom,ProlifF,DeathF,Cart,Bone,x1,y1,dec,CDom,InCart,InBone);
    end
    
    if mod(t,CTIME)==0 && NM>0                      % If the time-step is one where cell mechanisms occur and macrophages are in the system
        
        % Run proliferation function for macrophages
        [centresM,NM]=proliferation(centresM,centresF,NM,RM,RF,Dom,ProlifM,DeathM,Cart,Bone,x1,y1,dec,CDom,InCart,InBone);
    end
    
    %% Section 7c: Chemical evolution and degradation of cartilage and bone (see function at end of file)
    
    % Run chemical evolution function for degradation chemokine
    [DegChem]=chemokineequation(DegChem,DiffDegChem,DecayDegChem,DegFSecProdRate,DegMSecProdRate,dtpde,dx,dy,Nx,Ny,centresF,centresM,dec,dtcells,t);
    
    % Impose cartilage degradation dependent on chemokine and degradation rate
    Cart=max(0,Cart-DegCartRate.*DegChem);
    
    % Impose bone degradation dependent on chemokine and degradation rate
    Bone=max(0,Bone-DegBoneRate.*DegChem);
    
    %% Section 7d: Plot results (optional, note plotting frequently does slow down simulation time)
    
    % Plot only when t is a multiple of PlotTime and there is cartilage
    % present
    if mod(t,PlotTime)==0 && sum(sum(Cart>0))>0
        % Plot cell and cartilage/bone densities in figure 1
        g=figure(1);
        clf
        axis(Dom)
        box on
        title(['Number of macrophages (red=' num2str(NM),') and fibroblasts (purple=' num2str(NF),') at time = ' num2str(t*dtpde)],'FontSize',24)
        [X,Y] = meshgrid(x,y);
        plCart=Cart;
        plCart(plCart<=0) = NaN;
        plBone=Bone;
        plBone(plBone<=0) = NaN;
        ax1 = axes;
        surf(ax1,X,Y,plCart','EdgeColor','none');
        caxis([0 max(max(Cart))])
        view(2)
        axis(Dom)
        ax2 = axes;
        hold on
        surf(ax2,X,Y,plBone','EdgeColor','none');
        view(2)
        caxis([0 max(max(Bone))])
        viscircles(centresF, RF*ones(NF,1),'Color',[0.5 0 0.5],'LineWidth',1);  % Plot fibroblasts centresF and radius RF (times number of cells)
        viscircles(centresM, RM*ones(NM,1),'Color','r','LineWidth',1);
        hold off
        axis(Dom)
        ax2.Visible = 'off';
        ax2.XTick = [];
        ax2.YTick = [];
        ax1.Visible = 'off';
        ax1.XTick = [];
        ax1.YTick = [];
        colormap(ax1,flipud(winter))
        colormap(ax2,flipud(gray))
        set([ax1,ax2],'Position',[0.1300 0.1100 0.7750 0.8150]);
        cb1 = colorbar(ax1,'Position',[0.92 .11 .01 .815]);
        cb2 = colorbar(ax2,'Position',[0.05 .11 .01 .815]);
        cb1.Label.String = 'Cartilage density (>0)';
        cb1.Label.FontSize=24;
        cb2.Label.String = 'Bone density (>0)';
        cb2.Label.FontSize=24;
        drawnow
        
        % Plot chemokine concentration in figure 2
        h=figure(2);
        pDegChem=DegChem;
        surf(X,Y,pDegChem','EdgeColor','none');
        title(['Degradation chemokine at time = ' num2str(t*dtpde)],'FontSize',24)
        axis(Dom)
        caxis([0 max(max(DegChem))])
        view(2)
        cb3=colorbar;
        cb3.Label.FontSize=24;
        cb3.Label.String = 'Concentration';
        box on
        drawnow
        
    end
    
    % Save the total cell numbers, global concentrations, global densities
    % every Data_No_Time time-steps
    if mod(t,Data_No_Time)==0
        NFT(end+1)=NF;
        NMT(end+1)=NM;
        TSumDegChem(end+1)=sum(sum(DegChem));
        TSumCart(end+1)=sum(sum(Cart));
        TSumBone(end+1)=sum(sum(Bone));
    end
    
    % Save the spatial positions of cell numbers, concentrations, densities
    % every Data_Pos_Time time-steps    
    if mod(t,Data_Pos_Time)==0
        StorecentresF{end+1}=centresF;
        StorecentresM{end+1}=centresM;
        StoreCart{end+1}=Cart;
        StoreBone{end+1}=Bone;
        StoreDegChem{end+1}=DegChem;
    end
    
end

% Save specific data to the file Data.mat
save('Data.mat','-v7.3','StorecentresF','StorecentresM','StoreCart','StoreBone','StoreDegChem','NFT','NMT','TSumDegChem','TSumCart','TSumBone');

%% Optional: Plot total numbers over time for the run (comment out if note required)

figure('units','normalized','outerposition',[0 0 1 1])

t=linspace(1,T*dtpde,T+1);

subplot(1,5,1)
plot(t,NFT,'Color',[0.5 0 0.5],'LineWidth',2)
box on
xlim([1 T*dtpde])
ylim([0 max(NFT)*1.1])
axis square
title('Fibroblasts','FontSize',24,'Interpreter','latex')
xlabel('Time (days)','FontSize',24,'Interpreter','latex')
ylabel('Number of Cells','FontSize',24,'Interpreter','latex')

subplot(1,5,2)
plot(t,NMT,'Color',[1 0 0],'LineWidth',2)
box on
xlim([1 T*dtpde])
ylim([0 max(NMT)*1.1])
axis square
title('Macrophages','FontSize',24,'Interpreter','latex')
xlabel('Time (days)','FontSize',24,'Interpreter','latex')
ylabel('Number of Cells','FontSize',24,'Interpreter','latex')

subplot(1,5,3)
plot(t,TSumDegChem,'Color',[0 1 0],'LineWidth',2)
box on
xlim([1 T*dtpde])
ylim([0 max(TSumDegChem)*1.1])
axis square
title('MMPs','FontSize',24,'Interpreter','latex')
xlabel('Time (days)','FontSize',24,'Interpreter','latex')
ylabel('Global Concentration','FontSize',24,'Interpreter','latex')

subplot(1,5,4)
plot(t,TSumCart,'Color',[0 0 1],'LineWidth',2)
box on
xlim([1 T*dtpde])
ylim([0 max(TSumCart)*1.1])
axis square
title('Cartilage','FontSize',24,'Interpreter','latex')
xlabel('Time (days)','FontSize',24,'Interpreter','latex')
ylabel('Global Density','FontSize',24,'Interpreter','latex')

subplot(1,5,5)
plot(t,TSumBone,'Color',[0 0 0],'LineWidth',2)
box on
xlim([1 T*dtpde])
ylim([0 max(TSumBone)*1.1])
axis square
title('Bone','FontSize',24,'Interpreter','latex')
xlabel('Time (days)','FontSize',24,'Interpreter','latex')
ylabel('Global Density','FontSize',24,'Interpreter','latex')


%% End of Run File %


%% FUNCTIONS BELOW
%% Function 1: Initial condition of cells

function [centresA]=initialcondition(centresA,centresB,RA,RB,IXV,IYV)
% Function to determine the initial condition of cells of type A, ensures
% no overlap with other cell types (e.g other cells in population A or
% population B)
%
% INPUT:
% % centresA = centre coordinates of population A
% % centresB = centre coordinates of population B
% % RA = radius of cells of population A
% % RB = radius of cells of population B
% % IXV = vector that stores x limits of the boundary for initial condition
% % IXY = vector that stores y limits of the boundary for initial condition
%
% OUTPUT:
% % centresA = updated centre coordinates of population A
%% Function code time approx 0.008 seconds (efficient enough as not run each time-step)
%%

% Divide the initial population into 2 (so that half the cells begin on
% LHS of domain and half on the RHS)
Z1=round(length(centresA)/2);


% For first half of cells
% (NOTE: Within for loops is faster method than full vectors, as if full vectors all cells reset every time one is on same position)
for a=1:Z1
    
    % Set an index k
    k=1;
    
    % Set while loop
    while k==1
        
        % Set x coordinate to be a random value within the desired boundary
        centresA(a,1)= (IXV(1)+RA) + ((IXV(2)-RA)-(IXV(1)+RA))*rand;
        
        % Set y coordinate to be a random value within the desired boundary
        centresA(a,2)=(IYV(1)+RA) + ((IYV(2)-RA)-(IYV(1)+RA))*rand;
        
        % Set the new coordinate vector to be A
        A=centresA(a,:);
        
        % Set B to be the coordinates of population A omitting the current cell
        % under consideration
        B=centresA;
        B(a,:)=[];
        
        % Set B to be the coordnates of the other population
        C=centresB;
        
        % Find the distance between current cells position and all other cells
        % in population A
        dist1=sqrt((A(1)-B(:,1)).^2+(A(2)-B(:,2)).^2);
        
        % Find the distance between current cells position and all cells
        % in population B
        dist2=sqrt((A(1)-C(:,1)).^2+(A(2)-C(:,2)).^2);
        
        % Find number of population A cells that are within the distance 2RA of cell
        F1=length(find(dist1<(RA+RA)));
        
        % Find number of population B cells that are within the distance RA+RB of cell
        F2=length(find(dist2<(RA+RB)));
        
        % If their is no overlap of cells then this initial condition is fine
        % for this cell and k updates to stop while loop
        
        if F1==0 && F2==0
            k=0;
        end
        
    end
end

% For second half of cells
for b=(Z1+1):length(centresA)
    % Set an index k
    k=1;
    
    % Set while loop
    while k==1
        
        % Set x coordinate to be a random value within the desired boundary
        centresA(b,1)= (IXV(3)+RA) + ((IXV(4)-RA)-(IXV(3)+RA))*rand;
        
        % Set y coordinate to be a random value within the desired boundary
        centresA(b,2)=(IYV(1)+RA) + ((IYV(2)-RA)-(IYV(1)+RA))*rand;
        
        % Set the new coordinate vector to be A
        A=centresA(b,:);
        
        % Set B to be the coordinates of population A omitting the current cell
        % under consideration
        B=centresA;
        B(b,:)=[];
        
        % Set B to be the coordnates of the other population
        C=centresB;
        
        % Find the distance between current cells position and all other cells
        % in population A
        dist1=sqrt((A(1)-B(:,1)).^2+(A(2)-B(:,2)).^2);
        
        % Find the distance between current cells position and all cells
        % in population B
        dist2=sqrt((A(1)-C(:,1)).^2+(A(2)-C(:,2)).^2);
        
        % Find number of population A cells that are within the distance 2RA of cell
        F1=length(find(dist1<(RA+RA)));
        
        % Find number of population B cells that are within the distance RA+RB of cell
        F2=length(find(dist2<(RA+RB)));
        
        % If there is no overlap of cells then this initial condition is fine
        % for this cell and k updates to stop while loop
        if F1==0 && F2==0
            k=0;
        end
        
        
    end
end
end
% End of Function %

%% Function 2: MOvement of cells

function [centresA]=movement(centresA,centresB,NA,RA,RB,Dom,ProbMovA,Cart,Bone,x,y,dec,CDom,InCart,InBone,dx)
% Function to allow the movement of cells in population A without overlap
% with other cells in population A, cells in population B or areas where cartilage/bone
% density is non-zero
%
% INPUT:
% % centresA = centre coordinates of population A
% % centresB = centre coordinates of population B
% % NA = number of cells in population A
% % RA = radius of cells of population A
% % RB = radius of cells of population B
% % Dom = vector that stores the boundaries of the domain
% % ProbA = probability of movement for a cell of type A
% % Cart = density of cartilage
% % Bone = density of bone
% % x = x coordinates of lattice (cells off lattice, cartilage/bone on lattice)
% % y = y coordinates of lattice (cells off lattice, cartilage/bone on lattice)
% % dec = how many decimal places we require to find approx position of cells on lattice
% % CDom = domain of cartilage initial condition
% % InCart = total sum of initial cartilage denisty in full domain
% % InBone = total sum of initial bone denisty in full domain
% % dx = space-step on grid (for movement length)
% OUTPUT:
% % centresA = updated centre coordinates of population A


% Calculate a random value K for each cell in the population
K=rand(NA,1);

% Find the index values of the cells for which K is less than or equal to
% the probability of movement - These are the cells that will move in this
% time-step
Index=find(K<=ProbMovA);

% For all cells that are permitted to move in this time-step (Note doing
% this method rather than manipulating the vectors outwith the for loop
% appears faster)
if ~isempty(Index)==1
    for a=1:length(Index)
        
        % Set matrix Orig to store the original x and y coordinates of the
        % cells that are permitted
        Orig=centresA(Index(a),:);
        
        % Pick a random angle as the direction of movement
        A1 = (2*pi).*rand;
        
        % The desired new cell centre position is set as V
        V(1)=Orig(1)+dx * cos(A1);
        V(2)=Orig(2)+dx * sin(A1);
        
        % Check desired new position is within boundary of full domain (±RA to
        % ensure whole diameter of cell is within domain)
        if V(1)>(Dom(1)+RA) && V(2)>(Dom(3)+RA) && V(1)<=(Dom(2)-RA) && V(2)<=(Dom(4)-RA)
            
            % Set W to be the coordinates of population A
            W=centresA;
            W(Index(a),:)=NaN;
            
            % Find distance between desired position and all other cells in population A
            dist=sqrt((V(1)-W(:,1)).^2+(V(2)-W(:,2)).^2);
            
            % Find if any cells are within the distance of two radii of cell type A of desired
            % position
            F1=length(find(dist<(RA+RA)));
            
            if length(centresB)>0
                % Set Z to be coordinates of cells in population B
                Z=centresB;
                
                % Find distance between desired position and all cells in
                % population B
                dist=sqrt((V(1)-Z(:,1)).^2+(V(2)-Z(:,2)).^2);
                
                % Find if any cells are within the distance the sum of the radii of both cells of desired
                % position
                F2=length(find(dist<(RA+RB)));
            else
                F2=0;
            end
            % Set CVal to check cartilage density within the neighbourhood
            % of the desired cell position
            CVal=0;
            
            % Set BVal to check bone density within the neighbourhood
            % of the desired cell position
            BVal=0;
            
            % Set a new variable as a method of checking local levels of
            % cartilage/bone
            Check=0;
            
            % If the desired new cell position (including circumference of cell) is within the original area of
            % cartilage set Check to be 1
            if V(1)>=(CDom(1)-RA) && V(1)<=(CDom(2)+RA) && V(2)>=(CDom(3)-RA) && V(2)<=(CDom(4)+RA)
                Check=1;
            elseif V(1)>=(CDom(1)-RA) && V(1)<=(CDom(2)+RA) && V(2)>=(CDom(5)-RA) && V(2)<=(CDom(6)+RA)
                Check=1;
            end
            
            % If the desired new cell position is within the original domain
            % of cartilage
            if Check==1
                % If there has been degradation of cartilage in the domain
                if sum(sum(Cart))<InCart
                    % Value of n set according to the spatial scale chosen,
                    % this is used to link off-lattice to on-lattice
                    if dec==3
                        n=3;
                    else
                        n=100;
                    end
                    % Pick a range of angles to find points on circumference of the cell (see Test_theta file)
                    th = 0:pi/n:2*pi;
                    % Find x, y values of points on the circumference of the cell
                    xpol = RA * cos(th) + V(1);
                    ypol = RA * sin(th) + V(2);
                    % Round the x, y values of points on the circumference of
                    % the cell to become the nearest grid positions
                    xpol=round(xpol,dec);
                    ypol=round(ypol,dec);
                    % For all x, y positions on the circumference on the cell
                    for b=1:length(xpol)
                        % Find the x,y indexes of the x, y coordinates
                        IX=find(x==xpol(b));
                        IY=find(y==ypol(b));
                        % Find the density of cartilage at each point
                        IN=Cart(IX,IY);
                        % If cartilage density is non-zero then set CVal=1
                        % (note, this is not accumulative so if any (or all)
                        % positions are non-zero CVal=1 (binary)
                        if IN>0
                            CVal=1;
                        end
                    end
                    % If there has been no degradation of cartilage then the density here will be non-zero, so set CVal=1
                else
                    CVal=1;
                end
                
                % If there has been degradation of bone in the domain
                if sum(sum(Bone))<InBone
                    % Value of n set according to the spatial scale chosen,
                    % this is used to link off-lattice to on-lattice
                    if dec==3
                        n=3;
                    else
                        n=100;
                    end
                    % Pick a range of angles to find points on circumference of the cell (see Test_theta file)
                    th = 0:pi/n:2*pi;
                    % Find x, y values of points on the circumference of the cell
                    xpol = RA * cos(th) + V(1);
                    ypol = RA * sin(th) + V(2);
                    % Round the x, y values of points on the circumference of
                    % the cell to become the nearest grid positions
                    xpol=round(xpol,dec);
                    ypol=round(ypol,dec);
                    % For all x, y positions on the circumference on the cell
                    for b=1:length(xpol)
                        % Find the x,y indexes of the x, y coordinates
                        IX=find(x==xpol(b));
                        IY=find(y==ypol(b));
                        % Find the density of bone at each point
                        IN=Bone(IX,IY);
                        % If bone density is non-zero then set CVal=1
                        % (note, this is not accumulative so if any (or all)
                        % positions are non-zero BVal=1 (binary)
                        if IN>0
                            BVal=1;
                        end
                    end
                    % If there has been no degradation of bone then the density here will be non-zero, so set BVal=1
                else
                    BVal=1;
                end
            end
            
            
            % If there is no overlap with cells in population A, population B
            % or cartilage/bone densities the cell centre position updates to
            % the new desired position V
            if F1==0 && F2==0 && CVal==0 && BVal==0 % Note F1==1 as the original cell position should be overlap
                
                centresA(Index(a),:)=V;
                
            end
        end
    end
end
end
% End of Function %

%% Function 3: Proliferation of cells

function [centresA,NA]=proliferation(centresA,centresB,NA,RA,RB,Dom,ProlifA,DeathA,Cart,Bone,x,y,dec,CDom,InCart,InBone)
% Function to allow the death and division of cells in population A without overlap
% with other cells in population A, cells in population B or areas where cartilage/bone
% density is non-zero
%
% INPUT:
% % centresA = centre coordinates of population A
% % centresB = centre coordinates of population B
% % NA = number of cells in population A
% % RA = radius of cells of population A
% % RB = radius of cells of population B
% % Dom = vector that stores the boundaries of the domain
% % ProlifA = probability of division for a cell of type A
% % DeathA = probability of death for a cell of type A
% % Cart = density of cartilage
% % Bone = density of bone
% % x = x coordinates of lattice (cells off lattice, cartilage/bone on lattice)
% % y = y coordinates of lattice (cells off lattice, cartilage/bone on lattice)
% % dec = how many decimal places we require to find approx position of cells on lattice
% % CDom = domain of cartilage initial condition
% % InCart = total sum of initial cartilage denisty in full domain
% % InBone = total sum of initial bone denisty in full domain
%
% OUTPUT:
% % centresA = updated centre coordinates of population A
% % NA = updated number of cells in position A

%% Updated method

% Calculate a random value K for each cell in the population
K=rand(NA,1);

% Find the index values of the cells for which K is less than or equal to
% the probability of death - These are the cells that will die in this
% time-step
Index1=find(K<=DeathA);

% Find the index values of the cells for which K is between the probability of death and probability of division plus death - These are the cells that will divide in this
% time-step
Index2=find(K>DeathA & K<=(ProlifA+DeathA));

% Set to NaN the centre positions of the cells that will be removed through
% cell death

centresA(Index1,:)=NaN;

% For all cells that are permitted to divide in this time-step (Note doing
% this method rather than manipulating the vectors outwith the for loop
% appears faster)
if ~isempty(Index2)==1
    for a=1:length(Index2)
        
        % Set matrix Orig to store the original x and y coordinates of the
        % cells that are permitted to divide
        Orig=centresA(Index2(a),:);
        
        % Pick a random angle as the direction of division
        A1 = (2*pi).*rand;
        
        % The desired new cell centre position is set as V
        V(1)=Orig(1)+2*RA * cos(A1);
        V(2)=Orig(2)+2*RA * sin(A1);
        
        
        
        % Check desired new position is within boundary of full domain (±RA to
        % ensure whole diameter of cell is within domain)
        if V(1)>(Dom(1)+RA) && V(2)>(Dom(3)+RA) && V(1)<=(Dom(2)-RA) && V(2)<=(Dom(4)-RA)
            
            % Set W to be the coordinates of population A
            W=centresA;
            
            % Find distance between desired position and all other cells in population A
            dist=sqrt((V(1)-W(:,1)).^2+(V(2)-W(:,2)).^2);
            
            % Find if any cells are within the distance of two radii of cell type A of desired
            % position
            F1=length(find(dist<(RA+RA)));
            
            if length(centresB)>0
                % Set Z to be coordinates of cells in population B
                Z=centresB;
                
                % Find distance between desired position and all cells in
                % population B
                dist=sqrt((V(1)-Z(:,1)).^2+(V(2)-Z(:,2)).^2);
                
                % Find if any cells are within the distance the sum of the radii of both cells of desired
                % position
                F2=length(find(dist<(RA+RB)));
            else
                F2=0;
            end
            % Set CVal to check cartilage density within the neighbourhood
            % of the desired cell position
            CVal=0;
            
            % Set BVal to check bone density within the neighbourhood
            % of the desired cell position
            BVal=0;
            
            % Set a new variable as a method of checking local levels of
            % cartilage/bone
            Check=0;
            
            % If the desired new cell position (including circumference of cell) is within the original area of
            % cartilage set Check to be 1
            if V(1)>=(CDom(1)-RA) && V(1)<=(CDom(2)+RA) && V(2)>=(CDom(3)-RA) && V(2)<=(CDom(4)+RA)
                Check=1;
            elseif V(1)>=(CDom(1)-RA) && V(1)<=(CDom(2)+RA) && V(2)>=(CDom(5)-RA) && V(2)<=(CDom(6)+RA)
                Check=1;
            end
            
            % If the desired new cell position is within the original domain
            % of cartilage
            if Check==1
                % If there has been degradation of cartilage in the domain
                if sum(sum(Cart))<InCart
                    % Value of n set according to the spatial scale chosen,
                    % this is used to link off-lattice to on-lattice
                    if dec==3
                        n=3;
                    else
                        n=100;
                    end
                    % Pick a range of angles to find points on circumference of the cell (see Test_theta file)
                    th = 0:pi/n:2*pi;
                    % Find x, y values of points on the circumference of the cell
                    xpol = RA * cos(th) + V(1);
                    ypol = RA * sin(th) + V(2);
                    % Round the x, y values of points on the circumference of
                    % the cell to become the nearest grid positions
                    xpol=round(xpol,dec);
                    ypol=round(ypol,dec);
                    % For all x, y positions on the circumference on the cell
                    for b=1:length(xpol)
                        % Find the x,y indexes of the x, y coordinates
                        IX=find(x==xpol(b));
                        IY=find(y==ypol(b));
                        % Find the density of cartilage at each point
                        IN=Cart(IX,IY);
                        % If cartilage density is non-zero then set CVal=1
                        % (note, this is not accumulative so if any (or all)
                        % positions are non-zero CVal=1 (binary)
                        if IN>0
                            CVal=1;
                        end
                    end
                    % If there has been no degradation of cartilage then the density here will be non-zero, so set CVal=1
                else
                    CVal=1;
                end
                
                % If there has been degradation of bone in the domain
                if sum(sum(Bone))<InBone
                    % Value of n set according to the spatial scale chosen,
                    % this is used to link off-lattice to on-lattice
                    if dec==3
                        n=3;
                    else
                        n=100;
                    end
                    % Pick a range of angles to find points on circumference of the cell (see Test_theta file)
                    th = 0:pi/n:2*pi;
                    % Find x, y values of points on the circumference of the cell
                    xpol = RA * cos(th) + V(1);
                    ypol = RA * sin(th) + V(2);
                    % Round the x, y values of points on the circumference of
                    % the cell to become the nearest grid positions
                    xpol=round(xpol,dec);
                    ypol=round(ypol,dec);
                    % For all x, y positions on the circumference on the cell
                    for b=1:length(xpol)
                        % Find the x,y indexes of the x, y coordinates
                        IX=find(x==xpol(b));
                        IY=find(y==ypol(b));
                        % Find the density of bone at each point
                        IN=Bone(IX,IY);
                        % If bone density is non-zero then set CVal=1
                        % (note, this is not accumulative so if any (or all)
                        % positions are non-zero BVal=1 (binary)
                        if IN>0
                            BVal=1;
                        end
                    end
                    % If there has been no degradation of bone then the density here will be non-zero, so set BVal=1
                else
                    BVal=1;
                end
            end
            
            
            % If there is no overlap with cells in population A, population B
            % or cartilage/bone densities the new cell centre position is added to
            % the new desired position V
            if F1==0 && F2==0 && CVal==0 && BVal==0
                
                centresA(end+1,:)=V; % Update centres to include new cell
                
            end
        end
    end
end

% Update vectors and remove dead cells
P1=centresA(:,1);
P2=centresA(:,2);
P1(isnan(P1))=[];
P2(isnan(P2))=[];
centresA=[P1 P2];

% Update total number of cells
NA=~isempty(centresA)*size(centresA,1);

end
% End of Function %

%% Function 4: MMPs evolution

function [NChem]=chemokineequation(NChem,DiffNChem,DecayNChem,NSecFProdRate,NSecMProdRate,dtpde,dx,dy,Nx,Ny,centersF,centersM,dec,dtcells,t)
% Function to allow the evolution of the specific chemokine
%
% INPUT:
% % NChem = chemokine concentration on each lattice position
% % DiffNChem = diffusion rate of chemokine
% % DecayNChem = decay rate of chemokine
% % NSecFProdRate = production rate of chemokine by fibroblasts
% % NSecMProdRate = production rate of chemokine by macrophages
% % dt = tims-step of simulations
% % dx = lattice space length in the x direction
% % dy = lattice space length in the y direction
% % Nx = number of grid positions in x direction
% % Ny = number of grid positions in y direction
% % centersF = cell centre positions of fibroblasts
% % centersM = cell centre positions of macrophages
% % dec = how many decimal places we require to find approx position of cells on lattice
% % dtcells = time-step of cell mechanisms
% % t= current time-step

% OUTPUT:
% % NChem = updated MMP chemokine concentration on each lattice position

%% Original Method

% Set alpha to be a matrix to calculate production rates
Alpha=zeros(Nx,Ny);


if mod(t,dtcells/dtpde)==0
    
    if length(centersF)>0
        val2=round(centersF,dec);
        val2x=round((val2(:,1)/dx)+1);
        val2y=round((val2(:,2)/dy)+1);
        for a=1:length(val2x)
            Alpha(val2x(a),val2y(a))=NSecFProdRate/(dx^2);
        end
    end
    
    
    if length(centersM)>0
        % Find closest lattice positions to cell centres (note this infers cells
        % produce chemokine at cell centre and not on circumference of cells)
        val1=round(centersM,dec);
        
        
        % Find indexes of positions
        val1x=round((val1(:,1)/dx)+1);
        val1y=round((val1(:,2)/dy)+1);
        
        
        % For each index set Alpha to be the production rate times density of cells on that position (note 1/dx^2)
        % (Note, these are not probabilities but deterministic)
        for a=1:length(val1x)
            Alpha(val1x(a),val1y(a))=NSecMProdRate/(dx^2);
        end
        
    end
end
% Store old concentration at each position
uold=NChem;

% Define local interaction terms, production and decay
f=Alpha(2:Nx-1,2:Ny-1).*(1-uold(2:Nx-1,2:Ny-1))-DecayNChem.*uold(2:Nx-1,2:Ny-1);

% Define diffusion terms in x direction
Lapu=DiffNChem.*(dtpde/(dx*dx)).*(uold(1:Nx-2,2:Ny-1)+uold(3:Nx,2:Ny-1)-2*uold(2:Nx-1,2:Ny-1));

% Define diffusion terms in y direction
Lapyu=DiffNChem.*(dtpde/(dy*dy)).*(uold(2:Nx-1,1:Ny-2)+uold(2:Nx-1,3:Ny)-2*uold(2:Nx-1,2:Ny-1));

% Update vectors
NChem(2:Nx-1,2:Ny-1)=uold(2:Nx-1,2:Ny-1)+Lapu+Lapyu+(dtpde.*f);

% Impose no-fluxboundary conditions
NChem(1,2:Ny-1)=NChem(2,2:Ny-1);
NChem(Nx,2:Ny-1)=NChem(Nx-1,2:Ny-1);
NChem(2:Nx-1,1)=NChem(2:Nx-1,2);
NChem(2:Nx-1,Ny)=NChem(2:Nx-1,Ny-1);
NChem(1,1)=NChem(2,2);
NChem(Nx,1)=NChem(Nx-1,2);
NChem(1,Ny)=NChem(2,Ny-1);
NChem(Nx,Ny)=NChem(Nx-1,Ny-1);

end
% end of function %

%% END OF FILE %%
