clear all,close all
rng('shuffle');

%% Particles for ISW
H=20; % total depth
L=500; %domain length
% x0=350; % starting x position

%% Particle parameters %%
numsteps=200; % number of steps integrated over between outs
numouts=40; %default numouts = 1000, saved outputs
dt=0.1; % timestep
ts=zeros(numouts,1); % array of sample times

g = 9.81; %gravity constant
settling= -0.001; % settling velocity

% load in the ISW info
load isw_infonou_hires.mat
sz=size(uisw);Nz=sz(1);Nx=sz(2);
zmin=min(zisw(:));
if exist('cisw') == 0
    cisw = c;
end

%% Inline Functions: flow from the ISW data
% one trick I use is to zero out the velocity if we are below the bottom
u=@(xq,zq,t,x,z,uvals) (zq>zmin).*interp2(x,z,uvals,xq,zq,'cubic',0)+0;  
w=@(xq,zq,t,x,z,wvals) (zq>zmin).*interp2(x,z,wvals,xq,zq,'cubic',0)+0;
% x,z are the sample coordinates
% uvals/wvals are the function value at the sample points
% xq,zq are the coordinates to be interpolated at

%% initial conditions
initial_condition = 1;
if initial_condition == 1
    % initial conditions for a square
    x=[];z=[];
    for zstart=-0.5:-2:-17.5
        xg1=linspace(350,360,20);
        zg1=linspace(zstart-1,zstart,20);
        [dumxx dumzz]=meshgrid(xg1,zg1);
        x=[x;dumxx(:)];
        z=[z;dumzz(:)];
        numens=length(x);
    end
elseif initial_condition == 2
    % initial conditions for two squares horizontally separated
    x=[];z=[];
    deltax = 5;
    for zstart=[-0.5,-4.5,-8.5,-12.5, -16.5]
        xg1=[linspace(350,360,20),linspace(360+deltax,370+deltax,20)];
        zg1=linspace(zstart-1,zstart,20);
        [dumxx dumzz]=meshgrid(xg1,zg1);
        x=[x;dumxx(:)];
        z=[z;dumzz(:)];
        numens=length(x);
    end
end
    dz0=abs(zg1(2)-zg1(1)); % differential change in z
    dx0=abs(xg1(2)-xg1(1)); % differential change in x

%% red noise parameters
memt=1;
rho=exp(-dt/memt); rhoc=sqrt(1-rho*rho);
mysig = 0.01;
prevrand1=zeros(numens,1);
prevrand2=prevrand1;

xs=zeros(numens,numouts+1);
zs=zeros(numens,numouts+1);
t=0;
eps=0.5*min([dz0, dx0]); % sets a tolerance distance for particle interactions
mydistmat=zeros(numens); % initialize adjacency matrix
mydistmatcum=mydistmat; % initialize adjacency matrix with edge weightings
zero_diag = ones(numens) - eye(numens); % used to ensure A_ii = 0 for all i
%mydistmat_changes=zeros(numens,numens,numouts);

%% The double loop over outputs and timesteps
for ii=1:numouts
    % storage
    xs(:,ii)=x;
    zs(:,ii)=z;
    ts(ii)=t;  
    % inner loop
    for jj=1:numsteps
        % noise
        nowrand1=randn(numens,1)*mysig;
        pertnow1=(nowrand1*rhoc+rho*prevrand1); prevrand1=pertnow1;
        nowrand2=randn(numens,1)*mysig;
        pertnow2=(nowrand2*rhoc+rho*prevrand2); prevrand2=pertnow2;
        
        %extracting local data for interpolation for speed up
        ileft=max(floor(Nx*min(x(:))/L)-10,1); % x index of leftmost particle - 10
        iright=min(floor(Nx*max(x(:))/L)+10,Nx); % x index of rightmost particle + 10
        jbot=min(floor(-Nz*min(z(:))/H)+10,Nz); % z index of bottommost particle + 10
        jtop=max(floor(-Nz*max(z(:))/H)-10,1); % z index of topmost particle - 10
        x_interp=xisw(jtop:jbot,ileft:iright); % sample x coordinates
        z_interp=zisw(jtop:jbot,ileft:iright); % sample z coordinates
        u_interp=uisw(jtop:jbot,ileft:iright); % sample u velocity values
        w_interp=wisw(jtop:jbot,ileft:iright); % sample w velocity values
        rhx=u(x,z,t,x_interp,z_interp,u_interp)+(z>zmin).*pertnow1-cisw;
            % interpolated u velocity + velocity perturbation - constant
        x=mod(x+dt*rhx,L); % updated x position
        rhz=w(x,z,t,x_interp,z_interp,w_interp)+settling.*(z>zmin)+(z>zmin).*pertnow2;
            % interpolated w velocity + settling/gravity + velocity perturbation
        z=z+dt*rhz; % updated z position
        z=z.*(z<0)-z.*(z>=0); % ensure z is negative
        t=t+dt; % updated time value
    end
    % forming the adjacency matrix
    distv=pdist([x,z]); % find interparticle distance between each pair
    distvind=(distv<eps); % note when particles interact/come within tolerance distance of eachother
    mydistmat_now=squareform(distvind); % turn into matrix where 1 indicates particle interaction
    mydistmat=(mydistmat+mydistmat_now.*(mydistmat==0)).*zero_diag;
        % matrix stores 1 if particles interacted at any time, 0 otherwise
    mydistmatcum=(mydistmatcum+mydistmat_now).*zero_diag;
        % matrix stores how many times particles interacted
    mydegnow(ii,:)=sum(mydistmat,2);
        % matrix storing the degree of each particle at time ii
    mydegnowinst(ii,:)=sum(mydistmat_now,2);
        % matrix storing how many interactions occurred at time ii (the instantaneous degree)
    %mydistmat_changes(:,:,ii)=mydistmat_now;
end
% One last storage
xs(:,numouts+1)=x;
zs(:,numouts+1)=z;
ts(numouts+1)=t;


%% FIGURES


%% MAX DEGREE TIME SERIES
[mydegnowmax mydegnowmaxi]=max(mydegnow,[],2);
    % mydegnowmax is an array of the maximum degree at each timestep
    % mydegnowmaxi is an array of the indices corresponding to mydegnowmax
[mydeginstmax mydeginstmaxi]=max(mydegnowinst,[],2); % same for the instantaneous degree
figure(2)
clf
betterplots
plot(ts(2:end),mydegnowmax,'kp-',ts(2:end),mydeginstmax,'bo-')
ylabel('max degree')
xlabel('time')
legend('cumulative','instantaneous')
grid on
savefig('outputs\maxdegree.fig')


%% WAVE INDUCED CURRENTS AND PARTICLE TRAJECTORIES
figure(1)
clf
betterplots
colormap darkjet
pcolor(xisw,zisw,uisw),shading interp
%pcolor(xisw,zisw,wisw),shading interp
hold on
if initial_condition == 1
    plot(xs(:,1:5:end),zs(:,1:5:end),'k.')
elseif initial_condition == 2
    plot(xs(:,1:13:end),zs(:,1:13:end),'k.')
end
xlabel('x (m)')
ylabel('z (m)')
title('Horizontal Wave-Induced Velocities')
%title('Vertical Wave-Induced Velocities')
colorbar
%plot(xs(3601:7200,1:5:41),zs(3601:7200,1:5:41),'w.')
savefig('outputs\particle_transport.fig')
%print('wave_induced_velocities','-dtiff','-r650')


%% PARTICLE DEGREE REPRESENTED BY MARKERSIZE
figure(3)
clf
betterplots
colormap darkjet
pcolor(xisw,zisw,uisw),shading interp
ii=floor(numouts/2);
hold on
for jj=1:size(xs,1)
    plot(xs(jj,ii+1),zs(jj,ii+1),'o','Markersize',2*mydegnow(ii,jj)+1,'Color',[1 1 1]*mydegnow(ii,jj)/max(10,max(mydegnow(:))))
    % at timestep ii, plot particles during transport, but the size of the
    % circle is proportional to the degree of that particle
end
savefig('outputs\degree_during_transport.fig')


%% INSTANTANEOUS PARTICLE DEGREE REPRESENTED BY MARKERSIZE
figure(4)
clf
betterplots
colormap darkjet
pcolor(xisw,zisw,uisw),shading interp
ii=floor(numouts/2);
hold on
for jj=1:3600
    plot(xs(jj,ii+1),zs(jj,ii+1),'o','Markersize',2*mydegnowinst(ii,jj)+1,'Color',[1 1 1]*mydegnowinst(ii,jj)/max(10,max(mydegnowinst(:))))
    % at timestep ii, plot particles during transport, but the size of the
    % circle is proportional to the instantaneous degree of that particle
end
savefig('outputs\instdegree_during_transport.fig')


%% POSITIONS OF THREE LARGEST CONNECTED CLUSTERS
G=graph(mydistmat);
[bin binsize]=conncomp(G); % finds all connected graphs in G, stores degree of all graphs in binsize
idx = binsize(bin) == max(binsize); % finds all nodes of the connected graph of maximum degree
sg = subgraph(G, idx); % creates sugraph with max degree

both_init_final = 1;
figure(12)
clf
betterplots
colormap darkjet
pcolor(xisw,zisw,uisw),shading interp
hold on
[aa ai]=maxk(binsize,3);
    % aa stores 3 largest degree values of connected graphs
    % ai stores corresponding indices
if both_init_final == 1
    x0=xs(:,1);
    z0=zs(:,1);
    myind=find(bin==ai(1));
    plot(x(myind),z(myind),'k.','MarkerSize',12) % connected graph of largest degree
    plot(x0(myind),z0(myind),'ko','MarkerSize',12)
    myind=find(bin==ai(2));
    plot(x(myind),z(myind),'w.','MarkerSize',12) % second largest
    plot(x0(myind),z0(myind),'wo','MarkerSize',12)
    myind=find(bin==ai(3));
    plot(x(myind),z(myind),'m.','MarkerSize',12) % third largest
    plot(x0(myind),z0(myind),'mo','MarkerSize',12)
else
    myind=find(bin==ai(1));
    plot(x(myind),z(myind),'k.','MarkerSize',12) % connected graph of largest degree
    myind=find(bin==ai(2));
    plot(x(myind),z(myind),'w.','MarkerSize',12) % second largest
    myind=find(bin==ai(3));
    plot(x(myind),z(myind),'m.','MarkerSize',12) % third largest 
end
title('Initial & Final Positions of Largest Connected Clusters')
xlabel('x (m)')
ylabel('z (m)')
drawnow
colorbar
savefig('outputs\three_largest_clusters.fig')
%print('three_largest_clusters','-dtiff','-r650')


%% PLOTTING CALCULATIONS
% % find three largest clusters
% G=graph(mydistmat);
% [bin binsize]=conncomp(G);
% [aa ai]=maxk(binsize,3);
% array to store cluster sizes
cluster_size = zeros(3,numouts+1);
% cell array to store average positions of clusters
meanx = {};
meanz = {};
% cell array to store particle degree of three largest clusters at final time
cluster_degree = {};
% cell array to store adjacency matrices of three largest clusters
cluster_adjmat = {};
% cell array to store particle local clustering coefficient of three largest clusters at final time
cluster_coeff = {};
% cell array to store indices of three largest clusters
cluster_ind = {};
for jj_cluster = 1:3
    % isolate positions of particles in cluster
    myind=find(bin==ai(jj_cluster));
    cluster_ind{jj_cluster} = myind;
    xs_cluster = xs(myind,:);
    zs_cluster = zs(myind,:);
    % initialize adjacency matrix
    mydistmat_cluster=zeros(size(myind));
    for ii = 1:numouts+1
        % form adjacency matrix of first largest cluster
        distv=pdist([xs_cluster(:,ii),zs_cluster(:,ii)]);
        distvind=(distv<eps);
        mydistmat_now=squareform(distvind);
        mydistmat_cluster=mydistmat_cluster+mydistmat_now.*(mydistmat_cluster==0);
        % calculate cluster size
        G_cluster=graph(mydistmat_cluster);
        [bin_cluster, binsize_cluster]=conncomp(G_cluster);
        cluster_size(jj_cluster,ii)=max(binsize_cluster);
    end
    % calculate average positions of cluster
    meanx{jj_cluster} = mean(xs_cluster,1);
    meanz{jj_cluster} = mean(zs_cluster,1);

    % store vector of particle degree
    cluster_degree{jj_cluster} = mydegnow(numouts,myind);

    % store adjacency matrix
    cluster_adjmat{jj_cluster} = mydistmat_cluster;

    num_cluster = size(xs_cluster,1);
    local_coeff = zeros(num_cluster,1);
    for ii = 1:num_cluster
        % calculate number of neighbours by summing row of adjacency matrix
        num_neighbours = sum(mydistmat_cluster(ii,:));

        % the local clustering coefficient is zero when the number of neighbours is 0 or 1
        if num_neighbours == 0 || num_neighbours == 1
            continue
        end

        % calculate number of neighbours of ii that are connected by an edge
        % for edge e_jj,kk to be counted we require
        % ii and jj to be neighbours or A_ii,jj = 1
        % ii and kk to be neighbours or A_kk,ii = 1
        % jj and kk to be neighbours or A_jj,kk = 1
        for jj = 1:num_cluster
            for kk = 1:num_cluster
                local_coeff(ii) = local_coeff(ii) + mydistmat_cluster(ii,jj)*mydistmat_cluster(jj,kk)*mydistmat_cluster(kk,ii);
            end
        end
        local_coeff(ii) = local_coeff(ii)/(num_neighbours*(num_neighbours-1));
    end
    % store vector of local clustering coefficients
    cluster_coeff{jj_cluster} = local_coeff;
end


%% CLUSTER SIZE TIME SERIES
figure(13)
clf
betterplots
colormap darkjet
% plot average paths of clusters
subplot(2,2,[1,2]);
pcolor(xisw,zisw,uisw),shading interp
hold on
plot(meanx{1},meanz{1},'k.','MarkerSize',12) % largest cluster
plot(meanx{2},meanz{2},'w.','MarkerSize',12) % second largest
plot(meanx{3},meanz{3},'m.','MarkerSize',12) % third largest
title('Mean paths of three largest clusters')
xlabel('x')
ylabel('z')
drawnow
% plot size of the clusters over time
subplot(2,2,[3,4]);
hold on
patch([0 0 max(ts) max(ts)],[0 max(cluster_size,[],"all") max(cluster_size,[],"all") 0], [0.4 0.4 0.4])
plot(ts,cluster_size(1,:),'ko-',ts,cluster_size(2,:),'wo-',ts,cluster_size(3,:),'mo-')
xlabel('t')
ylabel('Number of particles')
title('Cluster sizes over time')
xlim([0,max(ts)])
ylim([0,max(cluster_size,[],"all")])
drawnow
savefig('outputs\cluster_size.fig')


%% DEGREE AND LOCAL COEFFICIENT HISTOGRAMS
% plot histograms of particle degree and local clustering coefficient of three largest clusters at final time
figure(14)
clf
betterplots
colormap darkjet
% calculate upper limit of histograms for particle degree
[frequencies_cluster1,edges] = histcounts(cluster_degree{1});
[frequencies_cluster2,edges] = histcounts(cluster_degree{2});
[frequencies_cluster3,edges] = histcounts(cluster_degree{3});
degree_limit = max([max(frequencies_cluster1),max(frequencies_cluster2),max(frequencies_cluster3)]);
degree_limit = 10*ceil(degree_limit/10);
% calculate upper limit of histograms for local clustering coefficient
[frequencies_cluster1,edges] = histcounts(cluster_coeff{1});
[frequencies_cluster2,edges] = histcounts(cluster_coeff{2});
[frequencies_cluster3,edges] = histcounts(cluster_coeff{3});
coeff_limit = max([max(frequencies_cluster1),max(frequencies_cluster2),max(frequencies_cluster3)]);
coeff_limit = 10*ceil(coeff_limit/10);
% plots
subplot(2,3,1);
histogram(cluster_degree{1})
xlabel('Degree')
ylabel('Number of Particles')
title('Largest Cluster')
ylim([0 degree_limit])
subplot(2,3,2);
histogram(cluster_degree{2})
xlabel('Degree')
title('Second Largest Cluster')
ylim([0 degree_limit])
subplot(2,3,3);
histogram(cluster_degree{3})
xlabel('Degree')
title('Third Largest Cluster')
ylim([0 degree_limit])
subplot(2,3,4);
histogram(cluster_coeff{1})
xlabel('Local Clustering Coefficient')
ylabel('Number of Particles')
ylim([0 coeff_limit])
subplot(2,3,5);
histogram(cluster_coeff{2})
xlabel('Local Clustering Coefficient')
ylim([0 coeff_limit])
subplot(2,3,6);
histogram(cluster_coeff{3})
xlabel('Local Clustering Coefficient')
ylim([0 coeff_limit])
savefig('outputs\cluster_degandcoeff_hist.fig')
%print('degree_and_localcoeff_a','-dtiff','-r650')


%% DEGREE AND LOCAL COEFFECIENT REPRESENTED BY MARKERSIZE IN CLUSTERS
% plot final positions of three largest clusters where Markersize is
% proportional to degree and proportional to local clustering coefficient
figure(15)
clf
betterplots
colormap darkjet
% largest cluster
subplot(2,3,1); % Markersize proportional to particle degree
pcolor(xisw,zisw,uisw),shading interp
hold on
myind=cluster_ind{1};
for ii = 1:size(myind,2)
    plot(x(myind(ii)),z(myind(ii)),'o','Markersize',2*cluster_degree{1}(ii)+1,'Color',[1 1 1]*cluster_degree{1}(ii)/max(cluster_degree{1}))
end
title('Largest Cluster Degree')
ylabel('z (m)')
xlim([min(x(myind))-10 max(x(myind))+10])
ylim([min(z(myind))-1 max(z(myind))+1])
subplot(2,3,4); % Markersize proportional to local clustering coefficient
pcolor(xisw,zisw,uisw),shading interp
hold on
for ii = 1:size(myind,2)
    plot(x(myind(ii)),z(myind(ii)),'o','Markersize',10*cluster_coeff{1}(ii)+1,'Color',[1 1 1]*cluster_coeff{1}(ii))
end
title('Largest Cluster Coefficient')
ylabel('z (m)')
xlabel('x (m)')
xlim([min(x(myind))-10 max(x(myind))+10])
ylim([min(z(myind))-1 max(z(myind))+1])
% second largest
subplot(2,3,2); % Markersize proportional to particle degree
pcolor(xisw,zisw,uisw),shading interp
hold on
myind=cluster_ind{2};
for ii = 1:size(myind,2)
    plot(x(myind(ii)),z(myind(ii)),'o','Markersize',2*cluster_degree{2}(ii)+1,'Color',[1 1 1]*cluster_degree{2}(ii)/max(cluster_degree{2}))
end
title('Second Largest Cluster Degree')
xlim([min(x(myind))-10 max(x(myind))+10])
ylim([min(z(myind))-1 max(z(myind))+1])
subplot(2,3,5); % Markersize proportional to local clustering coefficient
pcolor(xisw,zisw,uisw),shading interp
hold on
for ii = 1:size(myind,2)
    plot(x(myind(ii)),z(myind(ii)),'o','Markersize',10*cluster_coeff{2}(ii)+1,'Color',[1 1 1]*cluster_coeff{2}(ii))
end
title('Second Largest Cluster Coefficient')
xlabel('x (m)')
xlim([min(x(myind))-10 max(x(myind))+10])
ylim([min(z(myind))-1 max(z(myind))+1])
% third largest
subplot(2,3,3); % Markersize proportional to particle degree
pcolor(xisw,zisw,uisw),shading interp
hold on
myind=cluster_ind{3};
for ii = 1:size(myind,2)
    plot(x(myind(ii)),z(myind(ii)),'o','Markersize',2*cluster_degree{3}(ii)+1,'Color',[1 1 1]*cluster_degree{3}(ii)/max(cluster_degree{3}))
end
title('Third Largest Cluster Degree')
xlim([min(x(myind))-10 max(x(myind))+10])
ylim([min(z(myind))-1 max(z(myind))+1])
subplot(2,3,6); % Markersize proportional to local clustering coefficient
pcolor(xisw,zisw,uisw),shading interp
hold on
for ii = 1:size(myind,2)
    plot(x(myind(ii)),z(myind(ii)),'o','Markersize',10*cluster_coeff{3}(ii)+1,'Color',[1 1 1]*cluster_coeff{3}(ii))
end
title('Third Largest Cluster Coefficient')
xlabel('x (m)')
xlim([min(x(myind))-10 max(x(myind))+10])
ylim([min(z(myind))-1 max(z(myind))+1])
drawnow
savefig('outputs\cluster_degandcoeff_finalpos.fig')
%print('degree_and_localcoeff_b','-dtiff','-r650')


%% PARTICLE INTERACTION GRAPHS
figure(11)
clf
betterplots
subplot(2,2,1); % all subgraphs
plot(G,'NodeColor','k','MarkerSize',1)
% p=plot(G,'XData',x,'YData',z);
% p.Marker='s';
% p.NodeColor='r';
title('All Particle Interaction Graphs')
subplot(2,2,2); % largest cluster
G_cluster = graph(cluster_adjmat{1});
p = plot(G_cluster,'NodeColor','k','MarkerSize',1);
p.EdgeCData = biconncomp(G_cluster);
title('Biconnected Components of Largest Cluster')
subplot(2,2,3); % second largest cluster
G_cluster = graph(cluster_adjmat{2});
p = plot(G_cluster,'NodeColor','k','MarkerSize',1);
p.EdgeCData = biconncomp(G_cluster);
title('Biconnected Components of Second Largest Cluster')
subplot(2,2,4); % third largest cluster
G_cluster = graph(cluster_adjmat{3});
p = plot(G_cluster,'NodeColor','k','MarkerSize',1);
p.EdgeCData = biconncomp(G_cluster);
title('Biconnected Components of Third Largest Cluster')
savefig('outputs\graph.fig')
%print('graphs','-dtiff','-r650')


%% SAVE DATA TO FILE
save("outputs\plotting_data.mat","ts","mydegnowmax","xs","zs","mydegnow", ...
    "mydegnowinst","G","x","z","aa","ai","meanx","meanz","cluster_size","cluster_degree","cluster_adjmat","cluster_coeff")


%% ARCHIVED FIGURES

% figure(16)
% clf
% betterplots
% colormap darkjet
% for jj_cluster = 1:3
%     % plot bioconnected graph components
%     subplot(2,3,jj_cluster);
%     G_cluster = graph(cluster_adjmat{jj_cluster});
%     p = plot(G_cluster);
%     p.EdgeCData = biconncomp(G_cluster);
% 
%     % % plot condensed graph
%     % subplot(2,3,jj_cluster+3);
%     % G_cluster_directed = digraph(cluster_adjmat{jj_cluster});
%     % C = condensation(G_cluster_directed);
%     % plot(C,'b')
% 
%     % plot block-cut tree graph
%     subplot(2,3,jj_cluster+3);
%     tree = bctree(G_cluster);
%     plot(tree);
% end
% savefig('outputs\cluster_graphs.fig')


% The degree is the sum of all interactions in a given row so plot that
% mydeg=sum(mydistmat,2);
% figure(12)
% clf
% betterplots
% plot3(x,mydeg,z,'bo')
% grid on
% xlabel('x')
% zlabel('z')
% ylabel('degree')

% save data to file
%save("plotting_data.mat","ts","mydegnowmax","xs","zs","mydegnow","G","x","z","aa","ai","mydeg")