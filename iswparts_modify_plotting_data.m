%% DESCRIPTION: %%
% for an existing file 'plotting_data.mat', this calculates the adjacency
% matrix and all subsequent data according to the changes implemented on
% June 19, 2024 and saves to a new file 'plotting_data_modified.mat'. This
% code does not produce any figures as you can use the new data file as
% input for 'iswparts_figures.m'
clear all,close all
load ('outputs\top negative current\square initial condition\mysig = 0.01\Oct 24 2023\plotting_data.mat')

%% Particles for ISW
H=20; % total depth
L=500; %domain length

%% Particle parameters %%
numsteps=200; % number of steps integrated over between outs
numouts=40; %default numouts = 1000, saved outputs
dt=0.1; % timestep

g = 9.81; %gravity constant
settling= -0.001; % settling velocity


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
mysig = 0.01;

eps=0.5*min([dz0, dx0]); % sets a tolerance distance for particle interactions
mydistmat=zeros(numens); % initialize adjacency matrix
mydistmatcum=mydistmat; % initialize adjacency matrix with edge weightings
zero_diag = ones(numens) - eye(numens); % used to ensure A_ii = 0 for all i
%mydistmat_changes=zeros(numens,numens,numouts);

%% The double loop over outputs and timesteps
for ii=1:numouts
    % extract position data
    x = xs(:,ii);
    z = zs(:,ii);
    t = ts(ii);

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


%% PLOTTING CALCULATIONS
% MAX DEGREE TIME SERIES
[mydegnowmax, mydegnowmaxi]=max(mydegnow,[],2);
    % mydegnowmax is an array of the maximum degree at each timestep
    % mydegnowmaxi is an array of the indices corresponding to mydegnowmax
[mydeginstmax, mydeginstmaxi]=max(mydegnowinst,[],2); % same for the instantaneous degree

% POSITIONS OF THREE LARGEST CONNECTED CLUSTERS
G=graph(mydistmat);
[bin binsize]=conncomp(G); % finds all connected graphs in G, stores degree of all graphs in binsize
[aa ai]=maxk(binsize,3);
    % aa stores 3 largest degree values of connected graphs
    % ai stores corresponding indices

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

%% SAVE DATA TO FILE
save("outputs\plotting_data_modified.mat","ts","mydegnowmax","xs","zs","mydegnow", ...
    "mydegnowinst","G","x","z","aa","ai","meanx","meanz","cluster_size","cluster_degree","cluster_adjmat","cluster_coeff")