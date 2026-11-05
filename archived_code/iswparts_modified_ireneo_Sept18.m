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
load isw_bot
sz=size(uisw);Nz=sz(1);Nx=sz(2);
zmin=min(zisw(:));

%% Inline Functions: flow from the ISW data
% one trick I use is to zero out the velocity if we are below the bottom
u=@(xq,zq,t,x,z,uvals) (zq>zmin).*interp2(x,z,uvals,xq,zq,'cubic',0)+0;  
w=@(xq,zq,t,x,z,wvals) (zq>zmin).*interp2(x,z,wvals,xq,zq,'cubic',0)+0;
% x,z are the sample coordinates
% uvals/wvals are the function value at the sample points
% xq,zq are the coordinates to be interpolated at

%% initial conditions for a square
% from x = 350 to 360, z = -1 to -18
x=[];z=[];
for zstart=-1:-2:-18
    xg1=linspace(350,360,20);
    zg1=linspace(zstart-1,zstart,20);
    [dumxx dumzz]=meshgrid(xg1,zg1);
    x=[x;dumxx(:)];
    z=[z;dumzz(:)];
    numens=length(x);
end
    dz0=abs(zg1(2)-zg1(1)); % differential change in z
    dx0=abs(xg1(2)-xg1(1)); % differential change in x
%% red noise parameters
memt=1;
rho=exp(-dt/memt); rhoc=sqrt(1-rho*rho);
mysig=0.01;
prevrand1=zeros(numens,1);
prevrand2=prevrand1;

xs=zeros(numens,numouts+1);
zs=zeros(numens,numouts+1);
t=0;
eps=0.5*min([dz0, dx0]); % sets a tolerance distance for particle interactions
mydistmat=zeros(numens); % initialize adjacency matrix
mydistmatcum=mydistmat; % initialize adjacency matrix with edge weightings
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
        x=mod(x+dt*rhx,L);
            % updated x position
        rhz=w(x,z,t,x_interp,z_interp,w_interp)+settling.*(z>zmin)+(z>zmin).*pertnow2;
            % interpolated w velocity + settling/gravity + velocity perturbation
        z=z+dt*rhz;
            % updated z position
        z=z.*(z<0)-z.*(z>=0);
            % ensure z is negative
        t=t+dt;
            % updated time value
    end
    % forming the adjacency matrix
    distv=pdist([x,z]); % find interparticle distance between each pair
    distvind=(distv<eps); % note when particles interact/come within tolerance distance of eachother
    mydistmat_now=squareform(distvind); % turn into matrix where 1 indicates particle interaction
    mydistmat=mydistmat+mydistmat_now.*(mydistmat==0);
        % matrix stores 1 if particles interacted at any time, 0 otherwise
    mydistmatcum=mydistmatcum+mydistmat_now;
        % matrix stores how many times particles interacted
    %mydistmat_changes(:,:,ii)=mydistmat_now;
end
% One last storage
xs(:,numouts+1)=x;
zs(:,numouts+1)=z;
ts(numouts+1)=t;
%% Plots      
figure(1)
clf
betterplots
colormap darkjet
pcolor(xisw,zisw,uisw),shading interp
hold on
plot(xs(:,1:5:41),zs(:,1:5:41),'k.')
savefig('particle_transport.fig')

G=graph(mydistmat);
figure(11)
clf
betterplots
plot(G,'k')
% p=plot(G,'XData',x,'YData',z);
% p.Marker='s';
% p.NodeColor='r';
savefig('graph.fig')

% The degree is the sum of all interactions in a given row so plot that
mydeg=sum(mydistmat,2);
figure(12)
clf
betterplots
plot3(x,mydeg,z,'bo')
grid on
xlabel('x')
zlabel('z')
ylabel('degree')
savefig('degree.fig')

% save data to file
save("plotting_data.mat","xs","zs","G","x","z","mydeg")