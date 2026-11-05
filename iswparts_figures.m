%% Load Data
clear all,close all
load('outputs\NO current\square initial condition\mysig = 0.01\Oct 24 2023\plotting_data_modified.mat')
nou_data.ts = ts;
nou_data.mydegnowmax = mydegnowmax;
nou_data.xs = xs;
nou_data.zs = zs;
nou_data.mydegnow = mydegnow;
nou_data.mydegnowinst = mydegnowinst;
nou_data.G = G;
nou_data.x = x;
nou_data.z = z;
nou_data.aa = aa;
nou_data.ai = ai;
nou_data.meanx = meanx;
nou_data.meanz = meanz;
nou_data.cluster_size = cluster_size;
nou_data.cluster_degree = cluster_degree;
nou_data.cluster_adjmat = cluster_adjmat;
nou_data.cluster_coeff = cluster_coeff;
load isw_infonou_hires.mat
nou_data.uisw = uisw;
nou_data.wisw = wisw;
nou_data.xisw = xisw;
nou_data.zisw = zisw;
load('outputs\top negative current\square initial condition\mysig = 0.01\Oct 24 2023\plotting_data_modified.mat')
topneg_data.ts = ts;
topneg_data.mydegnowmax = mydegnowmax;
topneg_data.xs = xs;
topneg_data.zs = zs;
topneg_data.mydegnow = mydegnow;
topneg_data.mydegnowinst = mydegnowinst;
topneg_data.G = G;
topneg_data.x = x;
topneg_data.z = z;
topneg_data.aa = aa;
topneg_data.ai = ai;
topneg_data.meanx = meanx;
topneg_data.meanz = meanz;
topneg_data.cluster_size = cluster_size;
topneg_data.cluster_degree = cluster_degree;
topneg_data.cluster_adjmat = cluster_adjmat;
topneg_data.cluster_coeff = cluster_coeff;
load isw_infonegu.mat
topneg_data.uisw = uisw;
topneg_data.wisw = wisw;
topneg_data.xisw = xisw;
topneg_data.zisw = zisw;
clearvars ts mydegnowmax xs zs mydegnow mydegnowinst G x z aa ai meanx meanz 
clearvars cluster_size cluster_degree cluster_adjmat cluster_coeff
clearvars uisw wisw xisw zisw c cisw
degree_max = max([nou_data.cluster_degree{1}, nou_data.cluster_degree{2}, nou_data.cluster_degree{3}, ...
    topneg_data.cluster_degree{1}, topneg_data.cluster_degree{2}, topneg_data.cluster_degree{3}])

H=20; clw=sqrt(9.81*0.002*0.25*H*0.75);


%% Wave-Induced Currents and Particle Transport
figure(3)
clf
betterplots
colormap darkjet
pcolor(nou_data.xisw,nou_data.zisw,nou_data.uisw/clw),shading interp
cmocean('balance')
hold on
plot(nou_data.xs(:,1:5:41),nou_data.zs(:,1:5:41),'k.')
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
%title('Horizontal Wave-Induced Velocities for ISW')
colorbar
clim([-1 1])
text(450,-2,'({\it a})','fontname','Times New Roman','FontSize',42,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',42,'FontWeight','normal')
%print('wave_induced_velocities_a','-dtiff','-r650')

figure(5)
clf
betterplots
colormap darkjet
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw/clw),shading interp
cmocean('balance')
hold on
plot(topneg_data.xs(:,1:5:41),topneg_data.zs(:,1:5:41),'k.')
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
%title('Horizontal Wave-Induced Velocities for ISW with Top-Trapped Negative Current')
colorbar
clim([-1 1])
text(450,-2,'({\it b})','fontname','Times New Roman','FontSize',42,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',42,'FontWeight','normal')
%print('wave_induced_velocities_c','-dtiff','-r650')


%% Three Largest Connected Clusters
figure(6)
clf
betterplots
colormap darkjet
pcolor(nou_data.xisw,nou_data.zisw,nou_data.uisw/clw),shading interp
cmocean('balance')
hold on
x0=nou_data.xs(:,1);
z0=nou_data.zs(:,1);
[bin binsize]=conncomp(nou_data.G);
myind=find(bin==nou_data.ai(1));
plot(nou_data.x(myind),nou_data.z(myind),'k.','MarkerSize',12) % connected graph of largest degree
plot(x0(myind),z0(myind),'ko','MarkerSize',12)
myind=find(bin==nou_data.ai(2));
plot(nou_data.x(myind),nou_data.z(myind),'w.','MarkerSize',12) % second largest
plot(x0(myind),z0(myind),'wo','MarkerSize',12)
myind=find(bin==nou_data.ai(3));
plot(nou_data.x(myind),nou_data.z(myind),'m.','MarkerSize',12) % third largest
plot(x0(myind),z0(myind),'mo','MarkerSize',12)
%title('Initial & Final Positions of Largest Connected Clusters for ISW')
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
text(450,-2,'({\it a})','fontname','Times New Roman','FontSize',42,'FontWeight','normal')
drawnow
colorbar
clim([-1 1])
set(gca,'fontname','Times New Roman','FontSize',42,'FontWeight','normal')
%print('three_largest_clusters_a','-dtiff','-r650')

figure(8)
clf
betterplots
colormap darkjet
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw/clw),shading interp
cmocean('balance')
hold on
x0=topneg_data.xs(:,1);
z0=topneg_data.zs(:,1);
[bin binsize]=conncomp(topneg_data.G);
myind=find(bin==topneg_data.ai(1));
plot(topneg_data.x(myind),topneg_data.z(myind),'k.','MarkerSize',12) % connected graph of largest degree
plot(x0(myind),z0(myind),'ko','MarkerSize',12)
myind=find(bin==topneg_data.ai(2));
plot(topneg_data.x(myind),topneg_data.z(myind),'w.','MarkerSize',12) % second largest
plot(x0(myind),z0(myind),'wo','MarkerSize',12)
myind=find(bin==topneg_data.ai(3));
plot(topneg_data.x(myind),topneg_data.z(myind),'m.','MarkerSize',12) % third largest
plot(x0(myind),z0(myind),'mo','MarkerSize',12)
%title('Initial & Final Positions of Largest Connected Clusters for ISW with Top-Trapped Negative Current')
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
text(450,-2,'({\it b})','fontname','Times New Roman','FontSize',42,'FontWeight','normal')
drawnow
colorbar
clim([-1 1])
set(gca,'fontname','Times New Roman','FontSize',42,'FontWeight','normal')
%print('three_largest_clusters_c','-dtiff','-r650')


%% Particle Interaction Graphs
% figure(9)
% clf
% betterplots
% subplot(2,2,1); % all subgraphs
% plot(nou_data.G,'NodeColor','k','MarkerSize',1)
% title('All Particle Interaction Graphs for ISW')
% subplot(2,2,2); % largest cluster
% G_cluster = graph(nou_data.cluster_adjmat{1});
% p = plot(G_cluster,'NodeColor','k','MarkerSize',1);
% p.EdgeCData = biconncomp(G_cluster);
% title('Biconnected Components of Largest Cluster')
% subplot(2,2,3); % second largest cluster
% G_cluster = graph(nou_data.cluster_adjmat{2});
% p = plot(G_cluster,'NodeColor','k','MarkerSize',1);
% p.EdgeCData = biconncomp(G_cluster);
% title('Biconnected Components of Second Largest Cluster')
% subplot(2,2,4); % third largest cluster
% G_cluster = graph(nou_data.cluster_adjmat{3});
% p = plot(G_cluster,'NodeColor','k','MarkerSize',1);
% p.EdgeCData = biconncomp(G_cluster);
% title('Biconnected Components of Third Largest Cluster')
% %print('graphs_a','-dtiff','-r650')

% figure(10)
% clf
% betterplots
% subplot(2,2,1); % all subgraphs
% plot(topneg_data.G,'NodeColor','k','MarkerSize',1)
% %title('All Particle Interaction Graphs for ISW with Top-Trapped Negative Current')
% text(90,95,'({\it a})','fontname','Times New Roman','FontSize',23,'FontWeight','normal')
% subplot(2,2,2); % largest cluster
% G_cluster = graph(topneg_data.cluster_adjmat{1});
% p = plot(G_cluster,'NodeColor','k','MarkerSize',1);
% p.EdgeCData = biconncomp(G_cluster);
% %title('Biconnected Components of Largest Cluster')
% text(2.7,6.25,'({\it b})','fontname','Times New Roman','FontSize',23,'FontWeight','normal')
% subplot(2,2,3); % second largest cluster
% G_cluster = graph(topneg_data.cluster_adjmat{2});
% p = plot(G_cluster,'NodeColor','k','MarkerSize',1);
% p.EdgeCData = biconncomp(G_cluster);
% %title('Biconnected Components of Second Largest Cluster')
% text(3.25,4.7,'({\it c})','fontname','Times New Roman','FontSize',23,'FontWeight','normal')
% subplot(2,2,4); % third largest cluster
% G_cluster = graph(topneg_data.cluster_adjmat{3});
% p = plot(G_cluster,'NodeColor','k','MarkerSize',1);
% p.EdgeCData = biconncomp(G_cluster);
% %title('Biconnected Components of Third Largest Cluster')
% text(4.35,3.55,'({\it d})','fontname','Times New Roman','FontSize',23,'FontWeight','normal')
% %print('graphs_c','-dtiff','-r650')

figure(11)
clf
betterplots
axes('OuterPosition',[0 0.5 0.5 0.5]) % plot all subgraphs
plot(topneg_data.G,'NodeColor','k','MarkerSize',1)
text(100,95,'({\it a})','fontname','Times New Roman','FontSize',23,'FontWeight','normal')

axes('OuterPosition',[0.5 0 0.5 1]) % plot largest connected graph, highlight biconnected components
G_cluster = graph(topneg_data.cluster_adjmat{1});
p = plot(G_cluster,'NodeColor','k','EdgeColor','k','MarkerSize',1);
%p.EdgeCData = biconncomp(G_cluster);
text(2.7,6.85,'({\it b})','fontname','Times New Roman','FontSize',23,'FontWeight','normal')

% find largest biconnected components (10 or more nodes)
biconncell = biconncomp(G_cluster,'OutputForm','cell');
biconn_size = cellfun("length",biconncell);
count_biconn = sum(biconn_size>=10);
colour_biconn = ['r','g','b','c','m','y'];
[biconn6_size,biconn6_ind] = maxk(biconn_size,count_biconn);
for ii = 1:count_biconn % highlight these components
    highlight(p,biconncell{biconn6_ind(ii)},'NodeColor',colour_biconn(ii),'MarkerSize',2)
end

axes('OuterPosition',[0 0 0.5 0.5]) % plot 6 largest biconnected components
combined_ind = [];
for ii = 1:count_biconn
    combined_ind = [combined_ind, biconncell{biconn6_ind(ii)}];
end
p = plot(subgraph(G_cluster,combined_ind),'EdgeColor','k');
p.NodeLabel = {};
low_ind = 1;
high_ind = 0;
for ii = 1:count_biconn
    high_ind = high_ind + biconn6_size(ii);
    highlight(p,low_ind:high_ind,'NodeColor',colour_biconn(ii))
    low_ind = low_ind + biconn6_size(ii);
end
text(22,21,'({\it c})','fontname','Times New Roman','FontSize',23,'FontWeight','normal')


%% Degree and Local Coefficient Histograms
h = figure(12);
clf
betterplots
colormap darkjet
% calculate upper limit of histograms for particle degree
[frequencies_cluster1,~] = histcounts(nou_data.cluster_degree{1});
[frequencies_cluster2,~] = histcounts(nou_data.cluster_degree{2});
[frequencies_cluster3,~] = histcounts(nou_data.cluster_degree{3});
degree_limit = max([max(frequencies_cluster1),max(frequencies_cluster2),max(frequencies_cluster3)]);
[frequencies_cluster1,~] = histcounts(topneg_data.cluster_degree{1});
[frequencies_cluster2,~] = histcounts(topneg_data.cluster_degree{2});
[frequencies_cluster3,~] = histcounts(topneg_data.cluster_degree{3});
degree_limit = max([degree_limit,max(frequencies_cluster1),max(frequencies_cluster2),max(frequencies_cluster3)]);
degree_limit = 10*ceil(degree_limit/10);
% calculate upper limit of histograms for local clustering coefficient
[frequencies_cluster1,~] = histcounts(nou_data.cluster_coeff{1});
[frequencies_cluster2,~] = histcounts(nou_data.cluster_coeff{2});
[frequencies_cluster3,~] = histcounts(nou_data.cluster_coeff{3});
coeff_limit = max([max(frequencies_cluster1),max(frequencies_cluster2),max(frequencies_cluster3)]);
[frequencies_cluster1,~] = histcounts(topneg_data.cluster_coeff{1});
[frequencies_cluster2,~] = histcounts(topneg_data.cluster_coeff{2});
[frequencies_cluster3,~] = histcounts(topneg_data.cluster_coeff{3});
coeff_limit = max([coeff_limit,max(frequencies_cluster1),max(frequencies_cluster2),max(frequencies_cluster3)]);
coeff_limit = 10*ceil(coeff_limit/10);

% plots
%subplot(2,3,1);
axes('OuterPosition',[0 0.5 0.33 0.5])
hold on
histogram(nou_data.cluster_degree{1},'FaceColor','b')
histogram(topneg_data.cluster_degree{1},'FaceColor',[0.7 0.7 0.7])
xlabel('$d_G(v_i)$','Interpreter','latex')
ylabel('Frequency')
% title('Largest Cluster')
ylim([0 degree_limit])
text(8.5,120,'({\it a})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

%subplot(2,3,2);
axes('OuterPosition',[0.33 0.5 0.33 0.5])
hold on
histogram(nou_data.cluster_degree{2},'FaceColor','b')
histogram(topneg_data.cluster_degree{2},'FaceColor',[0.7 0.7 0.7])
xlabel('$d_G(v_i)$','Interpreter','latex')
% title('Second Largest Cluster')
ylim([0 degree_limit])
text(6.75,120,'({\it b})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

%subplot(2,3,3);
axes('OuterPosition',[0.66 0.5 0.33 0.5])
hold on
histogram(nou_data.cluster_degree{3},'FaceColor','b')
histogram(topneg_data.cluster_degree{3},'FaceColor',[0.7 0.7 0.7])
xlabel('$d_G(v_i)$','Interpreter','latex')
% title('Third Largest Cluster')
ylim([0 degree_limit])
text(7.5,120,'({\it c})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

%subplot(2,3,4);
axes('OuterPosition',[0 0 0.33 0.5])
hold on
histogram(nou_data.cluster_coeff{1},'FaceColor','b')
histogram(topneg_data.cluster_coeff{1},'FaceColor',[0.7,0.7,0.7])
xlabel('$C_G(v_i)$','Interpreter','latex')
ylabel('Frequency')
ylim([0 coeff_limit])
text(0.9,295,'({\it d})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')
% inset
axes('OuterPosition',[5*0.33/16 5*0.5/12 5*0.33/8 5*0.5/12])
box on
hold on
histogram(nou_data.cluster_coeff{1},'FaceColor','b')
histogram(topneg_data.cluster_coeff{1},'FaceColor',[0.7,0.7,0.7])
ylim([0 50])
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

%subplot(2,3,5);
axes('OuterPosition',[0.33 0 0.33 0.5])
hold on
histogram(nou_data.cluster_coeff{2},'FaceColor','b')
histogram(topneg_data.cluster_coeff{2},'FaceColor',[0.7,0.7,0.7])
xlabel('$C_G(v_i)$','Interpreter','latex')
ylim([0 coeff_limit])
text(0.9,295,'({\it e})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')
% inset
axes('OuterPosition',[0.33+5*0.33/16 5*0.5/12 5*0.33/8 5*0.5/12])
box on
hold on
histogram(nou_data.cluster_coeff{2},'FaceColor','b')
histogram(topneg_data.cluster_coeff{2},'FaceColor',[0.7,0.7,0.7])
ylim([0 50])
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

%subplot(2,3,6);
axes('OuterPosition',[0.66 0 0.33 0.5])
hold on
histogram(nou_data.cluster_coeff{3},'FaceColor','b')
histogram(topneg_data.cluster_coeff{3},'FaceColor',[0.7,0.7,0.7])
xlabel('$C_G(v_i)$','Interpreter','latex')
ylim([0 coeff_limit])
text(0.9,295,'({\it f})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')
% inset
axes('OuterPosition',[0.66+5*0.33/16 5*0.5/12 5*0.33/8 5*0.5/12])
box on
hold on
histogram(nou_data.cluster_coeff{3},'FaceColor','b')
histogram(topneg_data.cluster_coeff{3},'FaceColor',[0.7,0.7,0.7])
ylim([0 50])
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

%print('degree_and_localcoeff_hist','-dtiff','-r650')


%% Degree and Local Coefficient Table
%[NoU_Degree_Cluster1,~] =  histcounts(nou_data.cluster_degree{1});
%[NoU_Degree_Cluster2,~] =  histcounts(nou_data.cluster_degree{2});
%[NoU_Degree_Cluster3,~] =  histcounts(nou_data.cluster_degree{3});
%if length(NoU_Degree_Cluster1) < degree_max
%    NoU_Degree_Cluster1 = [NoU_Degree_Cluster1, zeros(1,degree_max-length(NoU_Degree_Cluster1))];
%end
%if length(NoU_Degree_Cluster2) < degree_max
%    NoU_Degree_Cluster2 = [NoU_Degree_Cluster2, zeros(1,degree_max-length(NoU_Degree_Cluster2))];
%end
%if length(NoU_Degree_Cluster3) < degree_max
%    NoU_Degree_Cluster3 = [NoU_Degree_Cluster3, zeros(1,degree_max-length(NoU_Degree_Cluster3))];
%end
%[TopNeg_Degree_Cluster1,~] =  histcounts(topneg_data.cluster_degree{1});
%[TopNeg_Degree_Cluster2,~] =  histcounts(topneg_data.cluster_degree{2});
%[TopNeg_Degree_Cluster3,~] =  histcounts(topneg_data.cluster_degree{3});
%if length(TopNeg_Degree_Cluster1) < degree_max
%    TopNeg_Degree_Cluster1 = [TopNeg_Degree_Cluster1, zeros(1,degree_max-length(TopNeg_Degree_Cluster1))];
%end
%if length(TopNeg_Degree_Cluster2) < degree_max
%    TopNeg_Degree_Cluster2 = [TopNeg_Degree_Cluster2, zeros(1,degree_max-length(TopNeg_Degree_Cluster2))];
%end
%if length(TopNeg_Degree_Cluster3) < degree_max
%    TopNeg_Degree_Cluster3 = [TopNeg_Degree_Cluster3, zeros(1,degree_max-length(TopNeg_Degree_Cluster3))];
%end
%DegreeLabels = linspace(1,11,11)';
%DegreeValues = table(DegreeLabels,NoU_Degree_Cluster1',NoU_Degree_Cluster2',NoU_Degree_Cluster3', ...
%    TopNeg_Degree_Cluster1',TopNeg_Degree_Cluster2',TopNeg_Degree_Cluster3')
%
%[NoU_Coeff_Cluster1,coeff_edges] =  histcounts(nou_data.cluster_coeff{1});
%[NoU_Coeff_Cluster2,~] =  histcounts(nou_data.cluster_coeff{2});
%[NoU_Coeff_Cluster3,~] =  histcounts(nou_data.cluster_coeff{3});
%[TopNeg_Coeff_Cluster1,~] =  histcounts(topneg_data.cluster_coeff{1});
%[TopNeg_Coeff_Cluster2,~] =  histcounts(topneg_data.cluster_coeff{2});
%[TopNeg_Coeff_Cluster3,~] =  histcounts(topneg_data.cluster_coeff{3});
%CoeffLabels = [];
%for ii = 1:length(coeff_edges)-1
%    if ii == length(coeff_edges)-1
%        CoeffLabels = [CoeffLabels; "["+num2str(coeff_edges(ii))+","+num2str(coeff_edges(ii+1))+"]"];
%    else
%        CoeffLabels = [CoeffLabels; "["+num2str(coeff_edges(ii))+","+num2str(coeff_edges(ii+1))+")"];
%    end
%end
%CoeffValues = table(CoeffLabels,NoU_Coeff_Cluster1',NoU_Coeff_Cluster2',NoU_Coeff_Cluster3', ...
%    TopNeg_Coeff_Cluster1',TopNeg_Coeff_Cluster2',TopNeg_Coeff_Cluster3')


%% Degree and Local Coefficient Proportional to MarkerSize in Clusters
figure(13)
clf
betterplots
colormap darkjet
% largest cluster in plain ISW
axes('OuterPosition',[0 0.1+0.9/2 0.95/3 0.9/2]) % Markersize proportional to particle degree
pcolor(nou_data.xisw,nou_data.zisw,nou_data.uisw/clw),shading interp
cmocean('balance')
clim([-1 1])
hold on
[bin,~]=conncomp(nou_data.G);
myind=find(bin==nou_data.ai(1));
for ii = 1:size(myind,2)
    plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'o','Markersize',2*nou_data.cluster_degree{1}(ii)+1,'Color',[1 1 1]*nou_data.cluster_degree{1}(ii)/degree_max)
    if nou_data.cluster_coeff{1}(ii) > 0
        plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'mx')
    end
end
% title('Largest Cluster Degree')
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
xlim([min(nou_data.x(myind))-10 max(nou_data.x(myind))+10])
ylim([min(nou_data.z(myind))-1 max(nou_data.z(myind))+1])
text(215,-2.5,'({\it a})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

% second largest in plain ISW
axes('OuterPosition',[0.95/3 0.1+0.9/2 0.95/3 0.9/2]) % Markersize proportional to particle degree
pcolor(nou_data.xisw,nou_data.zisw,nou_data.uisw/clw),shading interp
cmocean('balance')
clim([-1 1])
hold on
myind=find(bin==nou_data.ai(2));
for ii = 1:size(myind,2)
    plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'o','Markersize',2*nou_data.cluster_degree{2}(ii)+1,'Color',[1 1 1]*nou_data.cluster_degree{2}(ii)/degree_max)
    if nou_data.cluster_coeff{2}(ii) > 0
        plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'mx')
    end
end
% title('Second Largest Cluster Degree')
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
xlim([min(nou_data.x(myind))-10 max(nou_data.x(myind))+10])
ylim([min(nou_data.z(myind))-1 max(nou_data.z(myind))+1])
text(240,-1.2,'({\it b})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

% third largest in plain ISW
axes('OuterPosition',[2*0.95/3 0.1+0.9/2 0.95/3 0.9/2]) % Markersize proportional to particle degree
pcolor(nou_data.xisw,nou_data.zisw,nou_data.uisw/clw),shading interp
cmocean('balance')
clim([-1 1])
hold on
myind=find(bin==nou_data.ai(3));
for ii = 1:size(myind,2)
    plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'o','Markersize',2*nou_data.cluster_degree{3}(ii)+1,'Color',[1 1 1]*nou_data.cluster_degree{3}(ii)/degree_max)
    if nou_data.cluster_coeff{3}(ii) > 0
        plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'mx')
    end
end
% title('Third Largest Cluster Degree')
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
xlim([min(nou_data.x(myind))-10 max(nou_data.x(myind))+10])
ylim([min(nou_data.z(myind))-1 max(nou_data.z(myind))+1])
text(169,-4,'({\it c})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

% largest cluster in top neg current
axes('OuterPosition',[0 0.1 0.95/3 0.9/2]) % Markersize proportional to particle degree
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw/clw),shading interp
cmocean('balance')
clim([-1 1])
hold on
[bin,~]=conncomp(topneg_data.G);
myind=find(bin==topneg_data.ai(1));
for ii = 1:size(myind,2)
    plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'o','Markersize',2*topneg_data.cluster_degree{1}(ii)+1,'Color',[1 1 1]*topneg_data.cluster_degree{1}(ii)/degree_max)
    if topneg_data.cluster_coeff{1}(ii) > 0
        plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'mx')
    end
end
% title('Largest Cluster Degree')
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
xlim([min(topneg_data.x(myind))-10 max(topneg_data.x(myind))+10])
ylim([min(topneg_data.z(myind))-1 max(topneg_data.z(myind))+1])
text(202.5,-10.75,'({\it d})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

% second largest in top neg current
axes('OuterPosition',[0.95/3 0.1 0.95/3 0.9/2]) % Markersize proportional to particle degree
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw/clw),shading interp
cmocean('balance')
clim([-1 1])
hold on
myind=find(bin==topneg_data.ai(2));
for ii = 1:size(myind,2)
    plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'o','Markersize',2*topneg_data.cluster_degree{2}(ii)+1,'Color',[1 1 1]*topneg_data.cluster_degree{2}(ii)/degree_max)
    if topneg_data.cluster_coeff{2}(ii) > 0
        plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'mx')
    end
end
% title('Second Largest Cluster Degree')
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
xlim([min(topneg_data.x(myind))-10 max(topneg_data.x(myind))+10])
ylim([min(topneg_data.z(myind))-1 max(topneg_data.z(myind))+1])
text(219.25,-4.15,'({\it e})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

% third largest in top neg current
axes('OuterPosition',[2*0.95/3 0.1 0.95/3 0.9/2]) % Markersize proportional to particle degree
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw/clw),shading interp
cmocean('balance')
clim([-1 1])
hold on
myind=find(bin==topneg_data.ai(3));
for ii = 1:size(myind,2)
    plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'o','Markersize',2*topneg_data.cluster_degree{3}(ii)+1,'Color',[1 1 1]*topneg_data.cluster_degree{3}(ii)/degree_max)
    if topneg_data.cluster_coeff{3}(ii) > 0
        plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'mx')
    end
end
% title('Third Largest Cluster Degree')
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
xlim([min(topneg_data.x(myind))-10 max(topneg_data.x(myind))+10])
ylim([min(topneg_data.z(myind))-1 max(topneg_data.z(myind))+1])
text(239.25,-2.3,'({\it f})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

% colorbar on right side of figure
cb = colorbar;
clim([-1 1])
set(cb,'Position',[0.95 0.185 0.02 0.78])
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

% legend at bottom of figure
axes('OuterPosition',[0 0.02 0.9 0.075],'NextPlot','add','DataAspectRatio',[1 1 1],'XLim',[0 degree_max],'YLim',[0 0.0001])
sample_points = linspace(0,degree_max,degree_max+1);
hold on
for ii = 1:length(sample_points)
    plot(sample_points(ii), [0],'o','Markersize',2*sample_points(ii)+1,'Color',[1 1 1]*sample_points(ii)/degree_max)
end
xlabel('$d_G(v_i)$','Interpreter','latex')
xticks(1:degree_max)
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

%print('degree_and_localcoeff_markersize','-dtiff','-r650')


%% Graphical Abstract
figure(101)
clf
betterplots
colormap darkjet
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw/clw),shading interp
cmocean('balance')
clim([-1 1])
hold on
plot(topneg_data.xs(:,1:5:41),topneg_data.zs(:,1:5:41),'k.')
pbaspect([1.2 1 1])
%print('GraphicalAbstract','-djpeg','-r300')