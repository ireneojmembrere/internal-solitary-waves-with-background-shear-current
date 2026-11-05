%% Load Data
clear all,close all
load('outputs\bottom current\square initial condition\mysig = 0.01\Oct 24 2023\plotting_data.mat')
bottom_data.ts = ts;
bottom_data.mydegnowmax = mydegnowmax;
bottom_data.xs = xs;
bottom_data.zs = zs;
bottom_data.mydegnow = mydegnow;
bottom_data.mydegnowinst = mydegnowinst;
bottom_data.G = G;
bottom_data.x = x;
bottom_data.z = z;
bottom_data.aa = aa;
bottom_data.ai = ai;
bottom_data.meanx = meanx;
bottom_data.meanz = meanz;
bottom_data.cluster_size = cluster_size;
bottom_data.cluster_degree = cluster_degree;
bottom_data.cluster_adjmat = cluster_adjmat;
bottom_data.cluster_coeff = cluster_coeff;
load isw_bot.mat
bottom_data.uisw = uisw;
bottom_data.wisw = wisw;
bottom_data.xisw = xisw;
bottom_data.zisw = zisw;
load('outputs\NO current\square initial condition\mysig = 0.01\Oct 24 2023\plotting_data.mat')
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
load('outputs\top negative current\square initial condition\mysig = 0.01\Oct 24 2023\plotting_data.mat')
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


%% Legend for Particle Degree Represented by MarkerSize
figure(1)
clf
axes('NextPlot','add','DataAspectRatio',[1 1 1],'XLim',[0 11],'YLim',[0 0.0001])
sample_points = linspace(0,11,12);
hold on
for ii = 1:11
    plot(sample_points(ii), [0],'o','Markersize',2*sample_points(ii)+1,'Color',[1 1 1]*sample_points(ii)/11)
end
xlabel('Particle Degree')
%print('legend_degree','-dtiff','-r650')


%% Legend for Local Clustering Coefficient Represented by MarkerSize
figure(2)
clf
axes('NextPlot','add','DataAspectRatio',[1 1 1],'XLim',[0 1],'YLim',[0 0.00001])
sample_points = linspace(0,1,11);
hold on
for ii = 1:11
    plot(sample_points(ii), [0],'o','Markersize',10*sample_points(ii)+1,'Color',[1 1 1]*sample_points(ii))
end
xlabel('Local Clustering Coefficient')
%print('legend_localcoeff','-dtiff','-r650')


%% Wave-Induced Currents and Particle Transport
figure(3)
clf
betterplots
colormap darkjet
pcolor(nou_data.xisw,nou_data.zisw,nou_data.uisw),shading interp
hold on
plot(nou_data.xs(:,1:5:41),nou_data.zs(:,1:5:41),'k.')
xlabel('x (m)')
ylabel('z (m)')
%title('Horizontal Wave-Induced Velocities for ISW')
colorbar
text(475,-2,'(a)')
%print('wave_induced_velocities_a','-dtiff','-r650')

% figure(4)
% clf
% betterplots
% colormap darkjet
% pcolor(bottom_data.xisw,bottom_data.zisw,bottom_data.uisw),shading interp
% hold on
% plot(bottom_data.xs(:,1:5:41),bottom_data.zs(:,1:5:41),'k.')
% xlabel('x (m)')
% ylabel('z (m)')
% title('Horizontal Wave-Induced Velocities for ISW with Bottom Current')
% colorbar
% %print('wave_induced_velocities_b','-dtiff','-r650')

figure(5)
clf
betterplots
colormap darkjet
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw),shading interp
hold on
plot(topneg_data.xs(:,1:5:41),topneg_data.zs(:,1:5:41),'k.')
xlabel('x (m)')
ylabel('z (m)')
%title('Horizontal Wave-Induced Velocities for ISW with Top-Trapped Negative Current')
colorbar
text(475,-2,'(b)')
%print('wave_induced_velocities_c','-dtiff','-r650')


%% Three Largest Connected Clusters
figure(6)
clf
betterplots
colormap darkjet
pcolor(nou_data.xisw,nou_data.zisw,nou_data.uisw),shading interp
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
xlabel('x (m)')
ylabel('z (m)')
text(475,-2,'(a)')
drawnow
colorbar
%print('three_largest_clusters_a','-dtiff','-r650')

% figure(7)
% clf
% betterplots
% colormap darkjet
% pcolor(bottom_data.xisw,bottom_data.zisw,bottom_data.uisw),shading interp
% hold on
% x0=bottom_data.xs(:,1);
% z0=bottom_data.zs(:,1);
% [bin binsize]=conncomp(bottom_data.G);
% myind=find(bin==bottom_data.ai(1));
% plot(bottom_data.x(myind),bottom_data.z(myind),'k.','MarkerSize',12) % connected graph of largest degree
% plot(x0(myind),z0(myind),'ko','MarkerSize',12)
% myind=find(bin==bottom_data.ai(2));
% plot(bottom_data.x(myind),bottom_data.z(myind),'w.','MarkerSize',12) % second largest
% plot(x0(myind),z0(myind),'wo','MarkerSize',12)
% myind=find(bin==bottom_data.ai(3));
% plot(bottom_data.x(myind),bottom_data.z(myind),'m.','MarkerSize',12) % third largest
% plot(x0(myind),z0(myind),'mo','MarkerSize',12)
% title('Initial & Final Positions of Largest Connected Clusters for ISW with Bottom Current')
% xlabel('x (m)')
% ylabel('z (m)')
% drawnow
% colorbar
% %print('three_largest_clusters_b','-dtiff','-r650')

figure(8)
clf
betterplots
colormap darkjet
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw),shading interp
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
xlabel('x (m)')
ylabel('z (m)')
text(475,-2,'(b)')
drawnow
colorbar
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
% plot(bottom_data.G,'NodeColor','k','MarkerSize',1)
% title('All Particle Interaction Graphs for ISW with Bottom Current')
% subplot(2,2,2); % largest cluster
% G_cluster = graph(bottom_data.cluster_adjmat{1});
% p = plot(G_cluster,'NodeColor','k','MarkerSize',1);
% p.EdgeCData = biconncomp(G_cluster);
% title('Biconnected Components of Largest Cluster')
% subplot(2,2,3); % second largest cluster
% G_cluster = graph(bottom_data.cluster_adjmat{2});
% p = plot(G_cluster,'NodeColor','k','MarkerSize',1);
% p.EdgeCData = biconncomp(G_cluster);
% title('Biconnected Components of Second Largest Cluster')
% subplot(2,2,4); % third largest cluster
% G_cluster = graph(bottom_data.cluster_adjmat{3});
% p = plot(G_cluster,'NodeColor','k','MarkerSize',1);
% p.EdgeCData = biconncomp(G_cluster);
% title('Biconnected Components of Third Largest Cluster')
% %print('graphs_b','-dtiff','-r650')

figure(11)
clf
betterplots
subplot(2,2,1); % all subgraphs
plot(topneg_data.G,'NodeColor','k','MarkerSize',1)
%title('All Particle Interaction Graphs for ISW with Top-Trapped Negative Current')
text(95,95,'(a)')
subplot(2,2,2); % largest cluster
G_cluster = graph(topneg_data.cluster_adjmat{1});
p = plot(G_cluster,'NodeColor','k','MarkerSize',1);
p.EdgeCData = biconncomp(G_cluster);
%title('Biconnected Components of Largest Cluster')
text(3,6.25,'(b)')
subplot(2,2,3); % second largest cluster
G_cluster = graph(topneg_data.cluster_adjmat{2});
p = plot(G_cluster,'NodeColor','k','MarkerSize',1);
p.EdgeCData = biconncomp(G_cluster);
%title('Biconnected Components of Second Largest Cluster')
text(3.45,4.7,'(c)')
subplot(2,2,4); % third largest cluster
G_cluster = graph(topneg_data.cluster_adjmat{3});
p = plot(G_cluster,'NodeColor','k','MarkerSize',1);
p.EdgeCData = biconncomp(G_cluster);
%title('Biconnected Components of Third Largest Cluster')
text(4.6,3.55,'(d)')
%print('graphs_c','-dtiff','-r650')


%% Degree and Local Coefficient Histograms
figure(12)
clf
betterplots
colormap darkjet
% calculate upper limit of histograms for particle degree
[frequencies_cluster1,~] = histcounts(nou_data.cluster_degree{1});
[frequencies_cluster2,~] = histcounts(nou_data.cluster_degree{2});
[frequencies_cluster3,~] = histcounts(nou_data.cluster_degree{3});
degree_limit = max([max(frequencies_cluster1),max(frequencies_cluster2),max(frequencies_cluster3)]);
% [frequencies_cluster1,~] = histcounts(bottom_data.cluster_degree{1});
% [frequencies_cluster2,~] = histcounts(bottom_data.cluster_degree{2});
% [frequencies_cluster3,~] = histcounts(bottom_data.cluster_degree{3});
% degree_limit = max([degree_limit,max(frequencies_cluster1),max(frequencies_cluster2),max(frequencies_cluster3)]);
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
% [frequencies_cluster1,~] = histcounts(bottom_data.cluster_coeff{1});
% [frequencies_cluster2,~] = histcounts(bottom_data.cluster_coeff{2});
% [frequencies_cluster3,~] = histcounts(bottom_data.cluster_coeff{3});
% coeff_limit = max([coeff_limit,max(frequencies_cluster1),max(frequencies_cluster2),max(frequencies_cluster3)]);
[frequencies_cluster1,~] = histcounts(topneg_data.cluster_coeff{1});
[frequencies_cluster2,~] = histcounts(topneg_data.cluster_coeff{2});
[frequencies_cluster3,~] = histcounts(topneg_data.cluster_coeff{3});
coeff_limit = max([coeff_limit,max(frequencies_cluster1),max(frequencies_cluster2),max(frequencies_cluster3)]);
coeff_limit = 10*ceil(coeff_limit/10);
% plots
subplot(2,3,1);
hold on
histogram(nou_data.cluster_degree{1},'FaceColor',[0.7,0.7,0.7],'EdgeColor','r')
% histogram(bottom_data.cluster_degree{1},'FaceColor',[0.7,0.7,0.7],'EdgeColor','g')
histogram(topneg_data.cluster_degree{1},'FaceColor',[0.7,0.7,0.7],'EdgeColor','b')
% legend('ISW','Bottom Current','Top-Negative Current')
xlabel('Degree')
ylabel('Number of Particles')
% title('Largest Cluster')
ylim([0 degree_limit])
text(11,110,'(a)')

subplot(2,3,2);
hold on
histogram(nou_data.cluster_degree{2},'FaceColor',[0.7,0.7,0.7],'EdgeColor','r')
% histogram(bottom_data.cluster_degree{2},'FaceColor',[0.7,0.7,0.7],'EdgeColor','g')
histogram(topneg_data.cluster_degree{2},'FaceColor',[0.7,0.7,0.7],'EdgeColor','b')
% legend('ISW','Bottom Current','Top-Negative Current')
xlabel('Degree')
% title('Second Largest Cluster')
ylim([0 degree_limit])
text(10,110,'(b)')

subplot(2,3,3);
hold on
histogram(nou_data.cluster_degree{3},'FaceColor',[0.7,0.7,0.7],'EdgeColor','r')
% histogram(bottom_data.cluster_degree{3},'FaceColor',[0.7,0.7,0.7],'EdgeColor','g')
histogram(topneg_data.cluster_degree{3},'FaceColor',[0.7,0.7,0.7],'EdgeColor','b')
% legend('ISW','Bottom Current','Top-Negative Current')
xlabel('Degree')
% title('Third Largest Cluster')
ylim([0 degree_limit])
text(10,110,'(c)')

subplot(2,3,4);
hold on
histogram(nou_data.cluster_coeff{1},'FaceColor',[0.7,0.7,0.7],'EdgeColor','r')
% histogram(bottom_data.cluster_coeff{1},'FaceColor',[0.7,0.7,0.7],'EdgeColor','g')
histogram(topneg_data.cluster_coeff{1},'FaceColor',[0.7,0.7,0.7],'EdgeColor','b')
% legend('ISW','Bottom Current','Top-Negative Current')
xlabel('Local Clustering Coefficient')
ylabel('Number of Particles')
ylim([0 coeff_limit])
text(0.95,295,'(d)')

subplot(2,3,5);
hold on
histogram(nou_data.cluster_coeff{2},'FaceColor',[0.7,0.7,0.7],'EdgeColor','r')
% histogram(bottom_data.cluster_coeff{2},'FaceColor',[0.7,0.7,0.7],'EdgeColor','g')
histogram(topneg_data.cluster_coeff{2},'FaceColor',[0.7,0.7,0.7],'EdgeColor','b')
% legend('ISW','Bottom Current','Top-Negative Current')
xlabel('Local Clustering Coefficient')
ylim([0 coeff_limit])
text(0.95,295,'(e)')

subplot(2,3,6);
hold on
histogram(nou_data.cluster_coeff{3},'FaceColor',[0.7,0.7,0.7],'EdgeColor','r')
% histogram(bottom_data.cluster_coeff{3},'FaceColor',[0.7,0.7,0.7],'EdgeColor','g')
histogram(topneg_data.cluster_coeff{3},'FaceColor',[0.7,0.7,0.7],'EdgeColor','b')
% legend('ISW','Bottom Current','Top-Negative Current')
xlabel('Local Clustering Coefficient')
ylim([0 coeff_limit])
text(0.95,295,'(f)')
%print('degree_and_localcoeff_hist','-dtiff','-r650')


%% Degree and Local Coefficient Proportional to MarkerSize in Clusters
figure(13)
clf
betterplots
colormap darkjet
% largest cluster
subplot(2,3,1); % Markersize proportional to particle degree
pcolor(nou_data.xisw,nou_data.zisw,nou_data.uisw),shading interp
hold on
[bin,~]=conncomp(nou_data.G);
myind=find(bin==nou_data.ai(1));
for ii = 1:size(myind,2)
    plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'o','Markersize',2*nou_data.cluster_degree{1}(ii)+1,'Color',[1 1 1]*nou_data.cluster_degree{1}(ii)/11)
end
title('Largest Cluster Degree')
ylabel('z (m)')
xlim([min(nou_data.x(myind))-10 max(nou_data.x(myind))+10])
ylim([min(nou_data.z(myind))-1 max(nou_data.z(myind))+1])
subplot(2,3,4); % Markersize proportional to local clustering coefficient
pcolor(nou_data.xisw,nou_data.zisw,nou_data.uisw),shading interp
hold on
for ii = 1:size(myind,2)
    plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'o','Markersize',10*nou_data.cluster_coeff{1}(ii)+1,'Color',[1 1 1]*nou_data.cluster_coeff{1}(ii))
end
title('Largest Cluster Coefficient')
ylabel('z (m)')
xlabel('x (m)')
xlim([min(nou_data.x(myind))-10 max(nou_data.x(myind))+10])
ylim([min(nou_data.z(myind))-1 max(nou_data.z(myind))+1])
% second largest
subplot(2,3,2); % Markersize proportional to particle degree
pcolor(nou_data.xisw,nou_data.zisw,nou_data.uisw),shading interp
hold on
myind=find(bin==nou_data.ai(2));
for ii = 1:size(myind,2)
    plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'o','Markersize',2*nou_data.cluster_degree{2}(ii)+1,'Color',[1 1 1]*nou_data.cluster_degree{2}(ii)/11)
end
title('Second Largest Cluster Degree')
xlim([min(nou_data.x(myind))-10 max(nou_data.x(myind))+10])
ylim([min(nou_data.z(myind))-1 max(nou_data.z(myind))+1])
subplot(2,3,5); % Markersize proportional to local clustering coefficient
pcolor(nou_data.xisw,nou_data.zisw,nou_data.uisw),shading interp
hold on
for ii = 1:size(myind,2)
    plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'o','Markersize',10*nou_data.cluster_coeff{2}(ii)+1,'Color',[1 1 1]*nou_data.cluster_coeff{2}(ii))
end
title('Second Largest Cluster Coefficient')
xlabel('x (m)')
xlim([min(nou_data.x(myind))-10 max(nou_data.x(myind))+10])
ylim([min(nou_data.z(myind))-1 max(nou_data.z(myind))+1])
% third largest
subplot(2,3,3); % Markersize proportional to particle degree
pcolor(nou_data.xisw,nou_data.zisw,nou_data.uisw),shading interp
hold on
myind=find(bin==nou_data.ai(3));
for ii = 1:size(myind,2)
    plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'o','Markersize',2*nou_data.cluster_degree{3}(ii)+1,'Color',[1 1 1]*nou_data.cluster_degree{3}(ii)/11)
end
title('Third Largest Cluster Degree')
xlim([min(nou_data.x(myind))-10 max(nou_data.x(myind))+10])
ylim([min(nou_data.z(myind))-1 max(nou_data.z(myind))+1])
subplot(2,3,6); % Markersize proportional to local clustering coefficient
pcolor(nou_data.xisw,nou_data.zisw,nou_data.uisw),shading interp
hold on
for ii = 1:size(myind,2)
    plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'o','Markersize',10*nou_data.cluster_coeff{3}(ii)+1,'Color',[1 1 1]*nou_data.cluster_coeff{3}(ii))
end
title('Third Largest Cluster Coefficient')
xlabel('x (m)')
xlim([min(nou_data.x(myind))-10 max(nou_data.x(myind))+10])
ylim([min(nou_data.z(myind))-1 max(nou_data.z(myind))+1])
drawnow
%print('degree_and_localcoeff_markersize_a','-dtiff','-r650')


figure(14)
clf
betterplots
colormap darkjet
% largest cluster
subplot(2,3,1); % Markersize proportional to particle degree
pcolor(bottom_data.xisw,bottom_data.zisw,bottom_data.uisw),shading interp
hold on
[bin,~]=conncomp(bottom_data.G);
myind=find(bin==bottom_data.ai(1));
for ii = 1:size(myind,2)
    plot(bottom_data.x(myind(ii)),bottom_data.z(myind(ii)),'o','Markersize',2*bottom_data.cluster_degree{1}(ii)+1,'Color',[1 1 1]*bottom_data.cluster_degree{1}(ii)/11)
end
title('Largest Cluster Degree')
ylabel('z (m)')
xlim([min(bottom_data.x(myind))-10 max(bottom_data.x(myind))+10])
ylim([min(bottom_data.z(myind))-1 max(bottom_data.z(myind))+1])
subplot(2,3,4); % Markersize proportional to local clustering coefficient
pcolor(bottom_data.xisw,bottom_data.zisw,bottom_data.uisw),shading interp
hold on
for ii = 1:size(myind,2)
    plot(bottom_data.x(myind(ii)),bottom_data.z(myind(ii)),'o','Markersize',10*bottom_data.cluster_coeff{1}(ii)+1,'Color',[1 1 1]*bottom_data.cluster_coeff{1}(ii))
end
title('Largest Cluster Coefficient')
ylabel('z (m)')
xlabel('x (m)')
xlim([min(bottom_data.x(myind))-10 max(bottom_data.x(myind))+10])
ylim([min(bottom_data.z(myind))-1 max(bottom_data.z(myind))+1])
% second largest
subplot(2,3,2); % Markersize proportional to particle degree
pcolor(bottom_data.xisw,bottom_data.zisw,bottom_data.uisw),shading interp
hold on
myind=find(bin==bottom_data.ai(2));
for ii = 1:size(myind,2)
    plot(bottom_data.x(myind(ii)),bottom_data.z(myind(ii)),'o','Markersize',2*bottom_data.cluster_degree{2}(ii)+1,'Color',[1 1 1]*bottom_data.cluster_degree{2}(ii)/11)
end
title('Second Largest Cluster Degree')
xlim([min(bottom_data.x(myind))-10 max(bottom_data.x(myind))+10])
ylim([min(bottom_data.z(myind))-1 max(bottom_data.z(myind))+1])
subplot(2,3,5); % Markersize proportional to local clustering coefficient
pcolor(bottom_data.xisw,bottom_data.zisw,bottom_data.uisw),shading interp
hold on
for ii = 1:size(myind,2)
    plot(bottom_data.x(myind(ii)),bottom_data.z(myind(ii)),'o','Markersize',10*bottom_data.cluster_coeff{2}(ii)+1,'Color',[1 1 1]*bottom_data.cluster_coeff{2}(ii))
end
title('Second Largest Cluster Coefficient')
xlabel('x (m)')
xlim([min(bottom_data.x(myind))-10 max(bottom_data.x(myind))+10])
ylim([min(bottom_data.z(myind))-1 max(bottom_data.z(myind))+1])
% third largest
subplot(2,3,3); % Markersize proportional to particle degree
pcolor(bottom_data.xisw,bottom_data.zisw,bottom_data.uisw),shading interp
hold on
myind=find(bin==bottom_data.ai(3));
for ii = 1:size(myind,2)
    plot(bottom_data.x(myind(ii)),bottom_data.z(myind(ii)),'o','Markersize',2*bottom_data.cluster_degree{3}(ii)+1,'Color',[1 1 1]*bottom_data.cluster_degree{3}(ii)/11)
end
title('Third Largest Cluster Degree')
xlim([min(bottom_data.x(myind))-10 max(bottom_data.x(myind))+10])
ylim([min(bottom_data.z(myind))-1 max(bottom_data.z(myind))+1])
subplot(2,3,6); % Markersize proportional to local clustering coefficient
pcolor(bottom_data.xisw,bottom_data.zisw,bottom_data.uisw),shading interp
hold on
for ii = 1:size(myind,2)
    plot(bottom_data.x(myind(ii)),bottom_data.z(myind(ii)),'o','Markersize',10*bottom_data.cluster_coeff{3}(ii)+1,'Color',[1 1 1]*bottom_data.cluster_coeff{3}(ii))
end
title('Third Largest Cluster Coefficient')
xlabel('x (m)')
xlim([min(bottom_data.x(myind))-10 max(bottom_data.x(myind))+10])
ylim([min(bottom_data.z(myind))-1 max(bottom_data.z(myind))+1])
drawnow
%print('degree_and_localcoeff_markersize_b','-dtiff','-r650')


figure(15)
clf
betterplots
colormap darkjet
% largest cluster
subplot(2,3,1); % Markersize proportional to particle degree
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw),shading interp
hold on
[bin,~]=conncomp(topneg_data.G);
myind=find(bin==topneg_data.ai(1));
for ii = 1:size(myind,2)
    plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'o','Markersize',2*topneg_data.cluster_degree{1}(ii)+1,'Color',[1 1 1]*topneg_data.cluster_degree{1}(ii)/11)
end
title('Largest Cluster Degree')
ylabel('z (m)')
xlim([min(topneg_data.x(myind))-10 max(topneg_data.x(myind))+10])
ylim([min(topneg_data.z(myind))-1 max(topneg_data.z(myind))+1])
subplot(2,3,4); % Markersize proportional to local clustering coefficient
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw),shading interp
hold on
for ii = 1:size(myind,2)
    plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'o','Markersize',10*topneg_data.cluster_coeff{1}(ii)+1,'Color',[1 1 1]*topneg_data.cluster_coeff{1}(ii))
end
title('Largest Cluster Coefficient')
ylabel('z (m)')
xlabel('x (m)')
xlim([min(topneg_data.x(myind))-10 max(topneg_data.x(myind))+10])
ylim([min(topneg_data.z(myind))-1 max(topneg_data.z(myind))+1])
% second largest
subplot(2,3,2); % Markersize proportional to particle degree
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw),shading interp
hold on
myind=find(bin==topneg_data.ai(2));
for ii = 1:size(myind,2)
    plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'o','Markersize',2*topneg_data.cluster_degree{2}(ii)+1,'Color',[1 1 1]*topneg_data.cluster_degree{2}(ii)/11)
end
title('Second Largest Cluster Degree')
xlim([min(topneg_data.x(myind))-10 max(topneg_data.x(myind))+10])
ylim([min(topneg_data.z(myind))-1 max(topneg_data.z(myind))+1])
subplot(2,3,5); % Markersize proportional to local clustering coefficient
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw),shading interp
hold on
for ii = 1:size(myind,2)
    plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'o','Markersize',10*topneg_data.cluster_coeff{2}(ii)+1,'Color',[1 1 1]*topneg_data.cluster_coeff{2}(ii))
end
title('Second Largest Cluster Coefficient')
xlabel('x (m)')
xlim([min(topneg_data.x(myind))-10 max(topneg_data.x(myind))+10])
ylim([min(topneg_data.z(myind))-1 max(topneg_data.z(myind))+1])
% third largest
subplot(2,3,3); % Markersize proportional to particle degree
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw),shading interp
hold on
myind=find(bin==topneg_data.ai(3));
for ii = 1:size(myind,2)
    plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'o','Markersize',2*topneg_data.cluster_degree{3}(ii)+1,'Color',[1 1 1]*topneg_data.cluster_degree{3}(ii)/11)
end
title('Third Largest Cluster Degree')
xlim([min(topneg_data.x(myind))-10 max(topneg_data.x(myind))+10])
ylim([min(topneg_data.z(myind))-1 max(topneg_data.z(myind))+1])
subplot(2,3,6); % Markersize proportional to local clustering coefficient
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw),shading interp
hold on
for ii = 1:size(myind,2)
    plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'o','Markersize',10*topneg_data.cluster_coeff{3}(ii)+1,'Color',[1 1 1]*topneg_data.cluster_coeff{3}(ii))
end
title('Third Largest Cluster Coefficient')
xlabel('x (m)')
xlim([min(topneg_data.x(myind))-10 max(topneg_data.x(myind))+10])
ylim([min(topneg_data.z(myind))-1 max(topneg_data.z(myind))+1])
drawnow
%print('degree_and_localcoeff_markersize_c','-dtiff','-r650')


%% (ALTERNATE) Degree and Local Coefficient Proportional to MarkerSize in Clusters
figure(16)
clf
betterplots
colormap darkjet
subplot(2,3,1);
pcolor(nou_data.xisw,nou_data.zisw,nou_data.uisw),shading interp
hold on
[bin,~]=conncomp(nou_data.G);
myind=find(bin==nou_data.ai(1));
for ii = 1:size(myind,2)
    plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'o','Markersize',2*nou_data.cluster_degree{1}(ii)+1,'Color',[1 1 1]*nou_data.cluster_degree{1}(ii)/11)
end
lower_x = min(nou_data.x(myind));
upper_x = max(nou_data.x(myind));
lower_y = min(nou_data.z(myind));
upper_y = max(nou_data.z(myind));
myind=find(bin==nou_data.ai(2));
for ii = 1:size(myind,2)
    plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'o','Markersize',2*nou_data.cluster_degree{2}(ii)+1,'Color',[1 1 1]*nou_data.cluster_degree{2}(ii)/11)
end
lower_x = min([nou_data.x(myind); lower_x]);
upper_x = max([nou_data.x(myind); upper_x]);
lower_y = min([nou_data.z(myind); lower_y]);
upper_y = max([nou_data.z(myind); upper_y]);
myind=find(bin==nou_data.ai(3));
for ii = 1:size(myind,2)
    plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'o','Markersize',2*nou_data.cluster_degree{3}(ii)+1,'Color',[1 1 1]*nou_data.cluster_degree{3}(ii)/11)
end
lower_x = min([nou_data.x(myind); lower_x]);
upper_x = max([nou_data.x(myind); upper_x]);
lower_y = min([nou_data.z(myind); lower_y]);
upper_y = max([nou_data.z(myind); upper_y]);
title('Particle Degree in ISW','FontSize',8)
ylabel('z (m)')
xlabel('x (m)')
xlim([lower_x-10 upper_x+10])
ylim([lower_y-1 upper_y+1])
drawnow

subplot(2,3,4);
pcolor(nou_data.xisw,nou_data.zisw,nou_data.uisw),shading interp
hold on
[bin,~]=conncomp(nou_data.G);
myind=find(bin==nou_data.ai(1));
for ii = 1:size(myind,2)
    plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'o','Markersize',10*nou_data.cluster_coeff{1}(ii)+1,'Color',[1 1 1]*nou_data.cluster_coeff{1}(ii))
end
lower_x = min(nou_data.x(myind));
upper_x = max(nou_data.x(myind));
lower_y = min(nou_data.z(myind));
upper_y = max(nou_data.z(myind));
myind=find(bin==nou_data.ai(2));
for ii = 1:size(myind,2)
    plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'o','Markersize',10*nou_data.cluster_coeff{2}(ii)+1,'Color',[1 1 1]*nou_data.cluster_coeff{2}(ii))
end
lower_x = min([nou_data.x(myind); lower_x]);
upper_x = max([nou_data.x(myind); upper_x]);
lower_y = min([nou_data.z(myind); lower_y]);
upper_y = max([nou_data.z(myind); upper_y]);
myind=find(bin==nou_data.ai(3));
for ii = 1:size(myind,2)
    plot(nou_data.x(myind(ii)),nou_data.z(myind(ii)),'o','Markersize',10*nou_data.cluster_coeff{3}(ii)+1,'Color',[1 1 1]*nou_data.cluster_coeff{3}(ii))
end
lower_x = min([nou_data.x(myind); lower_x]);
upper_x = max([nou_data.x(myind); upper_x]);
lower_y = min([nou_data.z(myind); lower_y]);
upper_y = max([nou_data.z(myind); upper_y]);
title('Local Clustering Coefficient in ISW','FontSize',8)
ylabel('z (m)')
xlabel('x (m)')
xlim([lower_x-10 upper_x+10])
ylim([lower_y-1 upper_y+1])
drawnow

subplot(2,3,2);
pcolor(bottom_data.xisw,bottom_data.zisw,bottom_data.uisw),shading interp
hold on
[bin,~]=conncomp(bottom_data.G);
myind=find(bin==bottom_data.ai(1));
for ii = 1:size(myind,2)
    plot(bottom_data.x(myind(ii)),bottom_data.z(myind(ii)),'o','Markersize',2*bottom_data.cluster_degree{1}(ii)+1,'Color',[1 1 1]*bottom_data.cluster_degree{1}(ii)/11)
end
lower_x = min(bottom_data.x(myind));
upper_x = max(bottom_data.x(myind));
lower_y = min(bottom_data.z(myind));
upper_y = max(bottom_data.z(myind));
myind=find(bin==bottom_data.ai(2));
for ii = 1:size(myind,2)
    plot(bottom_data.x(myind(ii)),bottom_data.z(myind(ii)),'o','Markersize',2*bottom_data.cluster_degree{2}(ii)+1,'Color',[1 1 1]*bottom_data.cluster_degree{2}(ii)/11)
end
lower_x = min([bottom_data.x(myind); lower_x]);
upper_x = max([bottom_data.x(myind); upper_x]);
lower_y = min([bottom_data.z(myind); lower_y]);
upper_y = max([bottom_data.z(myind); upper_y]);
myind=find(bin==bottom_data.ai(3));
for ii = 1:size(myind,2)
    plot(bottom_data.x(myind(ii)),bottom_data.z(myind(ii)),'o','Markersize',2*bottom_data.cluster_degree{3}(ii)+1,'Color',[1 1 1]*bottom_data.cluster_degree{3}(ii)/11)
end
lower_x = min([bottom_data.x(myind); lower_x]);
upper_x = max([bottom_data.x(myind); upper_x]);
lower_y = min([bottom_data.z(myind); lower_y]);
upper_y = max([bottom_data.z(myind); upper_y]);
title('Particle Degree in ISW with Bottom Current','FontSize',8)
ylabel('z (m)')
xlabel('x (m)')
xlim([lower_x-10 upper_x+10])
ylim([lower_y-1 upper_y+1])
drawnow

subplot(2,3,5);
pcolor(bottom_data.xisw,bottom_data.zisw,bottom_data.uisw),shading interp
hold on
[bin,~]=conncomp(bottom_data.G);
myind=find(bin==bottom_data.ai(1));
for ii = 1:size(myind,2)
    plot(bottom_data.x(myind(ii)),bottom_data.z(myind(ii)),'o','Markersize',10*bottom_data.cluster_coeff{1}(ii)+1,'Color',[1 1 1]*bottom_data.cluster_coeff{1}(ii))
end
lower_x = min(bottom_data.x(myind));
upper_x = max(bottom_data.x(myind));
lower_y = min(bottom_data.z(myind));
upper_y = max(bottom_data.z(myind));
myind=find(bin==bottom_data.ai(2));
for ii = 1:size(myind,2)
    plot(bottom_data.x(myind(ii)),bottom_data.z(myind(ii)),'o','Markersize',10*bottom_data.cluster_coeff{2}(ii)+1,'Color',[1 1 1]*bottom_data.cluster_coeff{2}(ii))
end
lower_x = min([bottom_data.x(myind); lower_x]);
upper_x = max([bottom_data.x(myind); upper_x]);
lower_y = min([bottom_data.z(myind); lower_y]);
upper_y = max([bottom_data.z(myind); upper_y]);
myind=find(bin==bottom_data.ai(3));
for ii = 1:size(myind,2)
    plot(bottom_data.x(myind(ii)),bottom_data.z(myind(ii)),'o','Markersize',10*bottom_data.cluster_coeff{3}(ii)+1,'Color',[1 1 1]*bottom_data.cluster_coeff{3}(ii))
end
lower_x = min([bottom_data.x(myind); lower_x]);
upper_x = max([bottom_data.x(myind); upper_x]);
lower_y = min([bottom_data.z(myind); lower_y]);
upper_y = max([bottom_data.z(myind); upper_y]);
title('Local Clustering Coefficient in ISW with Bottom Current','FontSize',8)
ylabel('z (m)')
xlabel('x (m)')
xlim([lower_x-10 upper_x+10])
ylim([lower_y-1 upper_y+1])
drawnow

subplot(2,3,3);
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw),shading interp
hold on
[bin,~]=conncomp(topneg_data.G);
myind=find(bin==topneg_data.ai(1));
for ii = 1:size(myind,2)
    plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'o','Markersize',2*topneg_data.cluster_degree{1}(ii)+1,'Color',[1 1 1]*topneg_data.cluster_degree{1}(ii)/11)
end
lower_x = min(topneg_data.x(myind));
upper_x = max(topneg_data.x(myind));
lower_y = min(topneg_data.z(myind));
upper_y = max(topneg_data.z(myind));
myind=find(bin==topneg_data.ai(2));
for ii = 1:size(myind,2)
    plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'o','Markersize',2*topneg_data.cluster_degree{2}(ii)+1,'Color',[1 1 1]*topneg_data.cluster_degree{2}(ii)/11)
end
lower_x = min([topneg_data.x(myind); lower_x]);
upper_x = max([topneg_data.x(myind); upper_x]);
lower_y = min([topneg_data.z(myind); lower_y]);
upper_y = max([topneg_data.z(myind); upper_y]);
myind=find(bin==topneg_data.ai(3));
for ii = 1:size(myind,2)
    plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'o','Markersize',2*topneg_data.cluster_degree{3}(ii)+1,'Color',[1 1 1]*topneg_data.cluster_degree{3}(ii)/11)
end
lower_x = min([topneg_data.x(myind); lower_x]);
upper_x = max([topneg_data.x(myind); upper_x]);
lower_y = min([topneg_data.z(myind); lower_y]);
upper_y = max([topneg_data.z(myind); upper_y]);
title('Particle Degree in ISW with Top-Trapped Negative Current','FontSize',8)
ylabel('z (m)')
xlabel('x (m)')
xlim([lower_x-10 upper_x+10])
ylim([lower_y-1 upper_y+1])
drawnow

subplot(2,3,6);
pcolor(topneg_data.xisw,topneg_data.zisw,topneg_data.uisw),shading interp
hold on
[bin,~]=conncomp(topneg_data.G);
myind=find(bin==topneg_data.ai(1));
for ii = 1:size(myind,2)
    plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'o','Markersize',10*topneg_data.cluster_coeff{1}(ii)+1,'Color',[1 1 1]*topneg_data.cluster_coeff{1}(ii))
end
lower_x = min(topneg_data.x(myind));
upper_x = max(topneg_data.x(myind));
lower_y = min(topneg_data.z(myind));
upper_y = max(topneg_data.z(myind));
myind=find(bin==topneg_data.ai(2));
for ii = 1:size(myind,2)
    plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'o','Markersize',10*topneg_data.cluster_coeff{2}(ii)+1,'Color',[1 1 1]*topneg_data.cluster_coeff{2}(ii))
end
lower_x = min([topneg_data.x(myind); lower_x]);
upper_x = max([topneg_data.x(myind); upper_x]);
lower_y = min([topneg_data.z(myind); lower_y]);
upper_y = max([topneg_data.z(myind); upper_y]);
myind=find(bin==topneg_data.ai(3));
for ii = 1:size(myind,2)
    plot(topneg_data.x(myind(ii)),topneg_data.z(myind(ii)),'o','Markersize',10*topneg_data.cluster_coeff{3}(ii)+1,'Color',[1 1 1]*topneg_data.cluster_coeff{3}(ii))
end
lower_x = min([topneg_data.x(myind); lower_x]);
upper_x = max([topneg_data.x(myind); upper_x]);
lower_y = min([topneg_data.z(myind); lower_y]);
upper_y = max([topneg_data.z(myind); upper_y]);
title('Local Clustering Coefficient in ISW with Top-Trapped Negative Current','FontSize',8)
ylabel('z (m)')
xlabel('x (m)')
xlim([lower_x-10 upper_x+10])
ylim([lower_y-1 upper_y+1])
drawnow

%print('degree_and_localcoeff_markersize_onefig','-dtiff','-r650')