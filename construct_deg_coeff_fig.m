%% A FIGURE TO SHOW HOW TO CONSTRUCT THE DEGREE/LOCAL CLUSTERING COEFFICIENT PROPORTIONAL TO MARKERSIZE FIGURE %%
figure(1);
clf
betterplots
colormap darkjet


% plot paths of a sample cluster
axes('OuterPosition',[0 0.5 0.5 0.5])
p1 = polyfit([1,19],[-2,-10],2);
x1 = linspace(1,19,100);
y1 = polyval(p1,x1);

p2 = polyfit([1,10,15,18],[-4,-3.5,-5.5,-5],5);
x2 = linspace(1,18,100);
y2 = polyval(p2,x2);

p3 = polyfit([1,6.65,12.5],[-6,-5.8,-6.5],2);
x3 = linspace(1,12.5,100);
y3 = polyval(p3,x3);

p4 = polyfit([1,6,12,21],[-8,-6,-5,-4],3);
x4 = linspace(1,21,100);
y4 = polyval(p4,x4);

p5 = polyfit([1,10,17,20],[-10,-8.5,-5.9,-7],3);
x5 = linspace(1,20,100);
y5 = polyval(p5,x5);

plot(x1,y1,'k-', x2,y2,'-k', x3,y3,'k-', x4,y4,'k-', x5,y5,'k-',...
    [7.7,11.7,14,13.5,16.3,6.65,],[-3.3,-5,-6.4,-4.8,-5.8,-5.8,],'ro','MarkerSize',15)
hold on
plot([x1(1),x1(end),x2(1),x2(end),x3(1),x3(end),x4(1),x4(end),x5(1),x5(end)], ...
    [y1(1),y1(end),y2(1),y2(end),y3(1),y3(end),y4(1),y4(end),y5(1),y5(end)],'k.','MarkerSize',30)
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
text(1,-2.5,'1')
text(1,-4.5,'2')
text(1,-6.5,'3')
text(1,-8.5,'4')
text(1,-10.5,'5')
text(23,-3,'({\it a})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')


% plot interaction graph
axes('OuterPosition',[0 0 0.5 0.5])
adj_matrix = [0, 1, 0, 1, 1;
    1, 0, 0, 1, 1;
    0, 0, 0, 1, 0;
    1, 1, 1, 0, 0;
    1, 1, 0, 0, 0];
G = graph(adj_matrix);
plot(G,'k-')
text(1.5,0.9,'({\it b})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')


% plot the neighbourhood of node 1 and highlight existing edges to help
% calculate degree and local clustering coefficient
axes('OuterPosition',[0.5 4/5 0.25 1/5])
G = graph([0, 1, 0, 1, 1;
    1, 0, 0, 1, 1;
    0, 0, 0, 1, 0;
    1, 1, 1, 0, 1;
    1, 1, 0, 1, 0]);
h = plot(G,'k-');
highlight(h,[1],'Marker','square','MarkerSize',10) % highlight node of interest
highlight(h,[1 1 1],[2 4 5],'EdgeColor','blue') % highlight neighbouring edges
highlight(h,[2 2],[4 5],'EdgeColor','green') % highlight existing edges between neighbours
highlight(h,[4],[5],'LineStyle','--','EdgeColor','red') % highlight missing edges between neighbours
text(1.05,0.85,'({\it c})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')


% plot the neighbourhood of node 2 and highlight existing edges to help
% calculate degree and local clustering coefficient
axes('OuterPosition',[0.5 3/5 0.25 1/5])
G = graph([0, 1, 0, 1, 1;
    1, 0, 0, 1, 1;
    0, 0, 0, 1, 0;
    1, 1, 1, 0, 1;
    1, 1, 0, 1, 0]);
h = plot(G,'k-');
highlight(h,[2],'Marker','square','MarkerSize',10) % highlight node of interest
highlight(h,[2 2 2],[1 4 5],'EdgeColor','blue') % highlight neighbouring edges
highlight(h,[1 1],[4 5],'EdgeColor','green') % highlight existing edges between neighbours
highlight(h,[4],[5],'LineStyle','--','EdgeColor','red') % highlight missing edges between neighbours 
text(1.05,0.85,'({\it d})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')


% plot the neighbourhood of node 3 and highlight existing edges to help
% calculate degree and local clustering coefficient
axes('OuterPosition',[0.5 2/5 0.25 1/5])
G = graph([0, 1, 0, 1, 1;
    1, 0, 0, 1, 1;
    0, 0, 0, 1, 0;
    1, 1, 1, 0, 0;
    1, 1, 0, 0, 0]);
h = plot(G,'k-');
highlight(h,[3],'Marker','square','MarkerSize',10) % highlight node of interest
highlight(h,[3],[4],'EdgeColor','blue') % highlight neighbouring edges 
text(1.4,0.55,'({\it e})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')


% plot the neighbourhood of node 4 and highlight existing edges to help
% calculate degree and local clustering coefficient
axes('OuterPosition',[0.5 1/5 0.25 1/5])
G = graph([0, 1, 1, 1, 1;
    1, 0, 1, 1, 1;
    1, 1, 0, 1, 0;
    1, 1, 1, 0, 0;
    1, 1, 0, 0, 0]);
h = plot(G,'k-');
highlight(h,[4],'Marker','square','MarkerSize',10) % highlight node of interest
highlight(h,[4 4 4],[1 2 3],'EdgeColor','blue') % highlight neighbouring edges
highlight(h,[1],[2],'EdgeColor','green') % highlight existing edges between neighbours
highlight(h,[1 2],[3 3],'LineStyle','--','EdgeColor','red') % highlight missing edges between neighbours 
text(1.15,0.8,'({\it f})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')


% plot the neighbourhood of node 5 and highlight existing edges to help
% calculate degree and local clustering coefficient
axes('OuterPosition',[0.5 0 0.25 1/5])
G = graph([0, 1, 0, 1, 1;
    1, 0, 0, 1, 1;
    0, 0, 0, 1, 0;
    1, 1, 1, 0, 0;
    1, 1, 0, 0, 0]);
h = plot(G,'k-');
highlight(h,[5],'Marker','square','MarkerSize',10) % highlight node of interest
highlight(h,[1 2],[5 5],'EdgeColor','blue') % highlight neighbouring edges
highlight(h,[1],[2],'EdgeColor','green') % highlight existing edges between neighbours
text(1.35,0.55,'({\it g})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')


% calculate degree and local clustering coefficients
deg = sum(adj_matrix);
cluster_size = length(deg);
local_coeff = zeros(cluster_size,1);
for ii = 1:cluster_size
    num_neighbours = deg(ii);

    if num_neighbours == 0 || num_neighbours == 1
        continue
    end

    for jj = 1:cluster_size
        for kk = 1:cluster_size
            local_coeff(ii) = local_coeff(ii) + adj_matrix(ii,jj)*adj_matrix(jj,kk)*adj_matrix(kk,ii);
        end
    end

    local_coeff(ii) = local_coeff(ii)/(num_neighbours*(num_neighbours-1));
end


% plot final positions
axes('OuterPosition',[0.75 0.2 0.25 0.8])
final_x = [x1(end),x2(end),x3(end),x4(end),x5(end)];
final_y = [y1(end),y2(end),y3(end),y4(end),y5(end)];
hold on
for ii = 1:length(final_x)
    plot([final_x(ii)],[final_y(ii)],'o','MarkerSize',2*deg(ii)+5,'Color',[1 1 1]*deg(ii)/(cluster_size-1))
    if local_coeff(ii) > 0
        plot(final_x(ii),final_y(ii),'mx')
    end
end
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
text(19.5,-4.5,'({\it h})','fontname','Times New Roman','FontSize',18,'FontWeight','normal')
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

% legend underneath plot
degree_max = cluster_size-1;
axes('OuterPosition',[0.75 0.025 0.25 0.2],'NextPlot','add','DataAspectRatio',[1 1 1],'XLim',[0 degree_max],'YLim',[0 0.0001])
sample_points = linspace(0,degree_max,degree_max+1);
hold on
for ii = 1:length(sample_points)
    plot(sample_points(ii), [0],'o','Markersize',2*sample_points(ii)+5,'Color',[1 1 1]*sample_points(ii)/(cluster_size-1))
end
xlabel('$d_G(v_i)$','Interpreter','latex')
xticks(1:degree_max)
set(gca,'fontname','Times New Roman','FontSize',18,'FontWeight','normal')

% save figure
%print('Figure9','-depsc','-r650')