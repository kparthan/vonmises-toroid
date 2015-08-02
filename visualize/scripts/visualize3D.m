function [] = visualize3D(K)

  addpath('export_fig');

  % draw a unit sphere
  n = 25;
  r = ones(n, n); % radius is 1 
  [thetas, phis] = meshgrid(linspace(0, pi, n), linspace(0, 2*pi, n));

  X = zeros(n,n);
  Y = zeros(n,n);
  Z = zeros(n,n);
  for i=1:n
    for j=1:n
      theta = thetas(i,j);
      phi = phis(i,j);
      [x,y,z] = spherical2cartesian(theta,phi);
      X(i,j) = x;
      Y(i,j) = y;
      Z(i,j) = z;
    end
  end

  fig = figure();
  hold on;
  axis equal;
  set(gcf, 'Color', 'w');
  %box on;
  %set(gca,'xgrid','on');
  surface(X,Y,Z,'FaceColor','none','LineWidth',0.01,'linestyle',':');
  view([-55 -30]);
  xlabh = get(gca,'XLabel');
  ylabh = get(gca,'YLabel');
  zlabh = get(gca,'ZLabel');
  set(xlabh,'Position',[0 1.3 -1]);
  set(ylabh,'Position',[1.7 0 -1]);
  %set(zlabh,'Position',[-1 1.2 0]);

  c(1,:) = [1 0 0];
  c(2,:) = [0 0.5 0];
  c(3,:) = [0 0 1];
  % plot the sampled data
  for k = 1:K
     data_file = strcat('../sampled_data/comp',num2str(k),'.dat');
     M = load(data_file);
     x = M(:,1);
     y = M(:,2);
     z = M(:,3);
     plot3(x,y,z,'.','Color',c(k,:));
     %colors = rand(1,3);
     %plot3(x,y,z,'.','Color',colors);
  end  
  xlabel('X_{1}','fontsize',18);
  ylabel('X_{2}','fontsize',18);
  zlabel('X_{3}','fontsize',18);
  set(gca,'fontsize',12);

  %draw_box_lines();

  % create legend
  %N = [1:K];
  %legend_cell = cellstr('unit sphere');
  %legend_cell = [legend_cell ; cellstr(num2str(N','%d'))];
  %legend(legend_cell);

  outfile = 'mixture_sphere';
  output_fig = strcat('../figs/',outfile,'.fig');
  output_eps = strcat('../figs/',outfile,'.eps');
  output_pdf = strcat('../figs/',outfile,'.pdf');

  %saveas(gcf,output_fig);
  %export_fig(output_pdf,'-pdf');

end

function [] = draw_box_lines()

  px = [-1 -1];
  py = [1 -1];
  pz = [-1 -1];
  plot3(px,py,pz,'k');

  px = [1 -1];
  py = [-1 -1];
  pz = [-1 -1];
  plot3(px,py,pz,'k');

  px = [-1 -1];
  py = [-1 -1];
  pz = [1 -1];
  plot3(px,py,pz,'k');

  px = [-1 -1];
  py = [1 -1];
  pz = [1 1];
  plot3(px,py,pz,'k');

  px = [1 -1];
  py = [-1 -1];
  pz = [1 1];
  plot3(px,py,pz,'k');

  px = [1 1];
  py = [-1 -1];
  pz = [-1 1];
  plot3(px,py,pz,'k');

end

