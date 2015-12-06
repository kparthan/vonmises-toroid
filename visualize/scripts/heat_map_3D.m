% file_name = "mixture_density.dat" (N X 3) [t1, t2, density] 

function [] = heat_map_3D(file_name)

  addpath('export_fig');

  fig = figure();
  hold on;
  axis equal;
  set(gcf, 'Color', 'w');
  set(gcf,'defaultaxesfontname','Arial');

  % draw a torus 
  r1 = 2;
  r2 = 1;
  [u,v] = meshgrid(0:10:360);
  X = (r1 + r2*cosd(v)).*cosd(u);
  Y = (r1 + r2*cosd(v)).*sind(u);
  Z = r2*sind(v);
  surface(X,Y,Z,'FaceColor','none','LineWidth',0.1,'linestyle',':');

  % plot data
  M = load(file_name);
  theta1 = M(:,1);
  theta2 = M(:,2);

  x = (r1 + r2*cos(theta2)).*cos(theta1);
  y = (r1 + r2*cos(theta2)).*sin(theta1);
  z = r2*sin(theta2);

  density = M(:,3);
  max_density = max(density);
  caxis([0 max_density]); % colormap range

  h=scatter3(x,y,z,2,'cdata',density);

  % plot children means %

%  p1 = [45, 60]; 
%  p2 = [30, 45]; 
%  p3 = [50, 45];
%  plot_point_on_torus(p1,'k');
%  plot_point_on_torus(p2,'k');
%  plot_point_on_torus(p3,'k');
%  c11 = [47.977, 62.977]; c12 = [42.023, 57.023]; 
%  c21 = [33.839, 41.161]; c22 = [26.161, 48.839]; 
%  c31 = [57.759, 52.759]; c32 = [42.241, 37.241];
%  plot_point_on_torus(c11,'r'); plot_point_on_torus(c12,'r');
%  plot_point_on_torus(c21,'r'); plot_point_on_torus(c22,'r');
%  plot_point_on_torus(c31,'r'); plot_point_on_torus(c32,'r');

  % plot children means %

  %set(gca, 'visible', 'off');
  %view ([-139 40]);

  axis off;
  %xlabel('X','fontsize',20);
  %ylabel('Y','fontsize',20);
  %zlabel('Z');

  file_name = '../figs/torus_lambda_5';
  output_fig = strcat(file_name,'.fig');
  output_eps = strcat(file_name,'.eps');
  output_pdf = strcat(file_name,'.pdf');

  %saveas(gcf,output_fig);

  %print2eps(output_eps);
  %eps2pdf(output_eps,output_pdf,1);

  %export_fig(output_pdf,'-pdf');

end

% point in degrees
function [] = plot_point_on_torus(point, color)

  r1 = 2;
  r2 = 1;

  theta1 = point(1); theta2 = point(2);

  x = (r1 + r2*cosd(theta2))*cosd(theta1);
  y = (r1 + r2*cosd(theta2))*sind(theta1);
  z = r2*sind(theta2);

  format = strcat(color,'.');
  plot3(x,y,z,format,'markersize',10);

end

