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
  %surface(X,Y,Z,'FaceColor','none','LineWidth',0.1,'linestyle',':');

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

  h=scatter3(x,y,z,1,'cdata',density);

  %set(gca, 'visible', 'off');
  %view ([-139 40]);
  %view ([61 -90]);

  axis off;
  %xlabel('X','fontsize',20);
  %ylabel('Y','fontsize',20);
  %zlabel('Z');

  file_name = '../figs/torus_lambda_5';
  output_fig = strcat(file_name,'.fig');
  output_eps = strcat(file_name,'.eps');
  output_pdf = strcat(file_name,'.pdf');

  saveas(gcf,output_fig);

  %print2eps(output_eps);
  %eps2pdf(output_eps,output_pdf,1);

  export_fig(output_pdf,'-pdf');

end
