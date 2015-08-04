function [] = kludge()

  addpath('export_fig');

  r1 = 2;
  r2 = 1;

  [u,v] = meshgrid(0:10:360);

  X = (r1 + r2*cosd(v)).*cosd(u);
  Y = (r1 + r2*cosd(v)).*sind(u);
  Z = r2*sind(v);

  fig = figure();
  hold on;
  %box off;
  axis equal;
  set(gca, 'visible', 'off');
  set(gcf, 'Color', 'w');
  %surface(X,Y,Z,'FaceColor','none','LineWidth',0.1,'linestyle','--');
  surface(X,Y,Z,'FaceColor','none','LineWidth',0.01,'linestyle',':');
  %surfl(X,Y,Z);

  M1 = load('../sampled_data/comp1.dat');
  theta1 = M1(:,1);
  theta2 = M1(:,2);

  x = (r1 + r2*cos(theta2)).*cos(theta1);
  y = (r1 + r2*cos(theta2)).*sin(theta1);
  z = r2*sin(theta2);
  plot3(x,y,z,'.','Color',[1 0 0]);

  M2 = load('../sampled_data/comp2.dat');
  theta1 = M2(:,1);
  theta2 = M2(:,2);

  x = (r1 + r2*cos(theta2)).*cos(theta1);
  y = (r1 + r2*cos(theta2)).*sin(theta1);
  z = r2*sin(theta2);
  plot3(x,y,z,'.','Color',[0 0.5 0]);

  view([94 34]);

  xlabel('X');
  ylabel('Y');
  zlabel('Z');

  outfile = 'mixture2_torus';
  output_fig = strcat('../figs/',outfile,'.fig');
  output_eps = strcat('../figs/',outfile,'.eps');
  output_pdf = strcat('../figs/',outfile,'.pdf');

  saveas(gcf,output_fig);
  export_fig(output_pdf,'-pdf');

end
