% data_file: N X 2 angle pairs (in radians)
function [] = plot_torus_data(data_file)

  r1 = 2;
  r2 = 1;

  % plot the outline of the torus
  [u,v] = meshgrid(0:10:360);
  X = (r1 + r2*cosd(v)).*cosd(u);
  Y = (r1 + r2*cosd(v)).*sind(u);
  Z = r2*sind(v);

  fig = figure();
  hold on;
  axis equal;
  set(gcf, 'Color', 'w');
  %surface(X,Y,Z,'FaceColor','none','LineWidth',0.1,'linestyle','--');
  surface(X,Y,Z,'FaceColor','none','LineWidth',0.01,'linestyle',':');
  %surfl(X,Y,Z);

  M = load(data_file);
  theta1 = M(:,1);
  theta2 = M(:,2);

  x = (r1 + r2*cos(theta2)).*cos(theta1);
  y = (r1 + r2*cos(theta2)).*sin(theta1);
  z = r2*sin(theta2);
  plot3(x,y,z,'.','Color',[1 0 0]);

  view([-143 38]);

  xlabel('X');
  ylabel('Y');
  zlabel('Z');

end
