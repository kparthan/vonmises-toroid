function [] = heat_map_3D(file_name)

  % draw a torus 
  r1 = 2;
  r2 = 1;

  [u,v] = meshgrid(0:10:360);

  X = (r1 + r2*cosd(v)).*cosd(u);
  Y = (r1 + r2*cosd(v)).*sind(u);
  Z = r2*sind(v);

  fig = figure();
  hold on;
  axis equal;
  set(gcf, 'Color', 'w');
  surface(X,Y,Z,'FaceColor','none','LineWidth',0.1,'linestyle','--');
  %surfl(X,Y,Z);

  % plot data
  M = load(file_name);
  theta1 = M(:,1);
  theta2 = M(:,2);

  x = (r1 + r2*cos(theta2)).*cos(theta1);
  y = (r1 + r2*cos(theta2)).*sin(theta1);
  z = r2*sin(theta2);

  density = M(:,3);

  h=scatter3(x,y,z,2,'cdata',density);

  xlabel('X');
  ylabel('Y');
  zlabel('Z');

end
