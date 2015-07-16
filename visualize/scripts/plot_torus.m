function [] = plot_torus(data_file)

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
  %surface(X,Y,Z,'FaceColor','none','LineWidth',0.1,'linestyle','--');
  surfl(X,Y,Z);

  M = load(data_file);
  x = M(:,1);
  y = M(:,2);
  z = M(:,3);
  plot3(x,y,z,'.','Color',[1 0 0]);

%  M = load('../../bvm_cosine2.dat');
%  x = M(:,1);
%  y = M(:,2);
%  z = M(:,3);
%  plot3(x,y,z,'.','Color',[0 0 1]);

  xlabel('X');
  ylabel('Y');
  zlabel('Z');

end
