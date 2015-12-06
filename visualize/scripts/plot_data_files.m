% data_file: N X 2 angle pairs (in radians)
function [] = plot_data_files()

  fig = figure();
  hold on;
  axis equal;
  set(gcf, 'Color', 'w');

  % plot the outline of the torus
  [u,v] = meshgrid(0:10:360);
  r1 = 2;
  r2 = 1;

  X = (r1 + r2*cosd(v)).*cosd(u);
  Y = (r1 + r2*cosd(v)).*sind(u);
  Z = r2*sind(v);

  %surface(X,Y,Z,'FaceColor','none','LineWidth',0.1,'linestyle','--');
  surface(X,Y,Z,'FaceColor','none','LineWidth',0.01,'linestyle',':');

  plot_file('../sampled_data/bins_sine/comp1.dat',[1 0 0]);
  plot_file('../sampled_data/bins_sine/comp2.dat',[0 1 0]);
  plot_file('../sampled_data/bins_sine/comp3.dat',[0 0 1]);

  view([-143 38]);

  xlabel('X');
  ylabel('Y');
  zlabel('Z');

end

function [] = plot_file(data_file,color_vec)

  r1 = 2;
  r2 = 1;

  M = load(data_file);
  theta1 = M(:,1);
  theta2 = M(:,2);

  x = (r1 + r2*cos(theta2)).*cos(theta1);
  y = (r1 + r2*cos(theta2)).*sin(theta1);
  z = r2*sin(theta2);
  plot3(x,y,z,'.','Color',color_vec);

end
