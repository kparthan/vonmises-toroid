function [] = plot_torus()

  addpath('export_fig');

  r1 = 2;
  r2 = 1;

  % plot the outline of the torus
  [u,v] = meshgrid(0:10:360);
  X = (r1 + r2*cosd(v)).*cosd(u);
  Y = (r1 + r2*cosd(v)).*sind(u);
  Z = r2*sind(v);

  fig = figure();
  hold on;
  set(gcf, 'Color', 'w');
  axis equal;
  axis off;

  surface(X,Y,Z,'FaceColor','none','LineWidth',0.001,'linestyle','-');
  %surfl(X,Y,Z);

  x1 = [0 3]; y1 = [0 0];
  x2 = [0 3*cosd(50)]; y2 = [0 3*sind(50)];

%  drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0 );
%  drawArrow(x1,y1);
%  drawArrow = @(x,y,props) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,props{:});
%  drawArrow(x1,y1,{'MaxHeadSize',0.8,'Color','b','LineWidth',3});

  line(x1,y1,'LineWidth',2,'Color','b','linestyle',':');
  line(x2,y2,'LineWidth',2,'Color','b','linestyle',':');

  phi = '\phi';
  text(0.5,0.3,phi,'fontsize',15,'Color','b');
  plot_partial_circle(0.4);

  xlabel('X');
  ylabel('Y');
  zlabel('Z');

  outfile = 'torus_phi';
  output_fig = strcat('../figs/',outfile,'.fig');
  output_eps = strcat('../figs/',outfile,'.eps');
  output_pdf = strcat('../figs/',outfile,'.pdf');

  saveas(gcf,output_fig);
  export_fig(output_pdf,'-pdf');
  print2eps(output_eps);

end

function [] = plot_partial_circle(r)

  for theta=0:5:50 
    x = r * cosd(theta);
    y = r * sind(theta);
    plot(x,y,'.','Color','b','LineWidth',2,'linestyle',':');
  end

end

