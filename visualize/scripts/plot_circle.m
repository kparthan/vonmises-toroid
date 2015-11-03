function [] = plot_circle(r)

  addpath('export_fig');

  fig = figure();
  hold on;
  set(gcf, 'Color', 'w');
  axis equal;
  axis off;

  for theta=0:1:360 
    x = r * cosd(theta);
    y = r * sind(theta);
    plot(x,y,'.','Color','k','LineWidth',2,'linestyle',':');
  end

  x1 = [0 r]; y1 = [0 0];
  x2 = [0 r*cosd(50)]; y2 = [0 r*sind(50)];

  line(x1,y1,'LineWidth',2,'Color','b','linestyle',':');
  line(x2,y2,'LineWidth',2,'Color','b','linestyle',':');

  psi = '\psi';
  text(0.2,0.1,psi,'fontsize',15,'Color','b');
  plot_partial_circle(0.15);

  outfile = 'torus_psi';
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
