function [] = heat_map_2D(file_name)

  %addpath('export_fig');

  %[theta phi] = meshgrid(0:1:359.9,0:1:359.9);
  theta = 0:1:359.9;
  phi = 0:1:359.9;

  fig = figure();
  hold on;
  set(gcf, 'Color', 'w');
  xlabel('X');
  ylabel('Y');
  %view ([0 90]);

  set(gca,'Xlim',[0 360]);
  set(gca,'Ylim',[0 360]);

  M = load(file_name);
  mesh(M);
  %surf(M);
  
end
