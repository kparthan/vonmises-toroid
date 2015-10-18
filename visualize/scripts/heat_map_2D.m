% file_name = "prob_bins.dat" (360 X 360) && bins = 1
% file_name = "mixture_density.dat" (N X 3) [t1, t2, density] && bins = 0

function [] = heat_map_2D(file_name, bins)

  %addpath('export_fig');

  %[theta phi] = meshgrid(0:1:359.9,0:1:359.9);
  fig = figure();
  hold on;
  axis equal;
  set(gcf, 'Color', 'w');
  xlabel('\phi','fontsize',12);
  ylabel('\theta','fontsize',12);
  set(gca,'Xlim',[0 360]);
  set(gca,'Ylim',[0 360]);
  set(gca,'XTick',[0:45:360]);
  set(gca,'YTick',[0:45:360]);
  %view ([0 90]);

  M = load(file_name);

  if bins == 1
    mesh(M);
  else
    theta = M(:,1) .* 180/pi;
    phi = M(:,2) .* 180/pi;
    density = M(:,3);
    scatter(phi,theta,1,'cdata',density);
  end
  
end
