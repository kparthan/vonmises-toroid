function [] = visualize_mixture_contours_cdf(K)

  addpath('export_fig');

  bins_folder = '../sampled_data/';
 
  % figure properties
  fig = figure();
  caxis([0 1]); % colormap range
  axis equal
  hold on;
  set(gcf, 'Color', 'w');
  xlabel('\theta_1','fontsize',15);
  ylabel('\theta_2','fontsize',15);
  %xlabel('\theta','fontsize',12);
  %ylabel('\phi','fontsize',12);
  xlabh = get(gca,'XLabel');
  ylabh = get(gca,'YLabel');
  set(xlabh,'interpreter','tex');
  set(ylabh,'interpreter','tex');
  set(gca,'Xlim',[0 360]);
  set(gca,'Ylim',[0 360]);
  set(gca,'xtick',[0:60:360],'fontsize',10);
  set(gca,'ytick',[0:60:360],'fontsize',10);
  %set(gca,'xtick',[0:30:90],'fontsize',15);
  %set(gca,'ytick',[0:30:90],'fontsize',15);
  %view ([0 0]);

  % plot the contours 
  phi = 0:1:359.9;  % meshgrid columns
  theta = 0:1:359.9;  % meshgrid rows 
%  min_val = 0;
%  max_val = 0;
  for k = 1:K
    %k
    data_file = strcat(bins_folder,'comp',num2str(k),'_prob_bins2D.dat');

    Mt = load(data_file);
    prob_bins = Mt';
    [max_val_bins max_index] = max(prob_bins(:));
    prob_bins = prob_bins / sum(sum(prob_bins));
    [sorted_prob_bins, indexes] = sort(prob_bins(:),1,'descend');
    sorted_prob_bins_cumsum = cumsum(sorted_prob_bins);
    cdf_bins_list = ones(length(sorted_prob_bins),1);
    cdf_bins_list(indexes) = sorted_prob_bins_cumsum;
    cdf_bins = reshape(cdf_bins_list,size(prob_bins));

    min_val = min(cdf_bins(:));
    max_val = max(cdf_bins(:));
    range = max_val - min_val;
    cdf_bins = (cdf_bins - min_val) / range; % in [0,1]

    level = 0.8;
    norm_level = (level - min_val) / range;
    contour_levels = [norm_level norm_level];
    [C,h] = contour(cdf_bins,contour_levels,'LineWidth',1.5,'LineColor','black');
    %[C,h] = contour(cdf_bins,1,'LineWidth',2,'LineColor','black');
    %clabel(C,h);

%    [c1,h1] = contour(cdf_bins,contour_levels,'LineWidth',2,'LineColor','black');
%    % Take all the info from the contourline output argument:
%    i0 = 1;
%    i2 = 1;       
%    while i0 <  length(c1)
%        i1 = i0+[1:c1(2,i0)];
%        zLevel(i2) = c1(1,i0);
%        hold on
%        % And plot it with dashed lines:
%        %ph(i2) = plot(c1(1,i1),c1(2,i1),'k--','linewidth',2); 
%        ph(i2) = dashline(c1(1,i1),c1(2,i1),1.5,3,1.5,3,'Color','black','linewidth',2); 
%        i0 = i1(end)+1;
%        i2 = i2+1;
%    end
%    % Scrap the contourlines:
%    delete(h1);

    [row col] = ind2sub(size(prob_bins),max_index);
    cx = phi(col);
    cy = theta(row);
    %ht = text(cx,cy,num2str(k),'Color','red');

%    hcl = clabel(C,'Color','red');
%    for i=2:2:length(hcl)
%      old_label = get(hcl(i),'String');
%      new_label = num2str(k);
%      set(hcl(i),'String',new_label);
%    end
  end  

%  mus = [45.11755, 50.52585];
%  child1 = [58.41809, 45.31592];
%  child2 = [37.69036, 53.48310];
%  plot_init(mus,child1,child2);

  % plot the entire mixture density
  %data_file = strcat(bins_folder,'mixture_density.dat');
  data_file = strcat(bins_folder,'mixture_density.dat');
  M = load(data_file);

  density = M(:,3);
  min1 = min(density);
  max1 = max(density);
  range1 = max1 - min1;
  norm_density = (density - min1) / range1; % in [0,1]

  n = size(M,1);
  angles = zeros(n,2);  % in radians
  angles(:,1) = M(:,1) .* 180/pi;
  angles(:,2) = M(:,2) .* 180/pi;

  hs = scatter3(angles(:,1),angles(:,2),norm_density,0.1,'cdata',norm_density);

  %colorbar
  %outfile = 'original_mix';
  %outfile = 'iter1_parent';
  %outfile = 'iter1_init_c1';
  outfile = 'class_b';
  output_fig = strcat('../figs/',outfile,'.fig');
  output_eps = strcat('../figs/',outfile,'.eps');
  output_pdf = strcat('../figs/',outfile,'.pdf');

  %saveas(gcf,output_fig);
  export_fig(output_pdf,'-pdf');
  %print2eps(output_eps);
  %eps2pdf(output_eps,output_pdf);

end

function [] = plot_init(mus,child1,child2)

  plot(mus(1),mus(2),'r.','Markersize',35);
  plot(child1(1),child1(2),'k.','Markersize',35);
  plot(child2(1),child2(2),'k.','Markersize',35);

%  plot(mus(1),mus(2),'k.');
%  plot(child1(1),child1(2),'r.');
%  plot(child2(1),child2(2),'r.');
end

