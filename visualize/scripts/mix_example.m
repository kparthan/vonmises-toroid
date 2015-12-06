function [] = mix_example(K)

  addpath('export_fig');

  bins_folder = '../sampled_data/bins_sine/';
 
  % figure properties
  fig = figure();
  hold on;
  axis equal;
  set(gcf,'defaultaxesfontname','Arial');
  caxis([0 1]); % colormap range
  xlabh = get(gca,'XLabel');
  ylabh = get(gca,'YLabel');
  xlabel('\phi','fontsize',20);
  ylabel('\theta','fontsize',20);
  %xlabel('Longitude','fontsize',20);
  %ylabel('Co-latitude','fontsize',20);

  set(gcf, 'Color', 'w');
  set(gca,'Xlim',[0 90],'Ylim',[0 90]);
  set(gca,'xtick',[0:30:90],'fontsize',15);
  set(gca,'ytick',[0:30:90],'fontsize',15);

  % plot the contours 
  phi = 0:1:359.9;  % meshgrid columns
  theta = 0:1:179.9;  % meshgrid rows 
  for k = 1:K
    %k
    data_file = strcat(bins_folder,'comp',num2str(k),'_prob_bins2D.dat');

    prob_bins = load(data_file);
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

    level = 0.9;
    norm_level = (level - min_val) / range;
    contour_levels = [norm_level norm_level];

    % parent/ split post-EM/ kill post-EM/ merge post-EM
    [C,h] = contour(cdf_bins,contour_levels,'LineWidth',2,'LineColor','black');

    % split initialization
%    split = 3;
%    if (k == split) 
%      colour = 'red';
%    else
%      colour = 'black';
%    end
%    [C,h] = contour(cdf_bins,contour_levels,'LineWidth',2,'LineColor',colour);
%    m1 = [0.774, 0.512, 0.373];
%    m2 = [0.602, 0.737, 0.308];
%    [p1 t1] = cartesianVector2spherical(m1);  % p1 = phi1 (plotted on x-axis)
%    [p2 t2] = cartesianVector2spherical(m2);  % t1 = theta1 (plotted on y-axis)
%    plot(p1,t1,'o','Color',[0.0 0.0 0.0],'MarkerSize',10,'MarkerFaceColor',[0.0,0.0,0.0]);
%    plot(p2,t2,'o','Color',[0.0 0.0 0.0],'MarkerSize',10,'MarkerFaceColor',[0.0,0.0,0.0]);
%    mp = [0.695, 0.631, 0.344]; 
%    [pm tm] = cartesianVector2spherical(mp);
%    plot(pm,tm,'o','Color',[1.0 0.0 0.0],'MarkerSize',10,'MarkerFaceColor',[1.0,0.0,0.0]);

    % split optimized children/before adjustment
%    if (k == split || k == (split+1)) 
%      colour = 'red';
%    else
%      colour = 'black';
%    end
%    plot_dashed_contours(cdf_bins,contour_levels,colour);

    % kill operation
%    kill = 1;
%    if (k == kill) 
%      colour = [0 0.5 0];
%    else
%      colour = 'black';
%    end
%    [C,h] = contour(cdf_bins,contour_levels,'LineWidth',2,'LineColor',colour);

    % kill pre-EM / before adjustment 
%    if (k ~= kill)
%      plot_dashed_contours(cdf_bins,contour_levels,'black');
%    end

    % merge operation
%    merge1 = 3; merge2 = 2;
%    if (k == merge1 || k == merge2)
%      colour = 'blue';
%    else
%      colour = 'black';
%    end
%    [C,h] = contour(cdf_bins,contour_levels,'LineWidth',2,'LineColor',colour);

    % merge pre-EM / before adjustment
%    num_parent_components = 3;  % number of parent components
%    merged_position = num_parent_components - 1;
%    if (k ~= merged_position)
%      colour = 'black';
%    else
%      colour = 'blue';
%    end
%    plot_dashed_contours(cdf_bins,contour_levels,colour);

  end  

  % plot the entire mixture density
  data_file = strcat(bins_folder,'mixture_density.dat');
  M = load(data_file);

  density = M(:,4);
  min1 = min(density);
  max1 = max(density);
  range1 = max1 - min1;
  norm_density = (density - min1) / range1; % in [0,1]
  n = size(M,1);
  angles = zeros(n,2);
  for i=1:n
    x = M(i,1);
    y = M(i,2);
    z = M(i,3);
    [phi, theta] = cartesian2spherical(x,y,z);
    angles(i,1) = phi;
    angles(i,2) = theta;
  end
  hs = scatter3(angles(:,1),angles(:,2),norm_density,2,'cdata',norm_density);

  %outfile = 'original_mix';
  %outfile = 'iter1_parent';
  %outfile = 'iter1_init_c1';
  %outfile = 'iter2_parent';
  %outfile = 'iter2_init_c1';
  %outfile = 'iter2_split_c1_before_adjustment';
  %outfile = 'iter2_split_c1_post_EM';
  %outfile = 'iter2_init_c2';
  %outfile = 'iter2_split_c2_before_adjustment';
  %outfile = 'iter2_split_c2_post_EM';
  %outfile = 'iter2_delete_c1';
  %outfile = 'iter2_delete_c1_before_adjustment';
  %outfile = 'iter2_delete_c1_post_EM';
  %outfile = 'iter2_delete_c2';
  %outfile = 'iter2_delete_c2_before_adjustment';
  %outfile = 'iter2_delete_c2_post_EM';
  %outfile = 'iter2_merge_c1_c2';
  %outfile = 'iter2_merge_c1_c2_before_adjustment';
  %outfile = 'iter2_merge_c1_c2_post_EM';
  %outfile = 'iter3_init_c1';
  %outfile = 'iter3_split_c1_before_adjustment';
  %outfile = 'iter3_split_c1_post_EM';
  %outfile = 'iter3_delete_c1';
  %outfile = 'iter3_delete_c1_before_adjustment';
  %outfile = 'iter3_delete_c1_post_EM';
  %outfile = 'iter3_merge_c1_c2';
  %outfile = 'iter3_merge_c1_c2_before_adjustment';
  %outfile = 'iter3_merge_c1_c2_post_EM';
  %outfile = 'iter3_init_c2';
  %outfile = 'iter3_split_c2_before_adjustment';
  %outfile = 'iter3_split_c2_post_EM';
  %outfile = 'iter3_delete_c2';
  %outfile = 'iter3_delete_c2_before_adjustment';
  %outfile = 'iter3_delete_c2_post_EM';
  %outfile = 'iter3_merge_c2_c3';
  %outfile = 'iter3_merge_c2_c3_before_adjustment';
  %outfile = 'iter3_merge_c2_c3_post_EM';
  %outfile = 'iter3_init_c3';
  %outfile = 'iter3_split_c3_before_adjustment';
  %outfile = 'iter3_split_c3_post_EM';
  %outfile = 'iter3_delete_c3';
  %outfile = 'iter3_delete_c3_before_adjustment';
  %outfile = 'iter3_delete_c3_post_EM';
  %outfile = 'iter3_merge_c3_c2';
  %outfile = 'iter3_merge_c3_c2_before_adjustment';
  %outfile = 'iter3_merge_c3_c2_post_EM';
  output_fig = strcat('../figs/mix_example/',outfile,'.fig');
  output_pdf = strcat('../figs/mix_example/',outfile,'.pdf');

  %saveas(gcf,output_fig);
  %export_fig(output_pdf,'-pdf');

end

function [] = plot_dashed_contours(cdf_bins, contour_levels, colour)

  [C,h] = contour(cdf_bins,contour_levels,'LineWidth',2,'LineColor',colour,'LineStyle','--');
  % Take all the info from the contourline output argument:
  i0 = 1; i2 = 1;
  while i0 <  length(C)
      i1 = i0+[1:C(2,i0)];
      zLevel(i2) = C(1,i0);
      hold on
      % And plot it with dashed lines:
      ph(i2) = dashline(C(1,i1),C(2,i1),1.5,3,1.5,3,'Color',colour,'linewidth',2); 
      i0 = i1(end)+1;
      i2 = i2+1;
  end
  % Scrap the contourlines:
  delete(h);  

end

function [phi, theta] = cartesianVector2spherical(v)

  x = v(1);
  y = v(2);
  z = v(3);

  theta = acos(x);
  ratio = y/sin(theta);
  if (ratio > 1)
    ratio = 1;
  elseif ratio < -1
    ratio = -1;
  end
  angle = acos(ratio);
  phi = 0;
  if (y == 0 && z == 0)
    phi = 0;
  elseif (y == 0) 
    if (z > 0) 
      phi = angle;
    else
      phi = 2 * pi - angle;
    end
  elseif (z >= 0)
    phi = angle;
  elseif (z < 0)
    phi = 2 * pi - angle;
  end
  phi = phi * 180/pi;
  theta = theta * 180/pi;

end
