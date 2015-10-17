function [] = visualize_mixture_contours_cdf(K,pdf)

  addpath('export_fig');

  isvmf = strcmp(pdf,'vmf');
  bins_folder = '';
  outfile = '';
  if (isvmf == 1)
    bins_folder = '../sampled_data/bins_vmf/';
    %outfile = 'vmf_mix_reweighted_beta';
    %outfile = 'original';
    %outfile = 'b_vmf_37';
    outfile = 'b_vmf';
  elseif (isvmf == 0)
    bins_folder = '../sampled_data/bins_kent/';
    %outfile = 'kent_mix';
    %outfile = 'b_kent_mml_23';
    outfile = 'b_kent_mml_arun';
  end
  disp(strcat('bins_folder: ',{' '},bins_folder)); 
 
  % figure properties
  fig = figure();
  caxis([0 1]); % colormap range
  %axis equal
  hold on;
  set(gcf, 'Color', 'w');
  xlabel('Longitude\phi','fontsize',12);
  ylabel('Co-latitude\theta','fontsize',13);
  xlabh = get(gca,'XLabel');
  ylabh = get(gca,'YLabel');
  set(xlabh,'interpreter','tex');
  set(ylabh,'interpreter','tex');
  set(gca,'Ylim',[15 110]);
  %set(gca,'Ylim',[25 100]);
  set(gca,'Xlim',[0 360]);
  set(gca,'xtick',[0:60:360],'fontsize',10);
  set(gca,'ytick',[15:20:180],'fontsize',10);
  %set(gca,'ytick',[0:25:100],'fontsize',12);
  view ([0 90]);

  % plot the contours 
  phi = 0:1:359.9;  % meshgrid columns
  theta = 0:1:179.9;  % meshgrid rows 
%  min_val = 0;
%  max_val = 0;
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
 %   ht = text(cx,cy,num2str(k),'Color','red');
 %   [cx,cy,index] = number_component(isvmf,k);
 %   if (index > 0)
 %     ht = text(cx,cy,num2str(index),'Color','red','fontsize',8);
 %   end

%    hcl = clabel(C,'Color','red');
%    for i=2:2:length(hcl)
%      old_label = get(hcl(i),'String');
%      new_label = num2str(k);
%      set(hcl(i),'String',new_label);
%    end
  end  

  %uistack(hs,'bottom');
  %uistack(ht,'top');

  % plot the entire mixture density
  data_file = strcat(bins_folder,'mixture_density.dat');
  M = load(data_file);

  density = M(:,4);
  min1 = min(density);
  max1 = max(density);
  range1 = max1 - min1;
  %range = max_val - min_val;
  %norm_density = min_val + (((density - min1) / range1) * range);
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
%  min(norm_density)
%  max(norm_density)
  hs = scatter3(angles(:,1),angles(:,2),norm_density,0.1,'cdata',norm_density);

  %colorbar
  output_fig = strcat('../figs/protein_modelling/',outfile,'.fig');
  output_eps = strcat('../figs/protein_modelling/',outfile,'.eps');
  output_pdf = strcat('../figs/protein_modelling/',outfile,'.pdf');

  %saveas(gcf,output_fig);
  %print2eps(output_eps);
  %eps2pdf(output_eps,output_pdf,1);
  %export_fig(output_pdf,'-pdf');

end


function [cx,cy,index] = number_component(isvmf,k)

  cx = 1; cy = 1; index = -1;
  if (isvmf == 1) % vMF mixture
    if (k == 1)
      cx = 25; cy = 97; index = 1;
    elseif (k == 32)
      cx = 40; cy = 96; index = 2;
    elseif (k == 30)
      cx = 48; cy = 96; index = 3;
    elseif (k == 31)
      cx = 55; cy = 98; index = 4;
    elseif (k == 8)
      cx = 65; cy = 73; index = 5;
    elseif (k == 5)
      cx = 63; cy = 98; index = 6;
    elseif (k == 6)
      cx = 70; cy = 99; index = 7;
    elseif (k == 7)
      cx = 90; cy = 102; index = 8;
    elseif (k == 9)
      cx = 110; cy = 101; index = 9;
    elseif (k == 2)
      cx = 7; cy = 62; index = 10;
    elseif (k == 3)
      cx = 32; cy = 65; index = 11;
    elseif (k == 17)
      cx = 170; cy = 30; index = 12;
    elseif (k == 16)
      cx = 190; cy = 35; index = 13;
    elseif (k == 23)
      cx = 200; cy = 40; index = 14;
    elseif (k == 24)
      cx = 235; cy = 45; index = 15;
    elseif (k == 22)
      cx = 265; cy = 55; index = 16;
    elseif (k == 25)
      cx = 270; cy = 40; index = 17;
    elseif (k == 12)
      cx = 168; cy = 95; index = 18;
    elseif (k == 13)
      cx = 181; cy = 97; index = 19;
    elseif (k == 14)
      cx = 196; cy = 97; index = 20;
    elseif (k == 18)
      cx = 208; cy = 99; index = 21;
    elseif (k == 19)
      cx = 221; cy = 100; index = 22;
    elseif (k == 20)
      cx = 235; cy = 100; index = 23;
    elseif (k == 21)
      cx = 250; cy = 100; index = 24;
    end
  elseif (isvmf == 0) % Kent mixture
    if (k == 17)
      cx = 10; cy = 92; index = 1;
    elseif (k == 16)
      cx = 30; cy = 95; index = 2;
    elseif (k == 11)
      cx = 38; cy = 92; index = 3;
    elseif (k == 10)
      cx = 67; cy = 98; index = 4;
    elseif (k == 15)
      cx = 90; cy = 98; index = 5;
    elseif (k == 14)
      cx = 115; cy = 98; index = 6;
    elseif (k == 23)
      cx = 160; cy = 30; index = 7;
    elseif (k == 8)
      cx = 190; cy = 32; index = 8;
    elseif (k == 7)
      cx = 220; cy = 40; index = 9;
    elseif (k == 6)
      cx = 250; cy = 50; index = 10;
    elseif (k == 2)
      cx = 305; cy = 65; index = 11;
    elseif (k == 12)
      cx = 175; cy = 96; index = 12;
    elseif (k == 5)
      cx = 240; cy = 98; index = 13;
    end
  end

end
