function [] = visualize_mixture_contours_default(K)

  addpath('export_fig');

  bins_folder = '../sampled_data/bins_kent/';
  outfile = 'kent_mix';
  %bins_folder = '../sampled_data/bins_vmf/';
  %outfile = 'vmf_mix';
 
  % plot the entire mixture density
  data_file = strcat(bins_folder,'mixture_density.dat');
  M = load(data_file);
  density = M(:,4);
  min1 = min(density);
  max1 = max(density);
  range1 = max1 - min1;
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
  hs = scatter3(angles(:,1),angles(:,2),density,1,'cdata',density);

  hold on;

  % figure properties
  set(gcf, 'Color', 'w');
  xlabel('Longitude','fontsize',20);
  ylabel('Co-latitude','fontsize',20);
  set(gca,'Ylim',[0 120]);
  set(gca,'xtick',[0:60:360],'fontsize',12);
  set(gca,'ytick',[0:30:120],'fontsize',12);
  view ([0 90]);

  % plot the contours 
  theta = 0:1:179.9;
  phi = 0:1:359.9;
  for k = 1:K
    %k
    data_file = strcat(bins_folder,'comp',num2str(k),'_prob_bins2D.dat');
    prob_bins2 = load(data_file);
    min2 = min(prob_bins2(:));
    [max2 max_index] = max(prob_bins2(:));
    range2 = max2 - min2;

    factor = range1 / range2;
    prob_bins = (prob_bins2 - min2) * factor;
    prob_bins = prob_bins + min1;
    [C,h] = contour(prob_bins,1,'LineWidth',2,'LineColor','black');
    %clabel(C,h);

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

  %uistack(hs,'bottom');
  %uistack(ht,'top');

  output_fig = strcat('../figs/',outfile,'.fig');
  output_eps = strcat('../figs/',outfile,'.eps');
  output_pdf = strcat('../figs/',outfile,'.pdf');

  %saveas(gcf,output_fig);
  %print2eps(output_eps);
  %eps2pdf(output_eps,output_pdf,1);
  %export_fig(output_pdf,'-pdf');

end
