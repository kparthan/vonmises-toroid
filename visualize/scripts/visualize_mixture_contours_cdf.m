function [] = visualize_mixture_contours_cdf(K,pdf)

  addpath('export_fig');

  bins_folder = '';
  outfile = '';
  % IND    -- 0
  % SINE   -- 1
  % COSINE -- 2
  pdf_type = -1;   
  if (strcmp(pdf,'ind') == 1)
    pdf_type = 0;
    bins_folder = '../sampled_data/bins_ind/';
    outfile = 'bvm_ind_class_b';
%    outfile = 'bvm_ind_class_b_actual';
  elseif (strcmp(pdf,'sine') == 1)
    pdf_type = 1;
    bins_folder = '../sampled_data/bins_sine/';
    outfile = 'bvm_sine_class_b';
%    outfile = 'bvm_sine_class_b_actual';
  elseif (strcmp(pdf,'cosine') == 1)
    pdf_type = 2;
    bins_folder = '../sampled_data/bins_cosine/';
    outfile = 'bvm_cosine_class_b';
%    outfile = 'bvm_cosine_class_b_actual';
  end
  %disp(strcat('bins_folder: ',{' '},bins_folder)); 
  %disp(strcat('pdf type: ',num2str(pdf_type)));
 
  % figure properties
  fig = figure();
  caxis([0 1]); % colormap range
  axis equal
  hold on;
  set(gcf, 'Color', 'w');
  %xlabel('\psi_1','fontsize',15);
  %ylabel('\psi_2','fontsize',15);
  xlabel('\phi','fontsize',15);
  ylabel('\psi','fontsize',15);
  xlabh = get(gca,'XLabel');
  ylabh = get(gca,'YLabel');
  set(xlabh,'interpreter','tex');
  set(ylabh,'interpreter','tex');

 % set(gca,'Xlim',[0 360]);
 % set(gca,'Ylim',[0 360]);
 % set(gca,'xtick',[0:60:360],'fontsize',10);
 % set(gca,'ytick',[0:60:360],'fontsize',10);
 % phi = 0:1:359.9;  % meshgrid columns (X-axis)
 % psi = 0:1:359.9;  % meshgrid rows (Y-axis)

  set(gca,'Xlim',[-180 180]);
  set(gca,'Ylim',[-180 180]);
  set(gca,'xtick',[-180:60:180],'fontsize',10);
  set(gca,'ytick',[-180:60:180],'fontsize',10);
  phi = -180:1:179.9;  % meshgrid columns (X-axis)
  psi = -180:1:179.9;  % meshgrid rows (Y-axis)

  % plot the contours 
  for k = 1:K
    %k
    data_file = strcat(bins_folder,'comp',num2str(k),'_prob_bins2D.dat');
    Mt = load(data_file);
    Mt_trans = Mt';
    prob_bins = Mt_trans;

    %Mt_rearranged = rearrange(Mt);  % [-180 ... 180]
    %prob_bins = Mt_rearranged';

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
    line_width = 1.0;
    if (k==14 || k == 15 || k==18 || k==19 || k==20 || k==23 || k==24 || k==28)  % IND
    %if (k==3 || k==4)  % SINE
      line_width = 1.5;
    end 
    [C,h] = contour(phi,psi,cdf_bins,contour_levels,'LineWidth',line_width,'LineColor','black');
    %[C,h] = contour(cdf_bins,1,'LineWidth',2,'LineColor','black');
    %clabel(C,h);

%    [row col] = ind2sub(size(prob_bins),max_index);
%    cx = phi(col);
%    cy = psi(row);
%    ht = text(cx,cy,num2str(k),'Color','red');
    [cx,cy,index] = number_component(k,pdf_type);
    if (index > 0)
      ht = text(cx,cy,num2str(index),'Color','red','fontsize',10);
    end

%    hcl = clabel(C,'Color','red');
%    for i=2:2:length(hcl)
%      old_label = get(hcl(i),'String');
%      new_label = num2str(k);
%      set(hcl(i),'String',new_label);
%    end
  end  

  % plot the entire mixture density
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
  for i = 1:n
    if angles(i,1) > 180
      angles(i,1) = angles(i,1) - 360;
    end
    if angles(i,2) > 180
      angles(i,2) = angles(i,2) - 360;
    end
  end

  hs = scatter3(angles(:,1),angles(:,2),norm_density,0.1,'cdata',norm_density);

  %colorbar
  output_fig = strcat('../figs/protein_modelling/',outfile,'.fig');
  output_eps = strcat('../figs/protein_modelling/',outfile,'.eps');
  output_pdf = strcat('../figs/protein_modelling/',outfile,'.pdf');

  %saveas(gcf,output_fig);
  %export_fig(output_pdf,'-pdf');
  %print2eps(output_eps);
  %eps2pdf(output_eps,output_pdf);

end

function [cx,cy,index] = number_component(k,pdf_type)

  cx = 1; cy = 1; index = -1;

  if (pdf_type == 0)  % IND
    if (k == 27)
      cx = -179; cy = 145; index = 1;
    elseif (k == 26)
      cx = -152; cy = 125; index = 2;
    elseif (k == 25)
      cx = -148; cy = 110; index = 3;
    elseif (k == 24)
      cx = -40; cy = 145; index = 4;
    elseif (k == 23)
      cx = -50; cy = 158; index = 5;
    elseif (k == 21)
      cx = -65; cy = 178; index = 6;
    elseif (k == 28)
      cx = -110; cy = 178; index = 7;
    elseif (k == 12)
      cx = -163; cy = 60; index = 8;
    elseif (k == 9)
      cx = -130; cy = 75; index = 9;
    elseif (k == 22)
      cx = -70; cy = 60; index = 10;
    elseif (k == 10)
      cx = -23; cy = 103; index = 11;
    elseif (k == 13)
      cx = -150; cy = 28; index = 12;
    elseif (k == 14)
      cx = -83; cy = 15; index = 13;
    elseif (k == 15)
      cx = -66; cy = 5; index = 14;
    elseif (k == 19)
      cx = -55; cy = -10; index = 15;
    elseif (k == 18)
      cx = -50; cy = -25; index = 16;
    elseif (k == 20)
      cx = -50; cy = -50; index = 17;
    elseif (k == 16)
      cx = -50; cy = -73; index = 18;
    elseif (k == 17)
      cx = -100; cy = -68; index = 19;
    elseif (k == 11)
      cx = 12; cy = 65; index = 20;
    elseif (k == 4)
      cx = 30; cy = 53; index = 21;
    elseif (k == 2)
      cx = 73; cy = 35; index = 22;
    elseif (k == 3)
      cx = 50; cy = -3; index = 23;
    elseif (k == 5)
      cx = 65; cy = -20; index = 24;
    elseif (k == 6)
      cx = 120; cy = -50; index = 25;
    end
  elseif (pdf_type == 1)  % SINE
    if (k == 16)
      cx = -175; cy = 135; index = 1;
    elseif (k == 15)
      cx = -145; cy = 105; index = 2;
    elseif (k == 21)
      cx = -120; cy = 97; index = 3;
    elseif (k == 13)
      cx = -45; cy = 150; index = 4;
    elseif (k == 14)
      cx = -52; cy = 170; index = 5;
    elseif (k == 5)
      cx = -155; cy = 55; index = 6;
    elseif (k == 7)
      cx = -65; cy = 85; index = 7;
    elseif (k == 6)
      cx = -20; cy = 100; index = 8;
    elseif (k == 2)
      cx = -145; cy = 0; index = 9;
    elseif (k == 1)
      cx = -125; cy = 15; index = 10;
    elseif (k == 20)
      cx = -90; cy = -75; index = 11;
    elseif (k == 4)
      cx = -90; cy = -50; index = 12;
    elseif (k == 3)
      cx = -55; cy = -10; index = 13;
    elseif (k == 12)
      cx = 10; cy = 50; index = 14;
    elseif (k == 8)
      cx = 30; cy = 70; index = 15;
    elseif (k == 9)
      cx = 110; cy = -15; index = 16;
    end
  elseif (pdf_type == 2)  % COSINE
    cx = 1; cy = 1; index = -1;
  end

end

