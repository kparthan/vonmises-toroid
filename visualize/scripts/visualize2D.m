function [] = visualize2D(K)

  hold on;
  %c(1,:) = [1 0.5 0];

  % plot the sampled data
  for k = 1:K
   %k
   data_file = strcat('../sampled_data/comp',num2str(k),'.dat');
   M = load(data_file);

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

   colors = rand(1,3);
   %plot(angles(:,1),angles(:,2),'.','Color',c(k,:));
   plot(angles(:,1),angles(:,2),'.','Color',colors);
  end  

  % create legend
  %N = [1:K];
  %legend_cell = [cellstr(num2str(N','%d'))];
  %legend(legend_cell);

end
