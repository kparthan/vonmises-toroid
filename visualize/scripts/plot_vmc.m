function [] = plot_vmc(data_file)

% draw a unit circle
x = [-1:0.01:1];
y = sqrt(1 - x.^2);
plot(x,y,'.','Color',[0.5 0.5 0.5]);
hold on;
x2 = x;
plot(x2,-y,'.','Color',[0.5 0.5 0.5]);

%hold on;

% plot the sampled data
M = load(data_file);
x = M(:,1);
y = M(:,2);
%colors = rand(1,3);
plot(x,y,'.','Color',[1 0 0]);

xlabel('X');
ylabel('Y');

end