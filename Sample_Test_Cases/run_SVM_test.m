clc; clear all; close all;

%%{
data1 = [3 3 6 6; 1 -1 1 -1]'; % -1 category
data2 = [1 0 0 -1; 0 1 -1 0]'; % +1 category
%}

%{
rng(1);
r = sqrt(rand(100,1)); t = 2*pi*rand(100,1);
data1 = [r.*cos(t), r.*sin(t)];

r2 = sqrt(3*rand(100,1)+1); t2 = 2*pi*rand(100,1);
data2 = [r2.*cos(t2), r2.*sin(t2)];
%}

figure(1); hold on;
plot(data1(:,1),data1(:,2),'r.','MarkerSize',15);
plot(data2(:,1),data2(:,2),'b.','MarkerSize',15);
%ezpolar(@(x)1); ezpolar(@(x)2);
axis equal;

data3 = [data1; data2];
the_class = ones(length(data1) + length(data2), 1);
theclass(1:length(data1)) = -1;

cl = fitcsvm(data3,theclass,'KernelFunction','linear','BoxConstraint',Inf,'ClassNames',[-1,1]);

d = 0.02;
[x1Grid,x2Grid] = meshgrid(min(data3(:,1)):d:max(data3(:,1)),min(data3(:,2)):d:max(data3(:,2)));
xGrid = [x1Grid(:),x2Grid(:)];
[~,scores] = predict(cl,xGrid);

figure;
h(1:2) = gscatter(data3(:,1),data3(:,2),theclass,'rb','.');
hold on
%ezpolar(@(x)1);
h(3) = plot(data3(cl.IsSupportVector,1),data3(cl.IsSupportVector,2),'ko');
contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),[0 0],'k');
legend(h,{'-1','+1','Support Vectors'});
axis equal
hold off

cl2 = fitcsvm(data3,theclass,'KernelFunction','rbf');
[~,scores2] = predict(cl2,xGrid);

figure;
h(1:2) = gscatter(data3(:,1),data3(:,2),theclass,'rb','.');
hold on
%ezpolar(@(x)1);
h(3) = plot(data3(cl2.IsSupportVector,1),data3(cl2.IsSupportVector,2),'ko');
contour(x1Grid,x2Grid,reshape(scores2(:,2),size(x1Grid)),[0 0],'k');
legend(h,{'-1','+1','Support Vectors'});
axis equal
hold off
