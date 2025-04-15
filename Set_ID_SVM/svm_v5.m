function [w,b] = svm_v5(x,t,j)

x_cont = []; t_cont = []; xc = 0;
x_disc = []; t_disc = []; xd = 0;
if size(x,1) > size(x,2)
    x = x';
end
for i = 1:size(x,2)-2
    if j(i+1) == j(i)
        xc = xc + 1;
        x_cont(:,xc) = x(:,i);
        t_cont(xc) = t(i);
    elseif j(i+1) ~= j(i)
        xd = xd + 1;
        x_disc(:,xd) = x(:,i);
        t_disc(xd) = t(i);
    end
end

data1 = [x_cont(1,:); x_cont(2,:)]';
data2 = [x_disc(1,:); x_disc(2,:)]';

%{
figure(4);
plot(data1(:,1),data1(:,2),'r.','MarkerSize',15)
hold on
plot(data2(:,1),data2(:,2),'b.','MarkerSize',15)
%ezpolar(@(x)1);ezpolar(@(x)2);
axis equal
hold off
%}

data3 = [data1;data2];
theclass = ones(xc+xd,1);
theclass(1:xc) = -1;

cl = fitcsvm(data3,theclass,'KernelFunction','linear','BoxConstraint',Inf,'ClassNames',[-1,1]);

d = 0.02;
[x1Grid,x2Grid] = meshgrid(min(data3(:,1)):d:max(data3(:,1)),min(data3(:,2)):d:max(data3(:,2)));
xGrid = [x1Grid(:),x2Grid(:)];
[~,scores] = predict(cl,xGrid);

figure(5);
h(1:2) = gscatter(data3(:,1),data3(:,2),theclass,'br','.');
hold on
%ezpolar(@(x)1);
h(3) = plot(data3(cl.IsSupportVector,1),data3(cl.IsSupportVector,2),'ko');
contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),[0 0],'k');
legend(h,{'Flow','Jumps','Support Vectors'});
axis equal
hold off

w = [cl.Beta(1); cl.Beta(2)];
b = cl.Bias/cl.Beta(2);
