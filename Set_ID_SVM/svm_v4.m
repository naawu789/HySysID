function [w, b, dlim] = svm_v4(x0,t,J)

x_cont = []; t_cont = []; xc = 0;
x_disc = []; t_disc = []; xd = 0;
for i = 1:size(x0,2)-2
    if J(i+1) == J(i)
        xc = xc + 1;
        x_cont(:,xc) = x0(:,i);
        t_cont(xc) = t(i);
    elseif J(i+1) ~= J(i)
        xd = xd + 1;
        x_disc(:,xd) = x0(:,i);
        t_disc(xd) = t(i);
    end
end

s1 = [x_cont(1,:); x_cont(2,:)];
s2 = [x_disc(1,:); x_disc(2,:)];

C = 10000000;
x = [s1 s2];
y = [-ones(1,size(s1,2)) ones(1,size(s2,2))];
n = size(x,2);
ClassA = find( y == +1 );
ClassB = find( y == -1 );

dlim = [min(x(1,ClassA)), max(x(1,ClassA)), min(x(1,ClassB)), max(x(1,ClassB));
    min(x(2,ClassA)), max(x(2,ClassA)), min(x(2,ClassB)), max(x(2,ClassB))];

H = zeros(n,n);
for i=1:n
    for k=i:n
        H(i,k) = y(i)*y(k)*x(:,i)'*x(:,k);
        H(k,i) = H(i,k);
    end
end
f = -ones(n,1);
Aeq=y;
beq=0;
lb=zeros(n,1);
ub=C*ones(n,1);
%lb = [];
%ub = [];

Alg{1}='trust-region-reflective';
Alg{2}='interior-point-convex';
options=optimset('Algorithm',Alg{2},...
    'Display','off',...
    'MaxIter',20);
alpha=quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options)';

AlmostZero=(abs(alpha)<max(abs(alpha))/1e3);
alpha(AlmostZero)=0;
S=find(alpha>0 & alpha<C);

%{
C2 = 10;
x_A = x(1,ClassA); x_B = x(1,ClassB);
x_Am = mink(x_A,C2); x_Bm = mink(x_B,C2);

[x_m2,SA] = intersect(x(1,:),x_Am);
[x_m2,SB] = intersect(x(1,:),x_Bm);
S = [SA; SB]';
%}
%%

%S=find(alpha>0);
%{
w=0;
for i=S
    w=w+alpha(i)*y(i)*x(:,i);
end
b=mean(y(S)-w'*x(:,S));
%}
%%{
s3 = [x(:,S); ones(size(S))];
A = zeros(size(s3,2));
for i = 1:size(s3,2)
    for k = 1:size(s3,2)
       A(i,k) = s3(:,i)'*s3(:,k);
    end
end
X = linsolve(A,y(S)');
w = s3(1:end-1,:)*X;
b = s3(end,:)*X;
%}

%[dH D] = HausdorffDist(s1',s2');

%% Plot Results
Line = @(x1,x2) w(1)*x1+w(2)*x2+b;
LineA = @(x1,x2) w(1)*x1+w(2)*x2+b+1;
LineB = @(x1,x2) w(1)*x1+w(2)*x2+b-1;

figure(5);
hold off;
plot(x(1,ClassA),x(2,ClassA),'ro');
hold on;
plot(x(1,ClassB),x(2,ClassB),'bs');

plot(x(1,S),x(2,S),'ko','MarkerSize',12);

x1min = min(x(1,:));
x1max = max(x(1,:));
x2min = min(x(2,:));
x2max = max(x(2,:));

handle = ezplot(Line,[x1min x1max x2min x2max]);
set(handle,'Color','k','LineWidth',2);

handleA = ezplot(LineA,[x1min x1max x2min x2max]);
set(handleA,'Color','k','LineWidth',1,'LineStyle',':');

handleB = ezplot(LineB,[x1min x1max x2min x2max]);
set(handleB,'Color','k','LineWidth',1,'LineStyle',':');

legend('Data \in D','Data \in C');
%title(['w = [' num2str(w(1)) ', ' num2str(w(2)) '], b = ' num2str(b)]);
title('SVM applied to state x_1 and x_2');
xlabel('x_1'); ylabel('x_2');

end
