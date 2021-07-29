function [x_estimated,P_estimated] = unscentedkalmanfilter(dt,initial_state,q_measured,xi_measured,tau)
persistent firstRun
persistent xcap P
persistent Q R 
persistent Wm Wc
persistent alpha beta kappa
persistent gamma 
persistent nse
if isempty(firstRun)%Realiza alguns cálculos na primeiro loop
% Ruido do processo
Q=0.001*eye(7);
% Q(4,4)=0.002;
% Q(5,5)=0.002;
% Q(5,5)=0.00005;
% Ruido de medição
%R=diag([0.3 0.3 toRadian(3) 0.015 0.015]).^2;
%R=diag([0.02 0.015 0.01 0.0015 0.0005]);
R = diag([0.013 0.015 0.04 0.0015 0.0012]).^2;

xEst=initial_state;

% Parametros dos UKF
alpha=1;
beta =2;
kappa=0;

n=length(xEst);
lamda=alpha^2*(n+kappa)-n;

%Cálculo dos pesos
Wm=[lamda/(lamda+n)];
Wc=[(lamda/(lamda+n))+(1-alpha^2+beta)];
for i=1:2*n
    Wm=[Wm 1/(2*(n+lamda))];
    Wc=[Wc 1/(2*(n+lamda))];
end
gamma=sqrt(n+lamda);

P = 0.1*eye(7);
xcap= initial_state;%+0.001*eye(7)*randn(7,1);
firstRun=1;
nse=7;
end

X                              = sigmapoints(xcap,P,gamma); %Calculo dos sigma points
[x_pre,tX,P_pre,Xdev]   = utf(X,Wm,Wc,nse,Q,dt,tau); %Propagração dos sigma points em f[k]
[y_pre,~,Pzz,Xdev1] = uth(tX,Wm,Wc,nse,R); %Propagação dos sigma points em h[k]

Pxz  = Xdev'*diag(Wc)*Xdev1;   %Covariância x-z
K    = Pxz/Pzz;     %Ganho de Kalman
measures=[q_measured
          xi_measured];%Medidas
%y_pre(3)=wrapToPi(y_pre(3));
xcap = x_pre+(K*(measures-y_pre'))'; %Estimativa do estado em k
P    = P_pre-K*Pxz' ; %Estimativa da covariância em k
%xcap(3)=wrapToPi(xcap(3));
%Salva os valores
x_estimated = xcap';
P_estimated= P;
end


function [tmean,tsigmapoints,Pzz,tdeviations]=uth(X,Wm,Wc,n,R)
%Unscented Transformation
L     = size(X,1);
tmean = zeros(1,5);

Y     = X(:,1:5);
tmean = tmean+Wm*Y;%Predição das medidas
% tmean(3)=wrapToPi(tmean(3));
tsigmapoints = Y;
tdeviations  = Y-tmean(ones(1,L),:);
Pzz=tdeviations'*diag(Wc)*tdeviations+R;%Covariância z-z
end

function [tmean,tsigmapoints,P_pre,tdeviations]=utf(X,Wm,Wc,n,Q,dt,u)
%Unscented Transformation
L     = size(X,1);
tmean = zeros(1,n);

for k=1:L 
    X(k,:)    = f(X(k,:), u,dt); 
    tmean     = tmean+Wm(k)*X(k,:);       
end

tsigmapoints = X;
tdeviations  = X-tmean(ones(1,L),:);
P_pre=tdeviations'*diag(Wc)*tdeviations+Q;%Predição da covariância
end

function X=sigmapoints(x,P,gamma)
xaux=[x(1),x(2),x(3),x(4),x(5),x(6),x(7)];
A = gamma*chol(P)';
Y = xaux(ones(1,numel(xaux)),:);
X = [xaux; Y+A; Y-A]; 
end



function x = f(x, u,dt)
r=0.035;
b=0.075;
m=0.5;
I=0.000469;
% x(6)=0;
% x(7)=0;
il=x(6);
ir=x(7);
% x(3)=wrapToPi(x(3));

S=(1/(2*b))*[b*r*(1-il)*cos(x(3)) b*r*(1-ir)*cos(x(3))
            b*r*(1-il)*sin(x(3))  b*r*(1-ir)*sin(x(3))
               -2*r*(1-il)            2*r*(1-ir)]; 
           
B = [cos(x(3)) cos(x(3))
    sin(x(3))   sin(x(3))
    -b/2        b/2];

M = [m 0 0
     0 m 0
     0 0 I];


dot_xi=((S')*B)*u;

dot_xi=((S')*M*S)^(-1)*dot_xi; 


xi=[x(4)
    x(5)];

xi=xi+dot_xi*dt;

dot_q=S*xi;

q=[x(1)
   x(2)
   x(3)];
q=q+dot_q*dt;



x= [q
    xi
    x(6)
    x(7)];
end

function radian = toRadian(degree)
% degree to radian
radian = degree/180*pi;
end