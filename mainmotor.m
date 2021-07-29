close all;
clear all;
time = 0;
endtime = 300; % [sec]
global dt;
dt = 0.01; % [sec]
global dt2;
dt2=0.03;
% dt=0.01;
global r;
r=0.035;
global b;
b=0.075;
global m;
m=0.5;
global I;
I=0.000469;
global K1;
K1=1;
% K1=1;
global K2;
K2=20;
% K2=20;
global K3;
K3=0.01;
% K3=0.1;
global K4;
K4=15;
global K5;
K5=15;
global kx
kx=1;
global ky
ky=10;
global ktheta
ktheta=10;
global kdtheta_d 
kdtheta_d = 1;
global kdtheta_e 
kdtheta_e = 1;
global kptheta_d 
kptheta_d = 0.01;
global kptheta_e 
kptheta_e = 0.01;


global theta_roda
theta_roda=[0 0].';
global theta_roda_d
theta_roda_d=[0 0].';

nSteps = ceil((endtime - time)/dt);
result.time=[];
result.q=[];
result.q_r=[];
result.eta=[];
result.ukf=[];
result.eta_r=[];
result.e=[];
result.il=[];
result.ir=[];
result.measure=[];
result.WheelsSpeedReference=[];
result.xi_d=[];
result.dotxi_d=[];
result.xi=[];
result.xi_measured=[];
result.tau=[];
result.u=[];
result.einput=[];
result.wheelSpeedFromMotor=[];
q=[0 1 0].';
initial_state=[0 1 0 0 0 0 0].';
q_r=[0 0 0].';
e=[0 0 0].';
il=0;
ir=0;
xi_d=[0 0].';
xi_d_old=[0 0].';
dot_xi=[0 0].';
dot_xi_d=[0 0].';
xi=[0 0].';
tau=[0 0].';
eta_r=[0.3 0.0].';
q_measured = [0 0 0].';
xi_measured = [0 0 ].';
u=[0 0].';
x_e=initial_state;
counter = 0;
for i=1 : nSteps
    time = time + dt;
    counter = counter + 1;
    if(time>=100 && time<140)
        il=0.1;
        ir=0;
    elseif(time>=0 && time<60)
        il=0;
        ir=0.2;         
    elseif(time>=230 && time<260)
        il=0.1;
        ir=0;          
    else
        il=0;
        ir=0;
    end
   
    %Start reference trajectory block 
    dot_q_r = Calculate_dot_q_r(eta_r,q_r);
    q_r=dot_q_r*dt+q_r;
    
    %Mainting 30ms of control loop
   if(counter == 3) 
    
     %Start UKF
    [q_measured,xi_measured] = Observation(q, xi);
    [x_e,Pe]=unscentedkalmanfilter(dt2,initial_state,q_measured,xi_measured,tau);   
    %End UKF

    
    %Start tracking error block
    e=Calculate_e(q_r,x_e(1:3));
    %End tracking error block
    
    %Start auxiliary velocity block
    eta =   cinematicController(e,eta_r);
    %End auxiliary velocity block
    
    %Start Desired Velocity block 
    xi_d_old=xi_d;

    xi_d=ComputeWheelsSpeed(eta,x_e(6),x_e(7));%xi_d
    %End Desired Velocity block
    
    %Start Backstepping block
    dot_xi_d = (xi_d - xi_d_old)/dt2;
    u = PD(x_e(4:5),xi_d);
    %Start control input block

    tau = Control_Input(x_e(1:3),u,x_e(6),x_e(7));
    %End control input block
    counter =0;
   end
   
    [wheelSpeedFromMotor,e_motor] = tau2tensao(tau,xi_d);
    dot_xi = Dynamic(q,tau,il,ir);
    xi=dot_xi*dt+xi;
    %End Dynamic input block
    
    %Start kinematic model
    dot_q = Calculate_dot_q(xi,q,il,ir);
    q=dot_q*dt+q;

    %Start kinematic model
    result.ukf=[result.ukf x_e];
    result.xi_d=[result.xi_d,xi_d];
    result.xi=[result.xi,xi];
    result.eta_r=[result.eta_r,eta_r];
    result.time=[result.time,time];
    result.q=[result.q,q];
    result.q_r=[result.q_r,q_r];
    result.e=[result.e,e];
    result.u=[result.u,u];
    result.il=[result.il,il];
    result.ir=[result.ir,ir];
    result.einput=[result.einput,e_motor];
    result.dotxi_d=[result.dotxi_d,dot_xi_d];
    result.wheelSpeedFromMotor=[result.wheelSpeedFromMotor,wheelSpeedFromMotor];
    result.measure=[result.measure,q_measured];
    result.xi_measured=[result.xi_measured,xi_measured];
    result.tau=[result.tau,tau];
end


fig1=figure;
hold on
plot(result.q_r(1,:),result.q_r(2,:),'b','DisplayName','$Referencia$');
plot(result.q(1,:),result.q(2,:),'r','DisplayName','$Controlador$');
 plot(result.ukf(1,:),result.ukf(2,:),'g','DisplayName','$UKF$');
xlabel('$$\rho_x$$[m]','Interpreter','Latex','FontSize',20) ;
ylabel('$$\rho_y$$[m]','Interpreter','Latex','FontSize',20) ;
hl = legend('Location','none');
set(hl, 'Interpreter','latex');
grid on
hl.FontSize = 20;
hold off
saveas(fig1,'trajeto','fig');


fig2 = figure;
hold on
plot(result.time(1,:),result.einput(1,:),'r','DisplayName','$u_l$');
plot(result.time(1,:),result.einput(2,:),'g','DisplayName','$u_r$');
xlabel('Tempo[s]','Interpreter','Latex','FontSize',20) ;
ylabel('Tensao[V]','Interpreter','Latex','FontSize',20) ;
hl = legend('Location','best');
set(hl, 'Interpreter','latex');
hl.FontSize = 20;
grid on
hold off
saveas(fig2,'tensao','fig');






fig3 = figure;
hold on
plot(result.time(1,:),result.e(1,:),'r','DisplayName','$e_x $');
plot(result.time(1,:),result.e(2,:),'g','DisplayName','$e_y$ ');
xlabel('Tempo[s]','Interpreter','Latex','FontSize',20) ;
ylabel('Erro[m]','Interpreter','Latex','FontSize',20) ;
hl = legend('Location','best');
set(hl, 'Interpreter','latex');
hl.FontSize = 20;
grid on
hold off
saveas(fig3,'poserror','fig');


fig4 = figure;
hold on
plot(result.time(1,:),result.e(3,:),'b','DisplayName','$e_{\theta}$');
xlabel('Tempo$$[s]$$','Interpreter','Latex','FontSize',20) ;
ylabel('Erro[rad]','Interpreter','Latex','FontSize',20) ;
hl = legend('Location','best');
set(hl, 'Interpreter','latex');
hl.FontSize = 20;
grid on
hold off
saveas(fig4,'thetarror','fig');


fig5=figure;
hold on
plot(result.time(1,:),result.il(1,:),'red--','DisplayName','$i_{l}$');
plot(result.time(1,:),result.ukf(6,:),'g','DisplayName','$\hat{i}_{l}$');
xlabel('Tempo$$[s]$$','Interpreter','Latex','FontSize',20) 
ylabel('$$i_l$$','Interpreter','Latex','FontSize',30) 
hl = legend('Location','best');
set(hl, 'Interpreter','latex');
hl.FontSize = 20;
grid on
hold off
saveas(fig5,'il','fig');
% 
fig6=figure;
hold on
plot(result.time(1,:),result.ir(1,:),'red--','DisplayName','$i_{r}$');
plot(result.time(1,:),result.ukf(7,:),'g','DisplayName','$\hat{i}_{r}$');
xlabel('Tempo$$[s]$$','Interpreter','Latex','FontSize',20) ;
ylabel('$$i_r$$','Interpreter','Latex','FontSize',30) ;
hl = legend('Location','best');
set(hl, 'Interpreter','latex');
hl.FontSize = 20;
grid on
hold off
saveas(fig6,'ir','fig');

fig7 = figure;
hold on
plt1=plot(result.time(1,:),result.tau(1,:),'g','DisplayName','$\tau_l$');
plot(result.time(1,:),result.tau(2,:),'red','DisplayName','$\tau_r$');
xlabel('Tempo$$[s]$$','Interpreter','Latex','FontSize',20) ;
ylabel('$Torque$','Interpreter','Latex','FontSize',20) ;
hl = legend('show');
set(hl, 'Interpreter','latex');
hl.FontSize = 20;
grid on
hold off
saveas(fig7,'torque','fig');


%Calculo RMSE%

RMSEukf= sqrt(mean((result.ukf(1,:) - result.q(1,:)).^2));
RMSEmeasure= sqrt(mean((result.measure(1,:) -result.q(1,:)).^2));
fprintf("###############################################\n");
fprintf("RMSE UKF POS X: %d \n",RMSEukf);
fprintf("RMSE MEASURE POS X: %d \n",RMSEmeasure);

RMSEukf= sqrt(mean((result.ukf(2,:) - result.q(2,:)).^2));
RMSEmeasure= sqrt(mean((result.measure(2,:) -result.q(2,:)).^2));
fprintf("###############################################\n");
fprintf("RMSE UKF POS Y: %d \n",RMSEukf);
fprintf("RMSE MEASURE POS Y: %d \n",RMSEmeasure);


RMSEukf= sqrt(mean((result.ukf(3,:) - result.q(3,:)).^2));
RMSEmeasure= sqrt(mean((result.measure(3,:) -result.q(3,:)).^2));
fprintf("###############################################\n");
fprintf("RMSE UKF POS Theta: %d \n",RMSEukf);
fprintf("RMSE MEASURE POS Theta: %d \n",RMSEmeasure);


RMSEukf= sqrt(mean((result.ukf(4,:) - result.xi(1,:)).^2));
RMSEmeasure= sqrt(mean((result.xi_measured(1,:) -result.xi(1,:)).^2));
fprintf("###############################################\n");
fprintf("RMSE UKF Omega l: %d \n",RMSEukf);
fprintf("RMSE MEASURE Omega l: %d \n",RMSEmeasure);

RMSEukf= sqrt(mean((result.ukf(5,:) - result.xi(2,:)).^2));
RMSEmeasure= sqrt(mean((result.xi_measured(2,:) -result.xi(2,:)).^2));
fprintf("###############################################\n");
fprintf("RMSE UKF Omega r: %d \n",RMSEukf);
fprintf("RMSE MEASURE Omega r: %d \n",RMSEmeasure);


RMSEukf= sqrt(mean((result.ukf(6,:) - result.ir(1,:)).^2));
fprintf("###############################################\n");
fprintf("RMSE UKF i_r: %d \n",RMSEukf);

RMSEukf= sqrt(mean((result.ukf(7,:) - result.il(1,:)).^2));
fprintf("###############################################\n");
fprintf("RMSE UKF i_l: %d \n",RMSEukf);

function e=Calculate_e(q_r,q)%Eq 20 TCC

Aux=[cos(q(3)) sin(q(3)) 0;
    -sin(q(3)) cos(q(3)) 0;
    0               0      1];
e =Aux*(q_r-q);
end


function Eta = cinematicController(Pe,Eta_r)%Eq 22 TCC
global K1;
global K2;
global K3;

omega=Eta_r(2)+(Eta_r(1)/2)*(K2*(Pe(2)+K3*Pe(3))+(1/K3)*sin(Pe(3)));

if(omega>pi)
    omega=pi;
elseif(omega<-pi)
    omega=-pi;
end

v=K1*Pe(1)+Eta_r(1)*cos(Pe(3))-K3*Pe(3)*omega;

if(v>1.3)
    v=1.3;
elseif(v<-1.3)
    v=-1.3;
end

Eta=[v omega].';  
end

function wheelsSpeed = ComputeWheelsSpeed(Eta,ile,ire)%Eq 19 TCC
global r;
global b;
wheelsSpeed=(1/(2*r))*[(2*((1-ile)^-1)) (-b*((1-ile)^-1))
                        (2*((1-ire)^-1)) (b*((1-ire)^-1))];
wheelsSpeed=wheelsSpeed*Eta;

if(wheelsSpeed(1)>34.5575)
    wheelsSpeed(1)=34.5575;
end
if(wheelsSpeed(1)<-34.5575)
    wheelsSpeed(1)=-34.5575;
end
if(wheelsSpeed(2)>34.5575)
    wheelsSpeed(2)=34.5575;
end
if(wheelsSpeed(2)<-34.5575)
    wheelsSpeed(2)=-34.5575;
end

end



%Eq (33) e (34) TCC
function [wheelSpeed,tensao] = tau2tensao(tau,xi)
R=1.3/1.5;
n=100;
kt=0.1569/(100*1.5);
ke = 6/(100*(34.55752));
uL = (R*tau(1)/(2*n*kt))+n*ke*xi(1);
uR = (R*tau(2)/(2*n*kt))+n*ke*xi(2);
tensao = [uL uR].';
wl = (1/(n*ke))*(uL - ((R/(2*kt*n))*tau(1)));
wr = (1/(n*ke))*(uR - ((R/(2*kt*n))*tau(2)));
wheelSpeed = [wl wr].';
end

function tau = Control_Input(q,u,ile,ire) %Equation 4
global r;
global b;
global m;
global I;
S=(1/(2*b))*[b*r*(1-ile)*cos(q(3)) b*r*(1-ire)*cos(q(3))
            b*r*(1-ile)*sin(q(3))  b*r*(1-ire)*sin(q(3))
               -2*r*(1-ile)            2*r*(1-ire)];  
B = [cos(q(3)) cos(q(3))
    sin(q(3))   sin(q(3))
    -b/2        b/2];
M = [m 0 0
     0 m 0
     0 0 I];
tau=((((S.')*B))^-1)*((S.')*M*S)*u;       


end




function dot_xi = Dynamic(q,tau,il,ir) %Eq  5 TCC
global r;
global b;
global m;
global I;

S=(1/(2*b))*[b*r*(1-il)*cos(q(3)) b*r*(1-ir)*cos(q(3))
            b*r*(1-il)*sin(q(3))  b*r*(1-ir)*sin(q(3))
               -2*r*(1-il)            2*r*(1-ir)];  
B = [cos(q(3)) cos(q(3))
    sin(q(3))   sin(q(3))
    -b/2        b/2];
M = [m 0 0
     0 m 0
     0 0 I];


dot_xi=((S')*B)*tau;

dot_xi=((S')*M*S)^(-1)*dot_xi;



end



function u = PD(xi,xi_d)%Eq 23 TCC
global dt;
global theta_roda;
global theta_roda_d;
global kdtheta_d; 
global kdtheta_e ;
global kptheta_d;
global kptheta_e ;
KPD = [kdtheta_e 0 kptheta_e 0
       0       kdtheta_d  0 kptheta_d];
theta_roda   = theta_roda+xi*dt;
theta_roda_d = theta_roda_d+xi_d*dt;
xtiu=[(xi(1,1)-xi_d(1,1)) (xi(2,1)-xi_d(2,1)) (theta_roda(1,1)-theta_roda_d(1,1)) (theta_roda(2,1)-theta_roda_d(2,1))].';

u=-KPD*xtiu;
            
end


function dot_q_r = Calculate_dot_q_r(eta_r,q_r)
S_r=[cos(q_r(3)) 0
     sin(q_r(3)) 0
     0           1];
dot_q_r=S_r*eta_r;
end


function dot_q = Calculate_dot_q(xi,q,il,ir)%Eq 3 TCC
global r;
global b;
S=(1/(2*b))*[b*r*(1-il)*cos(q(3)) b*r*(1-ir)*cos(q(3))
            b*r*(1-il)*sin(q(3))  b*r*(1-ir)*sin(q(3))
               -2*r*(1-il)            2*r*(1-ir)];  
dot_q=S*xi;
end


%Simulate observation noise 
function [q_measured,xi_measured] = Observation(q, xi)
Rsigma_q=diag([0.01 0.01 0.04]).^2; 
Rsigma_xi=diag([0.001 0.001]).^2;
% Rsigma_q=diag([0.02 0.015 0.01]); 
% Rsigma_xi=diag([0.0015 0.0005]);
q_measured=q+Rsigma_q*randn(3,1);
xi_measured=xi+Rsigma_xi*randn(2,1);
end

