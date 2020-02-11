close all
%each loop a different reactant activity is step-changed
for zz=1:4
    %rate constant
k1=0.02;
k2=0.04;
k3=0.3;
k4=2;
%activities/pressures of reactants
P_A=1;
P_D=1;
P_E=1;
P_F=1;

%matrix A and vector b of the coupled ODE
A=[-(k1*P_A+k2*P_D),-k1*P_A,-k1*P_A;k2*P_D,-k3*P_E,0;0,k3*P_E,-k4*P_F];
b=[k1*P_A;0;0];

%V=eigenvectors, D=eigenvalues of matrix A
[V,D]=eig(A);

%initial condition: clean surface
theta0=[0;0;0];

%particular solution to the homogeneous ODE
thetap=A\-b;

%solve for arbitrary constants
c=V\(theta0-thetap);

%each function is the coverage of A*, B*, C*, as a function of time
f=@(t)real(((c(1)*V(1,1).*exp(D(1,1)*t)+c(2)*V(1,2).*exp(D(2,2)*t)+c(3)*V(1,3).*exp(D(3,3)*t))+thetap(1)));
g=@(t)real(((c(1)*V(2,1).*exp(D(1,1)*t)+c(2)*V(2,2).*exp(D(2,2)*t)+c(3)*V(2,3).*exp(D(3,3)*t))+thetap(2)));
h=@(t)real(((c(1)*V(3,1).*exp(D(1,1)*t)+c(2)*V(3,2).*exp(D(2,2)*t)+c(3)*V(3,3).*exp(D(3,3)*t))+thetap(3)));

%rate is k4 x activity of F x coverage of C*
r=@(t) k4*P_F*h(t);

%plot rate
figure(3*(zz-1)+1)
l=fplot(r,[0,100],'k')
data_r0=r(100);
set(l,'LineWidth',3)

%plot coverages
figure(3*(zz-1)+2)
data_cov0=[f(100),g(100),h(100)];
l=fplot(f,[0,100],'k')
set(l,'LineWidth',3)
hold on
l=fplot(g,[0,100],'b')
set(l,'LineWidth',3)
l=fplot(@(t) h(t)*10,[0,100],'g')
set(l,'LineWidth',3)
theta0=[f(100);g(100);h(100)];

%the value of zz determines which activity is step-changed
if zz==1
    P_A=3;
end
if zz==2
    P_D=3;
end
if zz==3
    P_E=3;
end
if zz==4
    P_F=3;
end
figure(3*(zz-1)+1)
hold on

%solve for the transient rate after the step change
A=[-(k1*P_A+k2*P_D),-k1*P_A,-k1*P_A;k2*P_D,-k3*P_E,0;0,k3*P_E,-k4*P_F];
b=[k1*P_A;0;0];

[V,D]=eig(A);

thetap=A\-b;

c=V\(theta0-thetap);

f=@(t)real(((c(1)*V(1,1).*exp(D(1,1)*t)+c(2)*V(1,2).*exp(D(2,2)*t)+c(3)*V(1,3).*exp(D(3,3)*t))+thetap(1)));
g=@(t)real(((c(1)*V(2,1).*exp(D(1,1)*t)+c(2)*V(2,2).*exp(D(2,2)*t)+c(3)*V(2,3).*exp(D(3,3)*t))+thetap(2)));
h=@(t)real(((c(1)*V(3,1).*exp(D(1,1)*t)+c(2)*V(3,2).*exp(D(2,2)*t)+c(3)*V(3,3).*exp(D(3,3)*t))+thetap(3)));
r=@(t) k4*P_F*h(t);

%plot rate after the perturbation
l=fplot(@(t) r(t-100),[100,140],'k')

%t=time after perturbation
t=0:0.01:40;
%store rate after perturbation in vector data_r
data_r=[r(t)'];
set(l,'LineWidth',3)
axis([0,140,0,100])
axis 'auto y'

%plot coverages fter perturbation
figure(3*(zz-1)+2)
hold on
l=fplot(@(t) f(t-100),[100,140],'k')
set(l,'LineWidth',3)
l=fplot(@(t) g(t-100),[100,140],'b')
set(l,'LineWidth',3)
l=fplot(@(t) h(t-100)*10,[100,140],'g')
set(l,'LineWidth',3)
axis([0,140,0,100])
axis 'auto y'
t=0:0.01:40;
%store transient coverages in data_cov
data_cov=[f(t)',g(t)',h(t)'];

%for loop for calculating sensitivities/DoRC
for a=1:4
    if a==1
        k1=k1*1.0001;
    elseif a==2
        k2=k2*1.0001;
    elseif a==3
        k3=k3*1.0001;
    elseif a==4
        k4=k4*1.0001;
    end
%solve for rate with perturbed rate constants
A=[-(k1*P_A+k2*P_D),-k1*P_A,-k1*P_A;k2*P_D,-k3*P_E,0;0,k3*P_E,-k4*P_F];
b=[k1*P_A;0;0];
%when calculating the sensntivities, we do not want to change the
%eigenvalues. We store the new eigenvalues in D_new. We will not use D_new
%when calculating the sensitivities--instead we use the old eigenvalues,
%D>. This is equivalent to treating \lambda*t (\lambda=eigenvalues) as
%dimensionsless time \tau_j and not differentitating \tau_j when
%calculating the sensitivity.
V_old=V;
[V,D_new]=eig(A);
thetap=A\-b;
c=V\(theta0-thetap);
f=@(t)real(((c(1)*V(1,1).*exp(D(1,1)*t)+c(2)*V(1,2).*exp(D(2,2)*t)+c(3)*V(1,3).*exp(D(3,3)*t))+thetap(1)));
g=@(t)real(((c(1)*V(2,1).*exp(D(1,1)*t)+c(2)*V(2,2).*exp(D(2,2)*t)+c(3)*V(2,3).*exp(D(3,3)*t))+thetap(2)));
h=@(t)real(((c(1)*V(3,1).*exp(D(1,1)*t)+c(2)*V(3,2).*exp(D(2,2)*t)+c(3)*V(3,3).*exp(D(3,3)*t))+thetap(3)));
%store value of r as a function of time for each perturbed rate. Undo
%perturbation in rate constant.
if a==1
    r_1=@(t) k4*P_F*h(t);
    k1=k1/1.0001;
elseif a==2
    r_2=@(t) k4*P_F*h(t);
    k2=k2/1.0001;
elseif a==3
    r_3=@(t) k4*P_F*h(t);
    k3=k3/1.0001;
elseif a==4
    r_4=@(t) k4*P_F*h(t);
    k4=k4/1.0001;
end
end
%calculate xrc_i
xrc1=@(t) k1./r(t).*(r_1(t)-r(t))./(0.0001*k1);
xrc2=@(t) k2./r(t).*(r_2(t)-r(t))./(0.0001*k2);
xrc3=@(t) k3./r(t).*(r_3(t)-r(t))./(0.0001*k3);
xrc4=@(t) k4./r(t).*(r_4(t)-r(t))./(0.0001*k4);
xtot=@(t) xrc1(t)+xrc2(t)+xrc3(t)+xrc4(t);

%plot xrc_i, format plots
figure
w=fplot(@(t) xtot(t-100),[100,140])
set(w,'LineWidth',3)
hold on
w=fplot(@(t) xrc1(t-100),[100,140])
set(w,'LineWidth',3)
get(w,'Color')
w=fplot(@(t) xrc2(t-100),[100,140])
set(w,'LineWidth',3)
get(w,'Color')
w=fplot(@(t) xrc3(t-100),[100,140])
set(w,'LineWidth',3)
get(w,'Color')
w=fplot(@(t) xrc4(t-100),[100,140])
set(w,'LineWidth',3)
get(w,'Color')
t=0:0.01:40;
data_xrc=[xrc1(t)',xrc2(t)',xrc3(t)',xrc4(t)'];
figure(3*(zz-1)+1)
xlabel('Time / a.u.')
ylabel('Rate / a.u.')
set(gca,'linewidth',2)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 14)
figure(3*(zz-1)+2)
xlabel('Time / a.u.')
ylabel('Coverage')
set(gca,'linewidth',2)
legend('C','D','E')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 14)
figure(3*(zz-1)+3)
legend('sum X_{RC,i}','X_{RC,1}','X_{RC,2}','X_{RC,3}','X_{RC,4}')
xlabel('Time / a.u.')
ylabel('X_{RC,i}')
set(gca,'linewidth',2)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 14)
legend('sum X_{RC,i}','X_{RC,1}','X_{RC,2}','X_{RC,3}','X_{RC,4}')
end
