clc;
clear;

%%defining the mass matrix and the stiffness matrix
L_1 = 1;
L_2 = 1;
L_3 = 1;
m_1 = 1;
m_2 = 1;
m_3 = 1;
g=9.8;
zeta=0.05;
nmodes=3;
theta_dot0=[1;0;0];
theta0=[0;0;0];
n=10;

M=[(m_1+m_2+m_3)*(L_1)^2 m_2*L_1*L_2+m_3*L_1*L_2 m_3*L_1*L_3; m_2*L_1*L_2+m_3*L_1*L_2 (m_2+m_3)*(L_2)^2 m_3*L_3*L_2; m_3*L_3*L_1 m_3*L_3*L_2 m_3*(L_3)^2];
K=[(m_1+m_2+m_3)*g*L_1 0 0; 0 (m_2+m_3)*g*L_2 0; 0 0 m_3*g*L_3];
%%
%%natural frequencies and modeshapes:
[phi, evals] = eig(K,M);
omega=sqrt(diag(evals));
fprintf ('the first natural frequency of the system in rad/sec are:%d\n',omega(1));
fprintf ('the second natural frequency of the system in rad/sec are:%d\n',omega(2));
fprintf ('the third natural frequency of the system in rad/sec are:%d\n',omega(3));
fprintf ('the first modeshape of the system in rad/sec are:');
fprintf(' %d\n ', phi(:,1).')
fprintf ('the second modeshape of the system in rad/sec are:');
fprintf(' %d\n ', phi(:,2)')
fprintf ('the third modeshape of the system in rad/sec are:');
fprintf(' %d\n ', phi(:,3).')
y=[zeros(1,nmodes);phi];
x=1:1:4;
%%plotting modeshapes
figure;
plot(x,y);
title('mode shapes')
xlabel('n')
ylabel('mode shapes')
hold on

%modeshapes animation
%%animation of first modeshape
figure;
   for i=1:4
    plot(x(i),y(i),'x');
    hold on
    plot(x(1:i),y(1:i));
    axis([1 4 -1.5 1]);
    grid on
    pause(0.2);
    drawnow
   end
   hold on
    for i=1:4
    plot(x(i),y(i+4),'x');
    hold on
    plot(x(1:i),y(5:i+4));
    axis([1 4 -1.5 1]);
    grid on
    pause(0.2);
    drawnow
    end
    hold on
   for i=1:4
    plot(x(i),y(i+8),'x');
    hold on
    plot(x(1:i),y(9:i+8));
    axis([1 4 -1.5 1]);
    grid on
    title('modeshapes')
    xlabel('phi')
    ylabel('modeshapes')
    pause(0.2);
    drawnow
   end
%%
%%getting the eigenvectors, and solving the equation using modal analysis
%%model
syms t
Kbar=(inv(sqrtm(M)))*(K)*(inv(sqrtm(M)));
[V, D]=eig(Kbar); 
omegga=sqrt(diag(D));
omegadamped=omegga*sqrt(1-zeta);
P=V/norm(V);  %%P is the normalized eigenvectors
S= (inv(sqrtm(M)))*P;  Sinv=P'*sqrtm(M);  %%getting S function
r_dot0=Sinv*theta_dot0;   %%getting r(t)
r0=Sinv*theta0;
r1= (sqrt(((zeta*omegga(1)*r0(1)+r_dot0(1))^2+(r0(1)*omegadamped(1))^2)/omegadamped(1)^2))+sin(omegadamped(1)*t+atan(omegadamped(1)*r0(1)/(r_dot0(1)+zeta*omegga(1)*r0(1))));
r2= (sqrt(((zeta*omegga(2)*r0(1)+r_dot0(2))^2+(r0(2)*omegadamped(2))^2)/omegadamped(2)^2))+sin(omegadamped(2)*t+atan(omegadamped(2)*r0(2)/r_dot0(2)+zeta*omegga(2)*r0(2)));
r3= (sqrt(((zeta*omegga(3)*r0(3)+r_dot0(3))^2+(r0(3)*omegadamped(3))^2)/omegadamped(3)^2))+sin(omegadamped(3)*t+atan(omegadamped(3)*r0(3)/(r_dot0(3)+zeta*omegga(3)*r0(3))));
x=S*[r1;r2;r3];
%%
%%plotting the response
figure;
fplot(t,x(1));
title('theta 1 vs t')
xlabel('t')
ylabel('theta 1(t)')
hold on
figure;
fplot(t,x(2));
title('theta2(t) vs t')
xlabel('t')
ylabel('theta 2(t)')
hold on
figure;
fplot(t,x(3));
title('theta3(t) vs t')
xlabel('t')
ylabel('theta 3(t)')
figure;
fplot(t,x(1));
hold on
fplot(t,x(2));
hold on
fplot(t,x(3));
title('theta1,theta2,theta3(t) vs t')
xlabel('t')
ylabel('theta1(t),theta2(t),theta3(t)')
%%
%%getting the forced response due to base excitation;
omega_b=1;
Y=0.03;
f1=2*zeta*omega(1)*omega_b*Y*cos(omega_b*t)+(omega(1))^2*Y*sin(omega_b*t);
F=P'*inv(sqrtm(M))*[f1;0;0];
%%particular solution

rp1=((omegga(1))^2*Y)/sqrt((omegga(1))^2-(omega_b)^2)*cos(omega_b*t-(pi/2));
rp2=((omegga(2))^2*Y)/sqrt((omegga(2))^2-(omega_b)^2)*cos(omega_b*t-(pi/2));
rp3=((omegga(3))^2*Y)/sqrt((omegga(3))^2-(omega_b)^2)*cos(omega_b*t-(pi/2));
xb=S*[rp1;rp2;rp3];
figure;
fplot(t,xb(1));
title('response of bottom due to base exitation')
xlabel('t')
ylabel('diplacement')
figure;
fplot(t,xb(2));
title('response of top due to base exitation')
xlabel('t')
ylabel('diplacement')
figure;
fplot(t,xb(3));
title('response of middle due to base exitation')
xlabel('t')
ylabel('diplacement')


 