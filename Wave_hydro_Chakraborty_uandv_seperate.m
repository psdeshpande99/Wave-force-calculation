D = 1.25; 
disp(['Diameter of cylinder is ',num2str(D),'m'])
T = 10;
disp(['Time period is ',num2str(T),'sec'])
d = 14;
disp(['depth of sea water is ',num2str(d),'m'])
H = 6;
disp(['Wave Height is ',num2str(H),'m'])
z = [-0*d, -0.2*d, -0.4*d, -0.6*d, -0.8*d, -1*d];
theta = 0:1:360;
Lo =1.56*T^2;
disp(['Deep water wavelength is ',num2str(Lo),'m']) %%display all parameters and then
L = d/0.1319;
disp(['Wavelength (from wave tables) is ',num2str(L),'m'])
k = 2*pi/L;
w = 2*pi/T;
rho = 1025;
disp(['Density of sea water is ',num2str(rho),'kg/m^3'])
umax = (pi*H/T)*(cosh(k*(d+z(1)))/sinh(k*d)).*sind(90);
disp(['Surface velocity of wave is',num2str(umax),'m/sec'])
Ft1c = zeros(length(z),length(theta));
Ft2c = zeros(length(z),length(theta));
Ft3c = zeros(length(z),length(theta));
Ft4c = zeros(length(z),length(theta));
Ft5c = zeros(length(z),length(theta));
Fi1c = zeros(length(z),length(theta));
Fi2c = zeros(length(z),length(theta));
Fi3c = zeros(length(z),length(theta));
Fi4c = zeros(length(z),length(theta));
Fi5c = zeros(length(z),length(theta));
Fd1c = zeros(length(z),length(theta));
Fd2c = zeros(length(z),length(theta));
Fd3c = zeros(length(z),length(theta));
Fd4c = zeros(length(z),length(theta));
Fd5c = zeros(length(z),length(theta));
Mt1c = zeros(length(z),length(theta));
Mt2c = zeros(length(z),length(theta));
Mt3c = zeros(length(z),length(theta));
Mt4c = zeros(length(z),length(theta));
Mt5c = zeros(length(z),length(theta));
Cd1c = 1.5;Cd2c = 1.5;Cd3c = 1.35;Cd4c = 1.5;Cd5c = 1.35;
Cm1c = 1.25;Cm2c = 1.4; Cm3c = 1.65; Cm4c = 1.4; Cm5c = 1.65;

%% Q1) a) no current       %%
KC_nocurrentc = umax*T/D;
str1c='Value of Keulegan–Carpenter(KC) number for no current situation is '+string(KC_nocurrentc);
disp(str1c)
% Cd1c = input('Please provide Cd value and hit enter  ');
% Cm1c = input('Please provide Cm value and hit enter  ');
for i =1:length(z)
    u = (pi*H/T)*(cosh(k*(d+z(i)))/sinh(k*d)).*sind(theta);
    apre = (-2*pi^2*H/T^2)*cosh(k*(d+z(i))/sinh(k*d));
    a = apre.*cosd(theta);
    Fd1c(i,:) = 0.5*Cd1c*rho*D.*u.*abs(u);
    Fi1c(i,:) = (Cm1c*rho*pi*(D^2).*a)./4;
    Ft1c(i,:) = Fd1c(i,:) + Fi1c(i,:);
    Mt1c(i,:) = Ft1c(i,:)*abs(z(length(z)+1-i));
end

%% Q1 b)uniform current (in the direction of wave)%%
u_current_uniformc = 1;
KC_uniform_supportingc = (umax+u_current_uniformc)*T/D;
str2c='Value of Keulegan–Carpenter(KC) number for uniform current in wave direction is'+string(KC_uniform_supportingc);
disp(str2c);
% Cd2c = input('Please provide Cd value and hit enter  ');
% Cm2c = input('Please provide Cm value and hit enter  ');
for i =1:length(z)
    u = ((pi*H/T)*(cosh(k*(d+z(i)))/sinh(k*d)).*sind(theta));
    apre = (-2*pi^2*H/T^2)*cosh(k*(d+z(i))/sinh(k*d));
    a = apre.*cosd(theta);
    Fd2c(i,:) = 0.5*Cd2c*rho*D.*u.*abs(u)+0.5*Cd2c*rho*D.*u_current_uniformc.*abs(u_current_uniformc);
    Fi2c(i,:) = (Cm2c*rho*pi*(D^2).*a)./4;
    Ft2c(i,:) = Fd2c(i,:) + Fi2c(i,:);
    Mt2c(i,:) = Ft2c(i,:)*abs(z(length(z)+1-i));
end

%% Q1 c)uniform current (in the opposite direction of wave)%%
u_current_uniform_oppoc = -1;
KC_uniform_opposingc = (umax+u_current_uniform_oppoc)*T/D;
str3c='Value of Keulegan–Carpenter(KC) number for uniform current in the opposite direction of wave'+string(KC_uniform_opposingc);
disp(str3c);
% Cd3c = input('Please provide Cd value and hit enter  ');
% Cm3c = input('Please provide Cm value and hit enter  ');
for i =1:length(z)
    u = ((pi*H/T)*(cosh(k*(d+z(i)))/sinh(k*d)).*sind(theta));
    apre = (-2*pi^2*H/T^2)*cosh(k*(d+z(i))/sinh(k*d));
    a = apre.*cosd(theta);
    Fd3c(i,:) = 0.5*Cd3c*rho*D.*u.*abs(u)+0.5*Cd3c*rho*D.*u_current_uniform_oppoc.*abs(u_current_uniform_oppoc);
    Fi3c(i,:) = (Cm3c*rho*pi*(D^2).*a)./4;
    Ft3c(i,:) = Fd3c(i,:) + Fi3c(i,:);
    Mt3c(i,:) = Ft3c(i,:)*abs(z(length(z)+1-i));
end

%% Q1 d)varying current (in the supporting direction of wave)%%
u_current_varying_maxc = 1;
KC_varying_supportingc = (umax+u_current_varying_maxc)*T/D;
str3c='Value of Keulegan–Carpenter(KC) number for varying current in wave direction is'+string(KC_varying_supportingc);
disp(str3c);
% Cd4c = input('Please provide Cd value and hit enter  ');
% Cm4c = input('Please provide Cm value and hit enter  ');
for i =1:length(z)
    u = ((pi*H/T)*(cosh(k*(d+z(i)))/sinh(k*d)).*sind(theta));
    apre = (-2*pi^2*H/T^2)*cosh(k*(d+z(i))/sinh(k*d));
    a = apre.*cosd(theta);
    Fd4c(i,:) = 0.5*Cd4c*rho*D.*u.*abs(u) + 0.5*Cd4c*rho*D.*(1*((d+z(i))/d)^(1/7)).*abs((1*((d+z(i))/d)^(1/7)));
    Fi4c(i,:) = (Cm4c*rho*pi*(D^2).*a)./4;
    Ft4c(i,:) = Fd4c(i,:) + Fi4c(i,:);
    Mt4c(i,:) = Ft4c(i,:)*abs(z(length(z)+1-i));
end

%% Q1 e)varying current (in the opposing direction of wave)%%
u_current_varying_oppo_maxc = -1;
KC_varying_opposingc = (umax+u_current_varying_oppo_maxc)*T/D;
str4c='Value of Keulegan–Carpenter(KC) number for varying opposing current is'+string(KC_varying_opposingc);
disp(str4c);
% Cd5c = input('Please provide Cd value and hit enter  ');
% Cm5c = input('Please provide Cm value and hit enter  ');
for i =1:length(z)
    u = ((pi*H/T)*(cosh(k*(d+z(i)))/sinh(k*d)).*sind(theta));
    apre = (-2*pi^2*H/T^2)*cosh(k*(d+z(i))/sinh(k*d));
    a = apre.*cosd(theta);
    Fd5c(i,:) = 0.5*Cd5c*rho*D.*u.*abs(u) + 0.5*Cd4c*rho*D.*(-1*((d+z(i))/d)^(1/7)).*abs((-1*((d+z(i))/d)^(1/7)));
    Fi5c(i,:) = (Cm5c*rho*pi*(D^2).*a)./4;
    Ft5c(i,:) = Fd5c(i,:) + Fi5c(i,:);
    Mt5c(i,:) = Ft5c(i,:)*abs(z(length(z)+1-i));
end

%% Q1 plots no current%%
figure('Name','No_Current','NumberTitle','off','Units','normalized','Position',[0 0 1 1]);
subplot(2,3,1)
plot(theta,Ft1c(1,:),'r','linewidth',1.1);
hold on;
plot(theta, Fi1c(1,:),'b','linewidth',1.1);
hold on;
plot(theta, Fd1c(1,:),'g','linewidth',1.1);
title('z/d=0')
grid on;
xlim([0 360]);
ylim([min(Ft1c(1,:)) max(Ft1c(1,:))]);
xlabel('theta');
ylabel('F total');
legend('Total','Inertia','Drag');

subplot(2,3,2)
plot(theta,Ft1c(2,:),'r','linewidth',1.1);
hold on;
plot(theta, Fi1c(2,:),'b','linewidth',1.1);
hold on;
plot(theta, Fd1c(2,:),'g','linewidth',1.1);
title('z/d=-0.2')
grid on;
xlim([0 360]);
ylim([min(Ft1c(2,:)) max(Ft1c(2,:))]);
xlabel('theta');
ylabel('F total');
legend('Total','Inertia','Drag');

subplot(2,3,3)
plot(theta,Ft1c(3,:),'r','linewidth',1.1);
hold on;
plot(theta, Fi1c(3,:),'b','linewidth',1.1);
hold on;
plot(theta, Fd1c(3,:),'g','linewidth',1.1);
title('z/d=-0.4')
grid on;
xlim([0 360]);
ylim([min(Ft1c(3,:)) max(Ft1c(3,:))]);
xlabel('theta');
ylabel('F total');
legend('Total','Inertia','Drag','linewidth',1.1);

subplot(2,3,4)
plot(theta,Ft1c(4,:),'r','linewidth',1.1);
hold on;
plot(theta, Fi1c(4,:),'b','linewidth',1.1);
hold on;
plot(theta, Fd1c(4,:),'g','linewidth',1.1);
title('z/d=-0.6')
grid on;
xlim([0 360]);
ylim([min(Ft1c(4,:)) max(Ft1c(4,:))]);
xlabel('theta');
ylabel('F total');
legend('Total','Inertia','Drag');

subplot(2,3,5)
plot(theta,Ft1c(5,:),'r','linewidth',1.1);
hold on;
plot(theta, Fi1c(5,:),'b','linewidth',1.1);
hold on;
plot(theta, Fd1c(5,:),'g','linewidth',1.1);
title('z/d=-0.8')
grid on;
xlim([0 360]);
ylim([min(Ft1c(5,:)) max(Ft1c(5,:))]);
xlabel('theta');
ylabel('F total');
legend('Total','Inertia','Drag');

subplot(2,3,6)
plot(theta,Ft1c(6,:),'r','linewidth',1.1);
hold on;
plot(theta, Fi1c(6,:),'b','linewidth',1.1);
hold on;
plot(theta, Fd1c(6,:),'g','linewidth',1.1);
title('z/d=-1')
grid on;
xlim([0 360]);
ylim([min(Ft1c(6,:)) max(Ft1c(6,:))]);
xlabel('theta');
ylabel('F total');
legend('Total','Inertia','Drag');

%% Q1 plots all cases variation with depth %%
figure('Name','All_cases_with_depth_variation','NumberTitle','off','Units','normalized','Position',[0 0 1 1]);

subplot(2,3,1)
plot(theta,Ft1c(1,:),'r','linewidth',1.1);
hold on;
plot(theta, Ft2c(1,:),'b','linewidth',1.1);
hold on;
plot(theta, Ft3c(1,:),'g','linewidth',1.1);
hold on;
plot(theta, Ft4c(1,:),'k','linewidth',1.1);
hold on;
plot(theta, Ft5c(1,:),'m','linewidth',1.1);
title('z/d=0')
grid on;
xlim([0 360]);
xlabel('theta');
ylabel('F total');
legend('No current','Uniform(+)','Uniform(-)','Varying(+)','Varying(-)');

subplot(2,3,2)
plot(theta,Ft1c(2,:),'r','linewidth',1.1);
hold on;
plot(theta, Ft2c(2,:),'b','linewidth',1.1);
hold on;
plot(theta, Ft3c(2,:),'g','linewidth',1.1);
hold on;
plot(theta, Ft4c(2,:),'k','linewidth',1.1);
hold on;
plot(theta, Ft5c(2,:),'m','linewidth',1.1);
title('z/d=-0.2')
grid on;
xlim([0 360]);
xlabel('theta');
ylabel('F total');
legend('No current','Uniform(+)','Uniform(-)','Varying(+)','Varying(-)');

subplot(2,3,3)
plot(theta,Ft1c(3,:),'r','linewidth',1.1);
hold on;
plot(theta, Ft2c(3,:),'b','linewidth',1.1);
hold on;
plot(theta, Ft3c(3,:),'g','linewidth',1.1);
hold on;
plot(theta, Ft4c(3,:),'k','linewidth',1.1);
hold on;
plot(theta, Ft5c(3,:),'m','linewidth',1.1);
title('z/d=-0.4')
grid on;
xlim([0 360]);
xlabel('theta');
ylabel('F total');
legend('No current','Uniform(+)','Uniform(-)','Varying(+)','Varying(-)');

subplot(2,3,4)
plot(theta,Ft1c(4,:),'r','linewidth',1.1);
hold on;
plot(theta, Ft2c(4,:),'b','linewidth',1.1);
hold on;
plot(theta, Ft3c(4,:),'g','linewidth',1.1);
hold on;
plot(theta, Ft4c(4,:),'k','linewidth',1.1);
hold on;
plot(theta, Ft5c(4,:),'m','linewidth',1.1);
title('z/d=-0.6')
grid on;
xlim([0 360]);
xlabel('theta');
ylabel('F total');
legend('No current','Uniform(+)','Uniform(-)','Varying(+)','Varying(-)');

subplot(2,3,5)
plot(theta,Ft1c(5,:),'r','linewidth',1.1);
hold on;
plot(theta, Ft2c(5,:),'b','linewidth',1.1);
hold on;
plot(theta, Ft3c(5,:),'g','linewidth',1.1);
hold on;
plot(theta, Ft4c(5,:),'k','linewidth',1.1);
hold on;
plot(theta, Ft5c(5,:),'m','linewidth',1.1);
title('z/d=-0.8')
grid on;
xlim([0 360]);
xlabel('theta');
ylabel('F total');
legend('No current','Uniform(+)','Uniform(-)','Varying(+)','Varying(-)');

subplot(2,3,6)
plot(theta,Ft1c(6,:),'r','linewidth',1.1);
hold on;
plot(theta, Ft2c(6,:),'b','linewidth',1.1);
hold on;
plot(theta, Ft3c(6,:),'g','linewidth',1.1);
hold on;
plot(theta, Ft4c(6,:),'k','linewidth',1.1);
hold on;
plot(theta, Ft5c(6,:),'m','linewidth',1.1);
title('z/d=-1')
grid on;
xlim([0 360]);
xlabel('theta');
ylabel('F total');
legend('No current','Uniform(+)','Uniform(-)','Varying(+)','Varying(-)');


%% F Total at maximum(theta) %%
Ftotal_nocurrentc = 0.2*d*(sum(Ft1c,1) + Ft1c(1,:) + Ft1c(end,:))/2;
Ftotal_uniform_current_supportingc = 0.2*d*(sum(Ft2c,1) + Ft2c(1,:) + Ft2c(end,:))/2;
Ftotal_uniform_current_oppoc = 0.2*d*(sum(Ft3c,1) + Ft3c(1,:) + Ft3c(end,:))/2;
Ftotal_varying_current_supportingc = 0.2*d*(sum(Ft4c,1) + Ft4c(1,:) + Ft4c(end,:))/2;
Ftotal_varying_current_opposingc = 0.2*d*(sum(Ft5c,1) + Ft5c(1,:) + Ft5c(end,:))/2;

Mtotal_nocurrentc = 0.2*d*(sum(Mt1c,1) + Mt1c(1,:) + Mt1c(end,:))/2;
Mtotal_uniform_current_supportingc = 0.2*d*(sum(Mt2c,1) + Mt2c(1,:) + Mt2c(end,:))/2;
Mtotal_uniform_current_oppoc = 0.2*d*(sum(Mt3c,1) + Mt3c(1,:) + Mt3c(end,:))/2;
Mtotal_varying_current_supportingc = 0.2*d*(sum(Mt4c,1) + Mt4c(1,:) + Mt4c(end,:))/2;
Mtotal_varying_current_opposingc = 0.2*d*(sum(Mt5c,1) + Mt5c(1,:) + Mt5c(end,:))/2;

[~,I1c] = find(Ftotal_nocurrentc == max(Ftotal_nocurrentc));
[~,I2c] = find(Ftotal_uniform_current_supportingc == max(Ftotal_uniform_current_supportingc));
[~,I3c] = find(Ftotal_uniform_current_oppoc == max(Ftotal_uniform_current_oppoc));
[~,I4c] = find(Ftotal_varying_current_supportingc == max(Ftotal_varying_current_supportingc));
[~,I5c] = find(Ftotal_varying_current_opposingc == max(Ftotal_varying_current_opposingc));
figure('Name','All_cases_F_total','NumberTitle','off','Units','normalized','Position',[0 0 1 1]);
plot(theta,Ftotal_nocurrentc,'r','linewidth',1.1);
hold on;
plot(theta,Ftotal_uniform_current_supportingc,'b','linewidth',1.1);
hold on;
plot(theta,Ftotal_uniform_current_oppoc,'g','linewidth',1.1);
hold on;
plot(theta,Ftotal_varying_current_supportingc,'k','linewidth',1.1);
hold on;
plot(theta,Ftotal_varying_current_opposingc,'m','linewidth',1.1);
title(['F total']);
xlabel('theta');
ylabel('F total');
xlim([0 360])
grid on;
legend('No current','Uniform(+)','Uniform(-)','Varying(+)','Varying(-)');


%% Moment plot %%
figure('Name','All_cases_M_total','NumberTitle','off','Units','normalized','Position',[0 0 1 1]);
plot(theta,Mtotal_nocurrentc,'r','linewidth',1.1);
hold on;
plot(theta,Mtotal_uniform_current_supportingc,'b','linewidth',1.1);
hold on;
plot(theta,Mtotal_uniform_current_oppoc,'g','linewidth',1.1);
hold on;
plot(theta,Mtotal_varying_current_supportingc,'k','linewidth',1.1);
hold on;
plot(theta,Mtotal_varying_current_opposingc,'m','linewidth',1.1);
title(['M total']);
xlabel('theta');
ylabel('M total');
xlim([0 360])
grid on;
legend('No current','Uniform(+)','Uniform(-)','Varying(+)','Varying(-)');

%% Ftotal sectional vs z/d %%
figure('Name','Sectional force depth variation','NumberTitle','off','Units','normalized','Position',[0 0 1 1]);

plot(Ft1c(:,I1c),z,'r','linewidth',1.1);
hold on;
plot(Ft2c(:,I1c),z,'b','linewidth',1.1);
hold on;
plot(Ft3c(:,I1c),z,'g','linewidth',1.1);
hold on;
plot(Ft4c(:,I1c),z,'k','linewidth',1.1);
hold on;
plot(Ft5c(:,I1c),z,'m','linewidth',1.1);
title(['F sectional']);
xlabel('F total');
ylabel('depth');
grid on;
legend('No current','Uniform(+)','Uniform(-)','Varying(+)','Varying(-)');

%% Mtotal sectional vs z/d %%
figure('Name','Sectional moment depth variation','NumberTitle','off','Units','normalized','Position',[0 0 1 1]);

plot(Mt1c(:,I1c),z,'r','linewidth',1.1);
hold on;
plot(Mt2c(:,I1c),z,'b','linewidth',1.1);
hold on;
plot(Mt3c(:,I1c),z,'g','linewidth',1.1);
hold on;
plot(Mt4c(:,I1c),z,'k','linewidth',1.1);
hold on;
plot(Mt5c(:,I1c),z,'m','linewidth',1.1);
title(['M sectional']);
xlabel('M total');
ylabel('depth');
grid on;
legend('No current','Uniform(+)','Uniform(-)','Varying(+)','Varying(-)');

%% Moment calculation
Ftmax_nocurrentc = Ftotal_nocurrentc(I1c);
Ftmax_uniform_current_supportingc = Ftotal_uniform_current_supportingc(I1c);
Ftmax_unform_oppoc = Ftotal_uniform_current_oppoc(I1c);
Ftmax_varying_current_supportingc = Ftotal_varying_current_supportingc(I1c);
Ftmax_varying_current_opposingc = Ftotal_varying_current_opposingc(I1c);
disp(['F total max_nocurrent (N) =                  ',num2str(Ftmax_nocurrentc)]);
disp(['F total max_uniform_current_supporting (N) = ',num2str(Ftmax_uniform_current_supportingc)]);
disp(['F total max_unform_oppo (N) =                ',num2str(Ftmax_unform_oppoc)]);
disp(['F total max_varying_current_supporting (N) = ',num2str(Ftmax_varying_current_supportingc)]);
disp(['F total max_varying_current_opposing (N) =   ',num2str(Ftmax_varying_current_opposingc)]);
disp('')
disp('')

Mtmax_nocurrentc = (2*Ft1c(1,I1c)*d+Ft1c(2,I1c)*0.8*d+Ft1c(3,I1c)*0.6*d+Ft1c(4,I1c)*0.4*d+Ft1c(5,I1c)*0.2*d)*0.2*d/2;
Mtmax_uniform_current_supportingc = (2*Ft2c(1,I1c)*d+Ft2c(2,I1c)*0.8*d+Ft2c(3,I1c)*0.6*d+Ft2c(4,I1c)*0.4*d+Ft2c(5,I1c)*0.2*d)*0.2*d/2;
Mtmax_unform_oppoc = (2*Ft3c(1,I1c)*d+Ft3c(2,I1c)*0.8*d+Ft3c(3,I1c)*0.6*d+Ft3c(4,I1c)*0.4*d+Ft3c(5,I1c)*0.2*d)*0.2*d/2;
Mtmax_varying_current_supportingc = (2*Ft4c(1,I1c)*d+Ft4c(2,I1c)*0.8*d+Ft4c(3,I1c)*0.6*d+Ft4c(4,I1c)*0.4*d+Ft4c(5,I1c)*0.2*d)*0.2*d/2;
Mtmax_varying_current_opposingc = (2*Ft5c(1,I1c)*d+Ft5c(2,I1c)*0.8*d+Ft5c(3,I1c)*0.6*d+Ft5c(4,I1c)*0.4*d+Ft5c(5,I1c)*0.2*d)*0.2*d/2;

disp(['M total max_nocurrent (N-m) =                  ',num2str(Mtmax_nocurrentc)]);
disp(['M total max_uniform_current_supporting (N-m) = ',num2str(Mtmax_uniform_current_supportingc)]);
disp(['M total max_unform_oppo (N-m) =                ',num2str(Mtmax_unform_oppoc)]);
disp(['M total max_varying_current_supporting (N-m) = ',num2str(Mtmax_varying_current_supportingc)]);
disp(['M total max_varying_current_opposing (N-m) =   ',num2str(Mtmax_varying_current_opposingc)]);
disp('')
disp('')


h1c = Mtmax_nocurrentc/Ftmax_nocurrentc;
h2c = Mtmax_uniform_current_supportingc/Ftmax_uniform_current_supportingc;
h3c = Mtmax_unform_oppoc/Ftmax_unform_oppoc;
h4c = Mtmax_varying_current_supportingc/Ftmax_varying_current_supportingc;
h5c = Mtmax_varying_current_opposingc/Ftmax_varying_current_opposingc;


disp(['lever arm nocurrent (m) =             ',num2str(h1c) ]);
disp(['lever arm uniform current (+ve) (m) = ',num2str(h2c) ]);
disp(['lever arm uniform current (-ve) (m) = ',num2str(h3c) ]);
disp(['lever arm varying current (+ve) (m) = ',num2str(h4c) ]);
disp(['lever arm varying current (-ve) (m) = ',num2str(h5c) ]);


figure
names = {'noCurr'; 'Uni(+)'; 'Uni(-)'; 'Vary(+)'; 'Vary(-)'};
X = [1,2,3,4,5];
Y = [h1c,h2c,h3c,h4c,h5c];
stem(X,Y,'LineStyle','-.',...
     'MarkerFaceColor','red',...
     'MarkerEdgeColor','green')
set(gca,'xtick',[1:5],'xticklabel',names)
xlim([0 6])
ylabel('Lever arm(m)')

figure
names = {'noCurr'; 'Uni(+)'; 'Uni(-)'; 'Vary(+)'; 'Vary(-)'};
X = [1,2,3,4,5];
Y = [Ftmax_nocurrentc, Ftmax_uniform_current_supportingc,Ftmax_unform_oppoc,Ftmax_varying_current_supportingc,Ftmax_varying_current_opposingc];
stem(X,Y,'LineStyle','-.',...
     'MarkerFaceColor','red',...
     'MarkerEdgeColor','green')
set(gca,'xtick',[1:5],'xticklabel',names)
xlim([0 6])
ylabel('Fmax(N)')

figure
names = {'noCurr'; 'Uni(+)'; 'Uni(-)'; 'Vary(+)'; 'Vary(-)'};
X = [1,2,3,4,5];
Y = [Mtmax_nocurrentc, Mtmax_uniform_current_supportingc,Mtmax_unform_oppoc,Mtmax_varying_current_supportingc,Mtmax_varying_current_opposingc];
stem(X,Y,'LineStyle','-.',...
     'MarkerFaceColor','red',...
     'MarkerEdgeColor','green')
set(gca,'xtick',[1:5],'xticklabel',names)
xlim([0 6])
ylabel('Mmax(N-m)')

