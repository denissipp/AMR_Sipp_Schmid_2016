clear all;
close all;
format long;

%% Constants
tol      = 0.75;	% Threshold for pseudo-inverse in control design (agressivity of controller)
finaldim = 14;      	% Final dimension of controller
resamp  = 50; 		% Number of time-steps defining sampling time

%%
fid=fopen('../Common/param.dat','rt');
L1=fscanf(fid,'%d',1)	% number of iterations in direct and adjoint simulations
dt_simu=fscanf(fid,'%f',1)	% time-step
fclose(fid);
Ts = dt_simu*resamp;		% sampling time

%% Load input-output data 
fid=fopen('../DNSLEARN/bfs.txt');	
str=char('%g %g %g %g %g %g\n');
D=fscanf(fid,str,[6 L1]);
for(i=1:L1)
    time(i)=rd(D(1,i),4);   % only consider 4 significant digits 
    yr(i)=rd(D(3,i),4);
    ur(i)=rd(D(6,i),4);
    zr(i)=rd(D(4,i),4);
end
datar = iddata(zr',[ur' yr'],dt_simu);
datar.InputName  = {'u'; 'y'};
datar.OutputName = 'z';
datar.TimeUnit   = 'Seconds';

%% Definition of filter to avoid aliasing when resampling
omc=2;			% cut-off frequency of filter
F=tf(omc^2,[1 2*0.8*omc omc^2]);
Fd=c2d(F,dt_simu);
% Filter data
uf = filter(Fd.num{1},Fd.den{1},ur);
yf = filter(Fd.num{1},Fd.den{1},yr);
zf = filter(Fd.num{1},Fd.den{1},zr);
dataf = iddata(zf',[uf' yf'],dt_simu);
dataf.InputName  = {'u'; 'y'};
dataf.OutputName = 'z';
dataf.TimeUnit   = 'Seconds';
% Represent unfiltered (red) and filtered (blue) unsampled data
figure(1);
subplot(3,1,1)
plot(time,yf);
hold on;
plot(time,yr,'r');
axis([200,250,-0.001,0.001]);
ylabel('y');
subplot(3,1,2)
plot(time,uf)
hold on;
plot(time,ur,'r');
axis([200,250,-0.0075,0.0075]);
ylabel('u');
subplot(3,1,3)
plot(time,zf)
hold on;
plot(time,zr,'r');
axis([200,250,-0.0075,0.0075]);
xlabel('t');
ylabel('z');

%% Resample data
datafs = resample(dataf,1,resamp,0); % No filter, just take one point out of resamp
datafs.InputName = dataf.InputName;
datafs.OutputName = dataf.OutputName;
% Define learning and validation datasets
datafs_learn  = datafs(500:3000);   % remove transient at beginning of simulation
datafs_valid  = datafs(3001:end);

%% Armax regression
na=5/Ts
nb=[6/Ts 6/Ts]
%na=0/Ts;
%nb=[30/Ts 30/Ts];
nk=[10/Ts 10/Ts];
nc=10;
model = armax(datafs_learn,[na nb nc nk],'Init','Zero');
Estim = ss(model);

[ylearn tlearn] = lsim(Estim, datafs_learn.u);	% simulate model with input stemming from learning dataset 
learning=norm(ylearn(300:end)-datafs_learn.y(300:end))/norm(datafs_learn.y(300:end))	% evaluate error; we remove the initial transient
[yvalid tvalid] = lsim(Estim, datafs_valid.u);	% simulate model with input stemming from validation dataset
valid=norm(yvalid(300:end)-datafs_valid.y(300:end))/norm(datafs_valid.y(300:end))	% evaluate error;

figure(2)	
subplot(2,1,1)  % represent impulse in u
impulse(Estim(1,1))
ylabel('u');
xlabel('t');
subplot(2,1,2)	% represent impulse in y
impulse(Estim(1,2))
ylabel('z');
xlabel('t');

figure(3);
subplot(2,1,1)
plot(tlearn,ylearn,tlearn,datafs_learn.y);
title('learning dataset')
ylabel('z');
subplot(2,1,2)
plot(tvalid,yvalid,tvalid,datafs_valid.y);  
title('validation dataset')
xlabel('t');
ylabel('z');

%% Feed-forward controller
Tmax = 100; % Time-horizon for feedforward controller
Nmax = floor(Tmax/Estim.Ts);

% Matrices Hu, Gu, Gs
Imp_u = impulse(Estim(:,1),Tmax); % Parametres de Markov_u
Imp_s = impulse(Estim(:,2),Tmax); % Parametres de Markov_s

%--------------- Hu ---------------
clear Hu; Hu = zeros(Nmax+1,Nmax+1);
for i=1:Nmax+1
  for j=1:i
	Hu(i,j) = Imp_u(i-j+1);
  end
end

%--------------- Gu ---------------
clear Gu; Gu = zeros(Nmax+1,Nmax+1);
for i=1:Nmax
  for j=1:Nmax-i+1
	Gu(i,j) = Imp_u(i+j);
  end
end

%--------------- Gs ---------------
clear Gs; Gs = zeros(Nmax+1,Nmax+1);
for i=1:Nmax+1
  for j=1:Nmax+1-i+1
	Gs(i,j) = Imp_s(i+j-1);
  end
end

% Pseudo-inverse 
Huplus = pinv(Hu,tol);  % tol is the tolerance. Modulates aggressivity of controller
RegBu = -Huplus*Gu;
RegBs = -Huplus*Gs;

% Impulse response of controller 
Ub = zeros(Nmax+1,1);
Sb = zeros(Nmax+1,1); Sb(1) = 1./Ts;
proj = zeros(1,Nmax+1); proj(1) = 1; 

u = zeros(Tmax/Estim.Ts,1);
for i=1:Nmax;
  u(i) = proj*(RegBu*Ub + RegBs*Sb);
  for k=Nmax+1:-1:2
	Ub(k) = Ub(k-1);
	Sb(k) = Sb(k-1);
  end
  Ub(1) = u(i);
  Sb(1) = 0; 
end

%% Represent controller as state-space model 

[a,b,c,d,totbnd,hsv]=imp2ss(u*Ts,-Ts,1,1,1.e-10);
% controller with dimension finaldim
[a,b,c,d,totbnd,hsv]=imp2ss(u*Ts,-Ts,1,1,2*norm(hsv(finaldim+1:end),1));
bregulator=ss(a,b,c,d,Ts);
bregulator.InputName = 'y';
bregulator.OutputName = 'u';

figure(4);  % impulse response of controller
plot(0:Ts:(Tmax/Estim.Ts-1)*Ts,u);
hold on
impulse(bregulator,'g') % in state-space form with 14 degrees of freedom
xlabel('t');
ylabel('u');
%% Closed-loop system
abclsys  = connect(Estim, bregulator,{'y'}, {'z','u'}); 

figure(5);  % Closed-loop system in time-domain (y->z)
impulse(abclsys(1,1),Estim(1,2));
figure(6);  % Closed-loop system in frequency domain (y->z)
bode(abclsys(1,1),Estim(1,2));
%% Run a simulation based on sampling time 
test = datafs_valid(1:end);

% uncontrolled simulation :
[zunc x ]  = lsim(Estim, [0*test.u(:,1) test.u(:,2) ]);
% controlled simulation
[zc x ]   = lsim(abclsys, [test.u(:,2)]);
% efficiency
reduction = (norm(zunc) - norm(zc(:,1)))/norm(zunc)

figure(7);
subplot(2,1,1)
plot(x,zunc,'r-',x,zc(:,1),'linewidth',1);legend('Uncontrolled','Controlled');
subplot(2,1,2)
plot(x,zc(:,2),'linewidth',1); legend('Control law');
%% Back to the DNS timestep
f_Estim            = d2d(Estim     , dt_simu  , 'tustin');
f_bregulator       = d2d(bregulator, dt_simu  , 'tustin');
f_abclsys          = connect(f_Estim,f_bregulator,{'y'}, {'z','u'});

figure(8)   % Bode diagram of regulator based on sampling time Ts and DNS time-step
bode(bregulator,f_bregulator)
xlabel('Frequency');
%% Run a simulation based on DNS time-step
clear test;
test = datar(1:end);

% uncontrolled simulation :
[zunc x ]  = lsim(f_Estim, [0*test.u(:,1) test.u(:,2) ]);
% controlled simulation
[zc x ]    = lsim(f_abclsys, [test.u(:,2)]);
% efficiency
reduction = (norm(zunc) - norm(zc(:,1)))/norm(zunc)

figure(9);
subplot(2,1,1)
plot(x,zunc,'r-',x,zc(:,1),'linewidth',1); legend('Uncontrolled','Controlled');
subplot(2,1,2)
plot(x,zc(:,2),'linewidth',1) ; legend('Control law');
%% Export regulator
bout = fopen('regulator.txt','w');
fprintf(bout,'%i \n',size(f_bregulator.a));
fprintf(bout,'%e \n', f_bregulator.a);
fprintf(bout,'%i \n',size(f_bregulator.b));
fprintf(bout,'%e \n', f_bregulator.b);
fprintf(bout,'%i \n',size(f_bregulator.c));
fprintf(bout,'%e \n', f_bregulator.c);
fprintf(bout,'%i \n',size(f_bregulator.d));
fprintf(bout,'%e \n', f_bregulator.d);
fclose(bout);
