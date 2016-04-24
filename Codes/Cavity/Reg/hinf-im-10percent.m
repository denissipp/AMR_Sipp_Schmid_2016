clear all
clf reset
format long;

link = './';

ngm=4;
fid=fopen('../Common/param.dat','rt');
num=fscanf(fid,'%d',1)
dt=fscanf(fid,'%f',1)
stp=fscanf(fid,'%d',1)
p=fscanf(fid,'%d',1)
nstab=fscanf(fid,'%d',1)
fclose(fid);
nt=2*ngm+p; % total size of stored state space model
ns=2*ngm+nstab; % actual size of state-space model

%% Read reduced state-space model
file=fopen('../ROM/A.txt','r+');
dat=fscanf(file,'%g',[nt,nt]);
fclose(file);
A=dat(1:ns,1:ns);

file=fopen('../ROM/B2.txt','r+');
dat=fscanf(file,'%g',[nt,1]);
fclose(file);
B=dat(1:ns,1);

file=fopen('../ROM/C.txt','r+');
dat=fscanf(file,'%g',[1,nt]);
fclose(file);
C=dat(1,1:ns);
%%
% Input-multiplicative perturbations (matrices Bd and Cd)
G=1.e-2;	% is G/W
%%
D=[0 0 1; 0 sqrt(G) 0];
sys1=ss(A,[B zeros(ns,1) B],[zeros(1,ns); C],D,'InputName',{'wd','g','u'},'OutputName',{'zd','y'});

%% Initialize algorithm with H2 solution
% Solve Riccati equation
H=[A -B*B'; zeros(ns,ns) -A'];
[U,H]=schur(H,'complex');
[Us,Hs]=ordschur(U,H,'lhp');
n=size(A,1);
Xc= Us(n+1:2*n,1:n)/Us(1:n,1:n); Xc=real(Xc);
% Solve Riccati equation
H=[A' -C'*C/G; -B*B' -A];
[U,H]=schur(H,'complex');
[Us,Hs]=ordschur(U,H,'lhp');
Xe= Us(n+1:2*n,1:n)/Us(1:n,1:n); Xe=real(Xe);
% Controller in state-space form
Cc=-B'*Xc;
Bc=Xe*C'/G;
Ac=A+B*Cc-Bc*C;
sys2=ss(Ac,Bc,Cc,0,'InputName','y','OutputName','u');
% Evaluation of closed-loop infinity norm
H=feedback(sys1,sys2,3,2,1);
normold=norm(H(1,1:2),inf)
%% Loop on gamma to obtain Hinf solution
gammaold=1.e30;
Acold=Ac;
Bcold=Bc;
Ccold=Cc;
jjinit=log10(norm(H(1,1:2),inf)*2);
for jj=jjinit:-0.01:-100
    gamma=10^jj
    % solve Riccati equation
    H=[A -B*B'*(1.-1./gamma^2); zeros(ns,ns) -A'];
    [U,H]=schur(H,'complex');
    [Us,Hs]=ordschur(U,H,'lhp');
    n=size(A,1);
    Xc= Us(n+1:2*n,1:n)/Us(1:n,1:n); Xc=real(Xc);
    % solve Riccati equation
    H=[A' -C'*C/G; -B*B' -A];
    [U,H]=schur(H,'complex');
    [Us,Hs]=ordschur(U,H,'lhp');
    Xe= Us(n+1:2*n,1:n)/Us(1:n,1:n); Xe=real(Xe);
    % exit criterion 
    if min(real(eig(eye(ns)-Xe*Xc/gamma^2))) <= 0
        break
    end
    % update controller 
    Cc=-B'*Xc;  
	Bc=(eye(ns)-Xe*Xc/gamma^2)\(Xe*C'/G);
	Ac=A+(1.-1./gamma^2)*B*Cc-Bc*C;
    sys2=ss(Ac,Bc,Cc,0,'InputName','y','OutputName','u');
    % update controller and norm
    H=feedback(sys1,sys2,3,2,1);
    normnew=norm(H(1,1:2),inf);
    if(normnew>gamma)
        break;
    end
    normold=normnew
    gammaold=gamma;
    Acold=Ac;
    Bcold=Bc;
    Ccold=Cc;
end
%% Final controller 10 percent above minimum
gamma=gammaold*1.1
H=[A -B*B'*(1.-1./gamma^2); zeros(ns,ns) -A'];
[U,H]=schur(H,'complex');
[Us,Hs]=ordschur(U,H,'lhp');
n=size(A,1);
Xc= Us(n+1:2*n,1:n)/Us(1:n,1:n); Xc=real(Xc);
% solve Riccati equation
H=[A' -C'*C/G; -B*B' -A];
[U,H]=schur(H,'complex');
[Us,Hs]=ordschur(U,H,'lhp');
Xe= Us(n+1:2*n,1:n)/Us(1:n,1:n); Xe=real(Xe);
% update controller 
Ccold=-B'*Xc;  
Bcold=(eye(ns)-Xe*Xc/gamma^2)\(Xe*C'/G);
Acold=A+(1.-1./gamma^2)*B*Ccold-Bcold*C;
%% Write controller in state-space form to file
AAfile=fopen('J.txt','w');
for i=1:size(Acold,1)
  for jj=1:size(Acold,1)
   fprintf(AAfile,'%.15g\n',Acold(jj,i));
  end
end
fclose(AAfile);
BBfile=fopen('L.txt','w');
CCfile=fopen('K.txt','w');
for i=1:size(Acold,1)
 fprintf(BBfile,'%.15g\n',Bcold(i));
 fprintf(CCfile,'%.15g\n',Ccold(i));
end
fclose(BBfile);
fclose(CCfile);
DDfile=fopen('M.txt','w');
fprintf(DDfile,'%.15g\n',0.);
fclose(DDfile);
%% Performance and robusteness evaluations
B1=[zeros(ns,1) B]; Bd=B; B2=B; C1=[C; zeros(1,ns)]; Cd=[zeros(1,ns)]; C2=C; D=[0 0 0 0; 0 0 0 1; 0 1 0 1; 1 0 0 0];
sys1=ss(A,[B1 Bd B2],[C1; Cd; C2],D,'InputName',{'g','w','wdu','u'},'OutputName',{'zy','zu','zd','y'});
sys2=ss(Acold,Bcold,Ccold,0,'InputName','y','OutputName','u');
H=feedback(sys1,sys2,4,4,1);
norm(H(3,3),inf)
gammaold

% GM+
sys=ss(A,B2,C2,0);
gmpmaxi=10;
f=@(x) max(real(zero(1-x*sys*sys2)));
while f(1)*f(gmpmaxi)>=0 && gmpmaxi < 10e25
   gmpmaxi=gmpmaxi*10;
end
if gmpmaxi > 10e20
        break
else
   gmp = fzero(f,[1 gmpmaxi],optimset('TolX',1.e-6))
end

% GM-
if f(0)*f(1)>=0
   gmm=0
else
   gmm=fzero(f,[0 1])
end

% PM    
g=@(x) max(real(zero(1-exp(-j*x)*sys*sys2)));
if g(0)*g(pi/2)>=0
   pm=pi/2
else
   pm=fzero(g,[0 pi/2])
end

%% Write to file
file=fopen(['perfo-im-10percent.txt'],'wt');
fprintf(file,'%s\n','VARIABLES= "N2zg",  "N2zw", "N2ug", "N2uw", "GM+", "GM-", "PM", "rho"');
fprintf(file,'%13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n',...
               norm(H(1,1),2),norm(H(1,2),2),norm(H(2,1),2),norm(H(2,2),2),20*log10(gmp),20*log10(gmm),pm*180/pi,1./norm(H(3,3),inf));
fclose(file);
