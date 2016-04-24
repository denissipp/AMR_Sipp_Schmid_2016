clear all
clf reset
format long;

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
fid=fopen('controlparam.dat','rt');
G=fscanf(fid,'%f',1)
l2=fscanf(fid,'%f',1)
fclose(fid);
%%
B2=B; C1=C;
% solve Riccati equation
H=[A -B2*B2'/l2; -C1'*C1 -A'];
[U,H]=schur(H,'complex');
[Us,Hs]=ordschur(U,H,'lhp');
X= Us(ns+1:2*ns,1:ns)/Us(1:ns,1:ns); X=real(X);
K=-(1./l2)*B2'*X;
% solve Riccati equation
H=[A' -C1'*C1/G; -B2*B2' -A];
[U,H]=schur(H,'complex');
[Us,Hs]=ordschur(U,H,'lhp');
X= Us(ns+1:2*ns,1:ns)/Us(1:ns,1:ns); X=real(X);
L=(1./G)*X*C1';
%% Write controller in state-space form to file
J=A-L*C1+B2*K;
AAfile=fopen('J.txt','w');
for i=1:size(J,1)
  for jj=1:size(J,1)
    fprintf(AAfile,'%.15g\n',J(jj,i));
  end
end
fclose(AAfile);
BBfile=fopen('L.txt','w');
CCfile=fopen('K.txt','w');
for i=1:size(J,1)
  fprintf(BBfile,'%.15g\n',L(i));
  fprintf(CCfile,'%.15g\n',K(i));
end
fclose(BBfile);
fclose(CCfile);
DDfile=fopen('M.txt','w');
fprintf(DDfile,'%.15g\n',0.);
fclose(DDfile);
%% Performance and robusteness evaluations
B1=[zeros(ns,1) B]; Bd=B; B2=B; C1=[C; zeros(1,ns)]; Cd=[zeros(1,ns)]; C2=C; D=[0 0 0 0; 0 0 0 1; 0 1 0 1; 1 0 0 0];
sys1=ss(A,[B1 Bd B2],[C1; Cd; C2],D,'InputName',{'g','w','wdu','u'},'OutputName',{'zy','zu','zd','y'});
sys2=ss(J,L,K,0,'InputName','y','OutputName','u');
H=feedback(sys1,sys2,4,4,1);

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
file=fopen(['perfo-lqg.txt'],'wt');
fprintf(file,'%s\n','VARIABLES= "N2zg",  "N2zw", "N2ug", "N2uw", "GM+", "GM-", "PM", "rho"');
fprintf(file,'%13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n',...
               norm(H(1,1),2),norm(H(1,2),2),norm(H(2,1),2),norm(H(2,2),2),20*log10(gmp),20*log10(gmm),pm*180/pi,1./norm(H(3,3),inf));
fclose(file);

