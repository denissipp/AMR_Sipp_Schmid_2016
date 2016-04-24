fid=fopen('../Common/param.dat','rt');
num=fscanf(fid,'%d',1)	% number of iterations in direct and adjoint simulations
dt=fscanf(fid,'%f',1)	% time-step
step=fscanf(fid,'%d',1)	% number of time-step between two snapshots
p=fscanf(fid,'%d',1)	% total number of balanced modes
nbre=(num-1)/step	% number of snapshots
fclose(fid);

% determine size of snapshots
fid=fopen(strcat('sol_dir/cbf_100001_.txt'),'rt');
n=fscanf(fid,'%d',1);
fclose(fid);

% read gramian
fid=fopen('gramian.txt','rt');
C=textscan(fid,'%f');
A=reshape(C{1},nbre+1,nbre+1);
clear C;
fclose(fid);

% perform singular value decomposition
[U,S,V]=svd(A);

% compute balanced modes
J=zeros(n,p);
for i=0:nbre	% loop on snapshots of direct simulation
    i
    if ((i==0)|(i==nbre)) 
            mult=sqrt(1.*step/3.*dt);
    elseif(mod(i,2)==1)
            mult=sqrt(4.*step/3.*dt);
	else
            mult=sqrt(2.*step/3.*dt);
    end

    ii=i*step;
    fid=fopen(strcat('sol_dir/cbf_',num2str(100000+ii+1),'_.txt'),'rt');
    num=fscanf(fid,'%d',1);
    b=textscan(fid,'%f',num);
    J=J+mult*b{1}*V(i+1,1:p)*S(1:p,1:p)^(-1/2);	% update balanced mode
    fclose(fid);
end

for j=1:p	% save direct balanced modes
    j
    fid=fopen(strcat('BPOD/mode_',num2str(j),'_.txt'),'wt');
    fprintf(fid,'%d\n',n);
    fprintf(fid,'%21.14e\n',J(:,j));
    fclose(fid);
end

% compute adjoint balanced modes
J=zeros(n,p);
for i=0:nbre	% loop on snapshots of adjoint simulation
    i
    if ((i==0)|(i==nbre)) 
            mult=sqrt(1.*step/3.*dt);
    elseif(mod(i,2)==1)
            mult=sqrt(4.*step/3.*dt);
	else
            mult=sqrt(2.*step/3.*dt);
    end

    ii=i*step;
    fid=fopen(strcat('sol_adj/adj_',num2str(100000+ii+1),'_.txt'),'rt');
    num=fscanf(fid,'%d',1);
    b=textscan(fid,'%f',num);
    J=J+mult*b{1}*U(i+1,1:p)*S(1:p,1:p)^(-1/2);	% update adjoint balanced mode
    fclose(fid);
end
    
for j=1:p	% save adjoint balanced modes
    j
    fid=fopen(strcat('BPOD/modeadj_',num2str(j),'_.txt'),'wt');
    fprintf(fid,'%d\n',n);
    fprintf(fid,'%21.14e\n',J(:,j));
    fclose(fid);
end
