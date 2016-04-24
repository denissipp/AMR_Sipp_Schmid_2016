fid=fopen('../Common/param.dat','rt');
N=fscanf(fid,'%d',1)	% number of iterations in direct and adjoint simulations
dt=fscanf(fid,'%f',1)	% time-step
fclose(fid);

fny=1./2/dt;       % Nyquist frequency 1/(2*dt)
fcut=5./2./pi;    % Cut-off frequency

finput=fcut/fny;    % Normalized cut-off frequency

u=0.0054*idinput(N,'rbs',[0 finput],[-1 1]);   % generate random-binary-signal

file=fopen('uexplor.txt','wt');	% write signal to file
fprintf(file,'%i\n',N);
for i=1:N
    fprintf(file,'%e\n',u(i));
end

