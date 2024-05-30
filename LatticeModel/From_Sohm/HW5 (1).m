%%%% Program EPM (based on the Cohen-Bergstresser 1966 paper)
%%%%
%%%% Read in G vectors included in the calculation
fid=fopen('g.txt','r');
G=fscanf(fid,'%d',[3,65]);
G';

%%%% read in K-pts
num_k = 100;
fid=fopen('k_new.txt','r');
K=fscanf(fid,'%f',[num_k,3]);
K=transpose(K);
K';
status=fclose('all');
fid=fopen('epm_out.txt','w');
fprintf(fid,' EPM OUPUT\n');

%%% set parameters
% form factors
% lattice constant
v3=-0.115 ;
v8=0.005 ;
v11=0.030 ;
 
% v3=0.0 ;
% v8=0.0 ;
% v11=0.0 ;

a= 5.43 /0.529;
con=(2*pi/a)^2;
energies = [];

% cut offs (k+G)^2 and potential
qcut=11.1;
gcut=11.1;


%%% Loop over k-vectors
for nk=1:num_k
xk=K(1,nk);
yk=K(2,nk);
zk=K(3,nk);

%%% Sort out the (k+G) space
n=0;
for j=1:65
test=(xk+G(1,j))^2+(yk+G(2,j))^2+(zk+G(3,j))^2;
if test < gcut
n=n+1;
GK(1,n)=xk+G(1,j);
GK(2,n)=yk+G(2,j);
GK(3,n)=zk+G(3,j);
H(n,n)=0.5*con*(GK(1,n)^2+GK(2,n)^2+GK(3,n)^2);
end
end
H';
n;

%%% Construct H matrix for off diagonal
%%
for j1=1:n
m=j1+1;
for j2=m:n
gx=GK(1,j1)-GK(1,j2);
gy=GK(2,j1)-GK(2,j2);
gz=GK(3,j1)-GK(3,j2);
sg=cos(0.25*pi*(gx+gy+gz));
gtest=gx^2+gy^2+gz^2;
v=0.;
if abs(gtest-3) < 0.01 v=v3; end
if abs(gtest-8) < 0.01 v=v8; end
if abs(gtest-11) <0.01 v=v11; end
H(j1,j2)=sg*v;
H(j2,j1)=H(j1,j2);
end
end

%% Find eigenvalues
E=eig(H);
E=27.2.*E-9.5623;   %What are these numbers?
energies = cat(2,energies,E(1:51));
fprintf(fid,'\n \n k-point %5.1f %5.1f %5.1f', xk,yk,zk);
for jj=1:8
fprintf(fid,'\n %8.3f',E(jj));
end
%clear H
end
status=fclose('all');


% Band structure:
x = length(K);
plot(transpose(energies),'k','LineWidth',2)
title("Problem 3: EPM", 'FontSize', 12)
ylabel("Energy", 'FontSize', 12)
xticks([0, 50, 100])
xticklabels({'L', 'Gamma', 'X'})
ylim([-14,31])


