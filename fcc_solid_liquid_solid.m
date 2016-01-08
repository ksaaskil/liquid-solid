% Create a solid-liquid-solid system
clear
filetowrite='121015g_sls.dat';
% mode='bulk';
mode='interface';

pbc=1;

ML=1;
MR=1;


alatt=1.56; % Lattice constant in units of sigma for solid argon at 0 K
%alatt=1.56;
% a=1; % Lattice constant in units of sigma for solid argon at 40 K
%alatt=1.579; % Multiplier for final multiplication of coordinates and box sizes
sigma=3.4; % Solid argon
%a=1; 

b1=[0,0,0];
b2=alatt*0.5*[0,1,1];
b3=alatt*0.5*[1,0,1];
b4=alatt*0.5*[1,1,0];

%b1=0; % The single atom

% The basis vectors of the cubic unit cell
xlatt=alatt;
ylatt=alatt;
zlatt=alatt;

% How many unit cells (each containing 8 atoms since we employ an
% orthogonal basis) are included for each bulk
%Nx=40;
Nx=30;
Ny=20;%Ny=1;
Nz=Ny;%Nz=1;

Ly=Ny*ylatt;
Lz=Nz*zlatt;

% Lengths of the bulk parts
L_bulk=Nx*xlatt;
L_spare=alatt/2-alatt/4-alatt/4;

% Length of the liquid part
% L_liquid=Nx*alatt;
L_liquid=0*alatt;
%L_liquid=0;
% Volume of liquid atoms
V_liquid=L_liquid*Ly*Lz;
% Liquid density
rho=0.7; % http://dx.doi.org/10.1016/j.ijheatmasstransfer.2006.03.005
% Number of particles
N_liquid=round(rho*V_liquid);
% N_liquid=0;

% Length of the total system
Lx=L_bulk+L_liquid+L_bulk+2*L_spare;
% Number of atoms in each bulk
N_bulk=4*Nx*Ny*Nz;

N=2*N_bulk+N_liquid;

% Atom positions
catom=zeros(N,3);
type=zeros(N,1);

%kmask=zeros(4*Nx*Ny*Nz,2);

lask=0;
lask2=1;

kmask=zeros(N,2);

for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            lask=lask+1;
            %clatt(lask,:)=i*a1+j*a2+k*a3;
            clatt=(i-1)*[xlatt,0,0]+(j-1)*[0,ylatt,0]+(k-1)*[0,0,zlatt];
            % First bulk
            catom(lask2,:)=clatt+b1;
            catom(lask2+1,:)=clatt+b2;
            catom(lask2+2,:)=clatt+b3;
            catom(lask2+3,:)=clatt+b4;
            type(lask2:lask2+3)=1;
            kmask(lask2:lask2+3,:)=repmat([j,k],4,1);
            % Second bulk
            Lshift=L_bulk+L_liquid+2*L_spare;
            catom(lask2+N_bulk+N_liquid,:)=[clatt(1)+Lshift,clatt(2),clatt(3)]+b1;
            catom(lask2+1+N_bulk+N_liquid,:)=[clatt(1)+Lshift,clatt(2),clatt(3)]+b2;
            catom(lask2+2+N_bulk+N_liquid,:)=[clatt(1)+Lshift,clatt(2),clatt(3)]+b3;
            catom(lask2+3+N_bulk+N_liquid,:)=[clatt(1)+Lshift,clatt(2),clatt(3)]+b4;
            type((lask2+N_bulk+N_liquid):(lask2+3+N_bulk+N_liquid))=3;
            lask2=lask2+4;
            
        end
    end
end

% Create the liquid atoms by random positioning
c_liquid=[rand(N_liquid,1)*L_liquid+L_bulk-xlatt/4+L_spare,rand(N_liquid,1)*Ly,rand(N_liquid,1)*Lz];

catom((N_bulk+1):(N_bulk+N_liquid),:)=c_liquid;
type((N_bulk+1):(N_bulk+N_liquid))=2;

% Some shifts of positions
catom(:,1)=catom(:,1)+xlatt/4;
catom(:,2)=catom(:,2)+ylatt/4;
%zdiff=c-(.5+u)*a3(3);
catom(:,3)=catom(:,3)+zlatt/4;
% N=size(catom,1);


figure(2222);clf;
set(gca,'fontsize',24)
%plot3(catom(inds1,1),catom(inds1,2),catom(inds1,3),'bo','markerfacecolor','b');
scale=1;
plot3(catom(type==1,1)*scale,catom(type==1,2)*scale,catom(type==1,3)*scale,'bo','markerfacecolor','b');
hold on
plot3(catom(type==2,1)*scale,catom(type==2,2)*scale,catom(type==2,3)*scale,'ro','markerfacecolor','r');
plot3(catom(type==3,1)*scale,catom(type==3,2)*scale,catom(type==3,3)*scale,'go','markerfacecolor','g');
%plot3(clatt(1:end,1),clatt(1:1:end,2),clatt(1:end,3),'bo','markerfacecolor','b');
%hold on
%plot3(catom(inds2,1),catom(inds2,2),catom(inds2,3),'ro','markerfacecolor','r')
axis equal
set(gca,'xlim',[0,Lx]*scale)
set(gca,'ylim',[0,Ly]*scale)
set(gca,'zlim',[0,Lz]*scale)
xlabel('x')
ylabel('y')
% return

fprintf('System contains %d atoms.\n',N);
fprintf('rho=%.2f, L_liq=%.2f.\n',rho,L_liquid);
N_latt=N;

tags=type;

% Start writing to file
% filepath='/proj/quantum-data/Kimmo/lammps/liquid-solid/';
filepath='/wrk/ksaaskil/lammps/liquid-solid/';

filename=strcat(filepath,filetowrite);
fid=fopen(filename,'w');

% The initializations
header=sprintf('# This is the fcc solid/liquid interface, rho=%.2f, L_liq=%.2f.\n',rho,L_liquid);
fprintf(fid,header);
fprintf(fid,'%d atoms\n',N);
fprintf(fid,'%d atom types\n',max(tags));
fprintf(fid,'\n');

% Box dimensions

% Additional box for safety, if no periodic boundary conditions are applied

if ~pbc
    boxedges=1;
else
    boxedges=0;
end

fprintf(fid,'%.2f %.2f xlo xhi\n',0,Lx);
fprintf(fid,'%.2f %.2f ylo yhi\n',-boxedges,Ly+boxedges);
fprintf(fid,'%.2f %.2f zlo zhi\n',-boxedges,Lz+boxedges);
fprintf(fid,'\n');

% The atom masses
fprintf(fid,'Masses\n\n');
M=ones(max(tags),1);

for i=1:max(tags)
    fprintf(fid,'%d %.3f\n',i,M(i));
end 
fprintf(fid,'\n');

% Write the atom coordinates and types
fprintf(fid,'Atoms\n\n');

for i=1:N
   fprintf(fid,'%d %d %.5f %.5f %.5f\n',i,tags(i),catom(i,1),catom(i,2),catom(i,3));
end
fprintf(fid,'\n');

fclose(fid);

fprintf('Wrote to file %s.\n',filename);

return
%%
% Start writing to file
% filepath='~/Documents/koodit/lammps/';

filename=strcat(filepath,filetowrite,'.kmask');
fid=fopen(filename,'w');

% The initializations

fprintf(fid,'# This is the wavevector mask file. id, kx, ky.');
fprintf(fid,'\n');

% Box dimensions

% Additional box for safety, if no periodic boundary conditions are applied

% For the two first unit cells
for i=1:8*Ny*Nz
    if tags(i)==1
        fprintf(fid,'%d %d %d',i,kmask(i,1),kmask(i,2));
        fprintf(fid,'\n');
    end
end

fclose(fid);

fprintf('Wrote to file %s.\n',filename);