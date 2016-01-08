% Create a solid-liquid-solid system
clear
filetowrite='081015f_C_H2O.dat';
% mode='bulk';
mode='interface';

pbc=1;

ML=1;
MR=1;


% alatt=1.5496; % Lattice constant in units of sigma for solid argon at 0 K
alatt=3.57;
% a=1; % Lattice constant in units of sigma for solid argon at 40 K
%alatt=1.579; % Multiplier for final multiplication of coordinates and box sizes
% sigma=3.4; % Solid argon
%a=1; 

%b1=0; % The single atom

% The basis vectors of the cubic unit cell
xlatt=3*sqrt(3)/3*alatt;
ylatt=alatt/sqrt(2);
zlatt=alatt*sqrt(3)/2/sqrt(2);

% How many unit cells (each containing 8 atoms since we employ an
% orthogonal basis) are included for each bulk
Nx=12;
Ny=12;%Ny=1;
Nz=Ny;%Nz=1;

Ly=Ny*ylatt;
Lz=Nz*zlatt;

% Lengths of the bulk parts
L_bulk=Nx*xlatt;

L_spare=2;
%L_spare=0;

% Length of the liquid part
% L_liquid=Nx*alatt;
L_liquid=50;
%L_liquid=0;
% Volume of liquid atoms
V_liquid=L_liquid*Ly*Lz;
% Liquid density
% rho=1.03e3; % kg/m^3
% Convert to the number of water molecules
% N_water=round(rho*V_liquid*1e-30/(16+2*1)/1.66e-27)
% return
% Number of particles
% N_liquid=3*N_water;
% Number of atoms in each bulk
% N_bulk=4*Nx*Ny*Nz;

if 1
   C_file=sprintf('diamond_%dx%d.xyz',Ny,Nz)
   file=strcat('/wrk/ksaaskil/lammps/liquid-solid/',C_file);
   fid=fopen(file,'r');
   A=textscan(fid,'%s%f%f%f','headerlines',2);
   fclose(fid);
   c_C=[A{2},A{3},A{4}];
   % c_liquid(:,1)=c_liquid(:,1)+L_bulk-xlatt/4;
   
   N_unit=size(c_C,1);
   N_bulk=Nx*N_unit;
   c_solid1=zeros(N_bulk,3);
   c_solid2=zeros(N_bulk,3);
   for k=1:Nx
       shiftx1=(k-1)*xlatt*repmat([1,0,0],N_unit,1);
       shiftx2=((k-1)*xlatt+L_bulk+L_spare+L_liquid+L_spare)*repmat([1,0,0],N_unit,1);
       c_solid1((k-1)*N_unit+(1:N_unit),:)=c_C+shiftx1;
       c_solid2((k-1)*N_unit+(1:N_unit),:)=c_C+shiftx2;
   end
   
   
end

if 1
   waterfile=sprintf('water_diamond_%dx%d_%dnm.xyz',Ny,Nz,L_liquid/10) 
   %waterfile=sprintf('water_fcc111_%dnm.xyz',1) 
   file=strcat('/wrk/ksaaskil/lammps/liquid-solid/',waterfile);
   fid=fopen(file,'r');
   A=textscan(fid,'%s%f%f%f','headerlines',2);
   fclose(fid);
   c_liquid=[A{2},A{3},A{4}];
   c_liquid(:,1)=c_liquid(:,1)+L_bulk+L_spare-0*xlatt/4;
   %c_liquid(:,2)=0.95*c_liquid(:,2)+ylatt/4;
   %c_liquid(:,3)=0.95*c_liquid(:,3)+zlatt/4;
   %c_liquid=[]
   N_liquid=size(c_liquid,1);
   N_water=size(c_liquid,1)/3;
   counter=1;
   mol_ids_liquid=zeros(N_water,1);
   bonds=zeros(2*N_water,2);
   angles=zeros(N_water,3);
   for i=1:N_water
        mol_ids_liquid(counter:counter+2)=i;
        bonds(2*(i-1)+1,:)=N_bulk+[counter,counter+1];
        bonds(2*(i-1)+2,:)=N_bulk+[counter,counter+2];
        angles(i,:)=N_bulk+[counter+1,counter,counter+2];
        counter=counter+3;
   end
end




% Length of the system
Lx=L_bulk+L_liquid+L_bulk+2*L_spare;


N=2*N_bulk+N_liquid;

% Atom positions
%catom=zeros(N,3);
catom=[c_solid1;c_liquid;c_solid2];

type=ones(N,1);

%kmask=zeros(4*Nx*Ny*Nz,2);

lask=0;
lask2=1;

kmask=zeros(N,2);


% Add the liquid
type((N_bulk+1):3:(N_bulk+N_liquid))=3; % O
type((N_bulk+2):3:(N_bulk+N_liquid))=4; % H
type((N_bulk+3):3:(N_bulk+N_liquid))=4; % H

type(N_bulk+N_liquid+1:end)=2;

% return

% Some shifts of positions
%catom(:,1)=catom(:,1)+xlatt/4;
%catom(:,2)=catom(:,2)+ylatt/4;
%zdiff=c-(.5+u)*a3(3);
%catom(:,3)=catom(:,3)+zlatt/4;
% Just for checks
%c_liquid(:,1)=c_liquid(:,1)+xlatt/4;


figure(2222);clf;
set(gca,'fontsize',24)
%plot3(catom(inds1,1),catom(inds1,2),catom(inds1,3),'bo','markerfacecolor','b');
scale=1;
plot3(catom(type==1,1)*scale,catom(type==1,2)*scale,catom(type==1,3)*scale,'bo','markerfacecolor','b');
hold on
plot3(catom(type==2,1)*scale,catom(type==2,2)*scale,catom(type==2,3)*scale,'go','markerfacecolor','g');
plot3(catom(type==3,1)*scale,catom(type==3,2)*scale,catom(type==3,3)*scale,'ro','markerfacecolor','r');
plot3(catom(type==4,1)*scale,catom(type==4,2)*scale,catom(type==4,3)*scale,'ko','markerfacecolor','k');
%plot3(clatt(1:end,1),clatt(1:1:end,2),clatt(1:end,3),'bo','markerfacecolor','b');
%hold on
%plot3(catom(inds2,1),catom(inds2,2),catom(inds2,3),'ro','markerfacecolor','r')
axis equal
set(gca,'xlim',[0,Lx]*scale)
%set(gca,'ylim',[0,Ly]*scale)
%set(gca,'zlim',[0,Lz]*scale)
xlabel('x')
ylabel('y')

% return
%fprintf('System contains %d atoms.\n',N);

N_latt=N;

c.x=catom(:,1);
c.y=catom(:,2);
c.z=catom(:,3);
c.xlo=0;c.xhi=Lx;
c.ylo=0;c.yhi=Ly;
c.zlo=0;c.zhi=Lz;

c.masses=[12,12,16,1];
c.types=type;
c.charges=-0.8472*(c.types==3)+0.4236*(c.types==4);
c.ids=(1:N)';
c.mol_ids=zeros(N,1);
c.mol_ids(c.types==3|c.types==4)=mol_ids_liquid;
c.bonds=bonds;
c.bondtypes=ones(size(c.bonds,1),1);
c.angles=angles;
c.angletypes=ones(size(c.angles,1),1);
filepath='/wrk/ksaaskil/lammps/liquid-solid/';

filename=strcat(filepath,filetowrite);

write_lammps_data(filename,c,'full');

return
