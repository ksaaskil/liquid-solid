% Create a solid-liquid-solid system
clear
filetowrite='071015a_H2O.dat';
% mode='bulk';
mode='interface';


padding=1;

waterfile=sprintf('071015a_H2O.xyz') 
%waterfile=sprintf('water_fcc111_%dnm.xyz',1) 
file=strcat('/wrk/ksaaskil/lammps/liquid-solid/',waterfile);
fid=fopen(file,'r');
A=textscan(fid,'%s%f%f%f','headerlines',2);
fclose(fid);
c_liquid=[A{2},A{3},A{4}];
   
c_liquid(:,1)=c_liquid(:,1)+padding;
c_liquid(:,2)=c_liquid(:,2)+padding;
c_liquid(:,3)=c_liquid(:,3)+padding;

N_liquid=size(c_liquid,1);
N_water=size(c_liquid,1)/3;
counter=1;
mol_ids_liquid=zeros(N_water,1);
bonds=zeros(2*N_water,2);
angles=zeros(N_water,3);
for i=1:N_water
    mol_ids_liquid(counter:counter+2)=i;
    bonds(2*(i-1)+1,:)=[counter,counter+1];
    bonds(2*(i-1)+2,:)=[counter,counter+2];
    angles(i,:)=[counter+1,counter,counter+2];
    counter=counter+3;
end
% return
% Length of the system


catom=c_liquid;

Lx=max(catom(:,1))+padding;
Ly=max(catom(:,2))+padding;
Lz=max(catom(:,3))+padding;
N=size(catom(:,1),1);

type=ones(N,1);
N_bulk=0;
% Add the liquid
type((N_bulk+1):3:(N_bulk+N_liquid))=1; % O
type((N_bulk+2):3:(N_bulk+N_liquid))=2; % H
type((N_bulk+3):3:(N_bulk+N_liquid))=2; % H

c.x=catom(:,1);
c.y=catom(:,2);
c.z=catom(:,3);
c.xlo=0;c.xhi=Lx;
c.ylo=0;c.yhi=Ly;
c.zlo=0;c.zhi=Lz;

c.masses=[16,1];
c.types=type;
c.charges=-0.8472*(c.types==1)+0.4236*(c.types==2);
c.ids=(1:N)';
c.mol_ids=zeros(N,1);
c.mol_ids(c.types==1|c.types==2)=mol_ids_liquid;
c.bonds=bonds;
c.bondtypes=ones(size(c.bonds,1),1);
c.angles=angles;
c.angletypes=ones(size(c.angles,1),1);
filepath='/wrk/ksaaskil/lammps/liquid-solid/';

filename=strcat(filepath,filetowrite);

write_lammps_data(filename,c,'full');

return
