
%clear
%file='220514c';
%file='220714b';
% file='210714c4';
%file='101014';
%file='211214b';
file='061015g';
%load(strcat('../koodit/lammps/',file,'_avepos_tar/',file,'_avepos.mat'));
%load(strcat('/wrk/ksaaskil/lammps/',file,'_anharm_tar/',file,'_anharm.mat'),...
%    'Jom_ave1','Jom_ave2','oms_fft','DoS_ave1');
%load(strcat('/wrk/ksaaskil/lammps/',file,'_tar/',file,'_2.mat'));
%load(strcat('/wrk/ksaaskil/lammps/',file,'.mat'));
%load(strcat('/proj/quantum-data/Kimmo/lammps/',file,'_tar/',file,'_Tom.mat'));
%load(strcat('/wrk/ksaaskil/lammps/',file,'.mat'));
path='/wrk/ksaaskil/lammps/liquid-solid/';
filename=strcat(path,file,'_gofr.dat');
%filename=strcat(path,file,'.Ti_equil.dat');
%filename=strcat(path,file,'.density_liquid.dat');
%filename=strcat(path,file,'.temp_liquid.dat');
%filename=strcat(path,file,'.temp_liquid.dat');
%load(strcat('/proj/quantum-data/Kimmo/lammps/',file,'_tar/',file,'_Tomega_alpha.mat'));

fid=fopen(filename,'r');

A=textscan(fid,'%f%f%f','headerlines',4);

fclose(fid);
xs=A{1};
ys=A{2};

% Turn kcal/mol to Kelvins
% ys=ys*4.2e3/6.022e23/1.38e-23;

% Shake correction
%ys(ys<245)=ys(ys<245)*3/2;

figure(4345);%clf;
hold on
plot(xs,ys/ys(end)+2,'b-','linewidth',2);
%set(gca,'xlim',[31,44])
% find(diff(xs)>
return
ind_turn=find(diff(A{1})<0,1,'first');
if ~isempty(ind_turn) 
    xs(ind_turn:ind_turn:end)=[];
    ys(ind_turn:ind_turn:end)=[];
    xs=xs(1:ind_turn-1);
    %return
    ys=reshape(ys,ind_turn-1,length(ys)/(ind_turn-1));
    for k=1:10:size(ys,2)
        plot(xs,ys(:,k),'bo-','linewidth',2);
        title(sprintf('k=%d/%d',k,size(ys,2)));
        pause(.1)
    end
end

%set(gca,'xlim',[42,60])

return
