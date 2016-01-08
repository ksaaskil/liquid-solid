
%clear
%file='220514c';
%file='220714b';
% file='210714c4';
%file='101014';
%file='211214b';
file='131115a';
%file='061015a';
%load(strcat('../koodit/lammps/',file,'_avepos_tar/',file,'_avepos.mat'));
%load(strcat('/wrk/ksaaskil/lammps/',file,'_anharm_tar/',file,'_anharm.mat'),...
%    'Jom_ave1','Jom_ave2','oms_fft','DoS_ave1');
%load(strcat('/wrk/ksaaskil/lammps/',file,'_tar/',file,'_2.mat'));
%load(strcat('/wrk/ksaaskil/lammps/',file,'.mat'));
%load(strcat('/proj/quantum-data/Kimmo/lammps/',file,'_tar/',file,'_Tom.mat'));
%load(strcat('/wrk/ksaaskil/lammps/',file,'.mat'));
path='/wrk/ksaaskil/lammps/liquid-solid/';
filename=strcat(path,file,'.params_steady.dat');
%filename=strcat(path,file,'.params_equil.dat');
%filename=strcat(path,file,'.density_liquid.dat');
%filename=strcat(path,file,'.temp_liquid.dat');
%filename=strcat(path,file,'.temp_liquid.dat');
%load(strcat('/proj/quantum-data/Kimmo/lammps/',file,'_tar/',file,'_Tomega_alpha.mat'));

fid=fopen(filename,'r');

A=textscan(fid,'%f%f%f%f%f%f','headerlines',1);

fclose(fid);
pxx=A{1};
pyy=A{2};
pzz=A{3};
xS=A{4}; % Left solid max
L=A{5};
rhos=A{6};

% Turn kcal/mol to Kelvins
% ys=ys*4.2e3/6.022e23/1.38e-23;

Px=mean(pxx(end/2:end));
Py=mean(pyy(end/2:end));
Pz=mean(pzz(end/2:end));
rho=mean(rhos(end/2:end))
% Shake correction
%ys(ys<245)=ys(ys<245)*3/2;
fprintf('Average pressure of liquid %.2f MPa.\n',Py*1.67e-21/(3.4e-10)^3/1e6);
fprintf('Average density %f (1/sigma^3).\n',rho);
fprintf('Average density of liquid %.2f kg/m^3.\n',rho*40*1.66e-27/(3.4e-10)^3);
fprintf('Average length of liquid %.2f LJ units.\n',mean(L(end/2:end)));
fprintf('Average length of liquid %.2f nm.\n',mean(L(end/2:end))*3.4e-10*1e9);
fprintf('Average position of right solid max %.2f.\n',mean(xS(end/2:end)));


figure(2345);clf;
hold on
plot(pxx,'b-','linewidth',2);
plot(pyy,'r-','linewidth',2);
plot(pzz,'g-','linewidth',2);
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
