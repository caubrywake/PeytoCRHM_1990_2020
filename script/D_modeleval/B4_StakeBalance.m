%% Mass balance stakes, with dates
close all
clear all

addpath(genpath('D:\PeytoCRHM_1990_2020\chrm'))
addpath('D:\PeytoCRHM_1990_2020\function')
savedir = 'D:\PeytoCRHM_1990_2020\data_process\modeleval\'
figdir = 'D:\PeytoCRHM_1990_2020\fig\modeleval\'
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\SWE.mat', 'SWE', 'timeCRHM')
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\HRUelev_area.mat', 'hruelev')
load('D:\PeytoCRHM_1990_2020\data_raw\massbalance\StakePeyto_2005_2020_withdates.mat')
Stake = table2array(StakePeyto20052020withdates(:, 1:8));
Stake_bwdate = table2array(StakePeyto20052020withdates(:, 9));
Stake_bsdate = table2array(StakePeyto20052020withdates(:, 10));

 % fix a error
 Stake(412, 5)=2165;
%% First, plot each stake to see its duration and elevation
 load('D:\PeytoCRHM_1990_2020\data_raw\massbalance\Stake_withcrrespondingHRU.mat')
Stake_CorrHRU(Stake_CorrHRU==0)=nan;
% grou stake in 3 vbatch: lower, north and south
stknorth = [190, 11, 180, 12, 13, 175, 171, 165];
stksouth =[220, 210, 140, 130, 120, 110, 100];
stklow = [160, 85, 805, 80, 82, 81, 60, 62, 51, 50, 52, 41, 40, 42];
    figure
    subplot(3,1,1)
for i = 1:length(stklow)
a = find(Stake(:, 2) == stklow(i));
plot(Stake_bwdate(a), Stake(a, 5), '-x'); hold on
end 
legend ('160', '85', '805','80',' 82', '81', '60', '62',' 51', '50',' 52',' 41', '40',' 42')
  ylim([2100 2350])
  
subplot(3,1,2)
for i = 1:length(stksouth)
a = find(Stake(:, 2) == stksouth(i));
plot(Stake_bwdate(a), Stake(a, 5), '-x'); hold on
end 
legend('220', '210', '140', '130', '120', '110',' 100');
  ylim([2100 2700])

    subplot(3,1,3)
    figure
for i = 1:length(stknorth)
a = find(Stake(:, 2) == stknorth(i));
plot(Stake_bwdate(a), Stake(a, 5), '-x'); hold on
end 
legend ('190', '11', '180', '12', '13', '175', '171', '165')
% 


 figure
 subplot(2,1,1)
a = find(Stake(:, 2) == 160);
plot(Stake_bwdate(a), Stake(a, 5), '-xk'); hold on
a = find(Stake(:, 2) == 165);
plot(Stake_bwdate(a), Stake(a, 5), '-xb'); hold on
a = find(Stake(:, 2) == 171);
plot(Stake_bwdate(a), Stake(a, 5), '-xr'); hold on
a = find(Stake(:, 2) == 175);
plot(Stake_bwdate(a), Stake(a, 5), '-xc'); hold on
subplot(2,1,2)
a = find(Stake(:, 2) == 160);
plot(Stake_bwdate(a), Stake(a, 6), '-xk'); hold on
a = find(Stake(:, 2) == 165);
plot(Stake_bwdate(a), Stake(a, 6), '-xb'); hold on
a = find(Stake(:, 2) == 171);
plot(Stake_bwdate(a), Stake(a, 6), '-xr'); hold on
a = find(Stake(:, 2) == 175);
plot(Stake_bwdate(a), Stake(a, 6), '-xc'); hold on


a = find(Stake(:, 2) == 190);
plot(Stake_bwdate(a),Stake(a, 5), '-xk'); hold on
a = find(Stake(:, 2) == 180);
plot(Stake_bwdate(a),Stake(a, 5), '-ok')
a = find(Stake(:, 2) == 175);
plot(Stake_bwdate(a),Stake(a, 5), '-*k')

%% corresponding HRU
% Per Stake
 load('D:\PeytoCRHM_1990_2020\data_raw\massbalance\Stake_withcrrespondingHRU.mat')
Stake_CorrHRU(Stake_CorrHRU==0)=nan;

%
for i = 1:length(Stake_CorrHRU);
    subplot(5,6,i)
a = find(Stake(:, 2) == Stake_CorrHRU(i,1));
plot(Stake_bwdate(a),Stake(a, 6), '-k'); hold on
hru = Stake_CorrHRU(i,2:3);
hru(isnan(hru))=[];
title (strcat('Stake', num2str(Stake_CorrHRU(i,1))))
for ii = 1:length(a)
b(ii) = find(timeCRHM== Stake_bwdate(a(ii)));
b= sort(b);
end 
sc3 = plot(timeCRHM(b), SWE(b, hru)./10,'-r');
title (strcat(num2str(Stake_CorrHRU(i,1)), ',', num2str(hru)))
end 

%% Per hru
load('D:\PeytoCRHM_1990_2020\data_raw\massbalance\HRU_correspondingStake.mat')
for i = 1:length(HRU_corrStake)
    hru = HRU_corrStake(i, 1);
    hrucorr_elev(i,1) = hruelev(hru)
end
[sortedhru,strid] = sort(hrucorr_elev, 'descend')

HRU_corrStake(HRU_corrStake==0)=nan;
HRU_corrStake=HRU_corrStake(strid, :)
% winter
%
close all
fig = figure('units','inches','outerposition',[0 0 8 6]);
for ihru = 1:length(HRU_corrStake)
subplot(4,3,ihru)
clear crhm_swe a crhm_idx stkelev
stk = HRU_corrStake(ihru, 2:end);
stk(isnan(stk))=[];
stk_numel=numel(stk);
for istake = 1:numel(stk)

    b = find(Stake(:, 2) ==  stk(istake));
    stkelev(istake) = nanmean(Stake(b, 5));
stkelev(isnan(stkelev))=[];  
stkelev = round(stkelev);
    a = find(Stake(:, 2) ==  HRU_corrStake(ihru, istake+1));
datea = Stake_bwdate(a);
for idate =1:length(datea) % fint the corresponding date for CRHM
 x = find(timeCRHM == datea(idate));
 if isempty(x) % if nat, replace by nan
     x = nan;
     crhm_swe(idate) = nan;
 else 
   crhm_swe(idate) = SWE(x, HRU_corrStake(ihru))   ;
 end 
end
p1 = plot(datea, Stake(a, 6)/100, '-dk', 'Markersize', 2, 'MarkerfaceColor','k'); hold on
p2 = plot(datea, crhm_swe/1000, '-bo', 'Markersize', 4', 'MarkerfaceColor','b');
 ylim([0 2.5]);
 yticks ([0:0.5:2])
end 

%text(datetime('01-Jun-2003'), 1900,strcat('HRU:',num2str(HRU_corrStake(ihru,1)),', Stakes:', num2str(stk)))
text(datetime('15-May-2003'), 2.25, strcat('HRU:',num2str(sortedhru(ihru,1)), 'm'), 'Color', 'b', 'Fontsize', 9)
text(datetime('15-May-2010'), 2.25, strcat('Stakes:', num2str(round(mean(stkelev))),'m'), 'Fontsize', 9);
%title(strcat('HRU:',num2str(HRU_corrStake(ihru,1)),', Stakes:', num2str(stk)), 'Fontsize', 8);

grid on
box on

end
% jointfig(fig,3,4)
% [ax1, h1] = suplabel({'Years'});
% ax1.Position = ax1.Position + [0 +0.04 0 0];
[ax2, h2] = suplabel('Winter Mass Balance (m. w.e.)', 'y');
ax2.Position = ax2.Position + [+0.04 0 0 0];

  leg = legend ([p1 p2], 'Measured', 'Modelled', ...
      'Orientation', 'horizontal', 'Location', [0.267 0.95 0 0] ); 
 %
 tightfig(fig)

figname ='Winter_Stake_PerHRU';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))


%% For summer
% 
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\HRUelev_area.mat', 'hruelev')
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\firn_ice.mat')
 HRU_corrStakesum =  HRU_corrStake;
  HRU_corrStakesum ([8,12], :)=[];
% yr = [2003:2009, 2012:2015];
ice = ice+firn;
fig = figure('units','inches','outerposition',[0 0 8 6]);
for ihru = 1:length(HRU_corrStakesum)
subplot(4,3,ihru)
clear crhm_spr crhm_fall x xvec a crhm_idx stkelev
stk = HRU_corrStakesum(ihru, 2:end);
stk(isnan(stk))=[];
stk_numel=numel(stk);
for istake = 1:numel(stk)
        b = find(Stake(:, 2) ==  stk(istake));
    stkelev(istake) = nanmean(Stake(b, 5));
stkelev(isnan(stkelev))=[];  
stkelev = round(stkelev);

a = find(Stake(:, 2) ==  HRU_corrStakesum(ihru, istake+1));
datea = Stake_bsdate(a);
for idate =1:length(datea) % fint the corresponding date for CRHM
 x = find(timeCRHM == datea(idate));
 xvec = datevec(Stake_bsdate(a(idate)));
 xvecyr = xvec(:, 1);

 if isempty(x) % if nat, replace by nan
     x = nan;
     crhm_fall(idate) = nan;
     crhm_spr(idate) = nan;
 else   
     crhm_fall(idate) = ice(x, HRU_corrStakesum(ihru)); % melt at that dat
     b = find(timeCRHM== datetime(strcat('1-Apr-', num2str(xvecyr))));
     crhm_spr(idate) = ice(b, HRU_corrStakesum(ihru));
end 
end

summer = crhm_fall- crhm_spr;
summer(summer==0) = nan;
p1 = plot(datea, Stake(a, 7)/100, '-dk', 'Markersize', 2, 'MarkerfaceColor','k'); hold on
p2 = plot(datea, summer/1000, '-ro', 'Markersize', 4', 'MarkerfaceColor','r');

% ylim([400 2100]);
end 
%text(datetime('01-Jun-2003'), 1900,strcat('HRU:',num2str(HRU_corrStake(ihru,1)),', Stakes:', num2str(stk)))
%title(strcat('HRU:',num2str(HRU_corrStake(ihru,1)),', Stakes:', num2str(stk)), 'Fontsize', 8);
grid on
box on
ylim([-6.5 1.5])

text(datetime('15-Oct-2003'), 0.9, strcat('HRU:',num2str(sortedhru(ihru,1)), 'm'), 'Color', 'r', 'Fontsize', 9)
text(datetime('15-Oct-2010'), 0.9, strcat('Stakes:', num2str(round(mean(stkelev))),'m'), 'Fontsize', 9);

end

% jointfig(fig,3,4)
% [ax1, h1] = suplabel({'Years'});
% ax1.Position = ax1.Position + [0 +0.04 0 0];
[ax2, h2] = suplabel('Summer Mass Balance (m. w.e.)', 'y');
ax2.Position = ax2.Position + [+0.04 0 0 0];

  leg = legend ([p1 p2], 'Measured', 'Modelled', ...
      'Orientation', 'horizontal', 'Location', [0.267 0.95 0 0] ); 
 %
 tightfig(fig)
figname ='Summer_Stake_PerHRU';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

%%
%%
clear crhm_swe a crhm_idx
a = find(Stake(:, 2) == 11);
for i =1:length(a)
 x = find(timeCRHM == Stake_bwdate(a(i)));
 if isempty(x)
     x = nan;
     crhm_swe(i) = nan;
 else 
   crhm_swe(i) = SWE(x, HRU_corrStake(1,i))   
 end 
end
plot(Stake_bwdate(a), Stake(a, 6)*10, '-xb'); hold on
plot(Stake_bwdate(a), crhm_swe, '-xk')

end
%%



clear crhm_swe a crhm_idx
a = find(Stake(:, 2) == 11);
for i =1:length(a)
crhm_idx(i,1) = find(timeCRHM == Stake_bwdate(a(i)));
end
crhm_swe = SWE(crhm_idx, HRU_corrStake(1,1))
plot(Stake_bwdate(a), Stake(a, 6)*10, '-db'); hold on
plot(timeCRHM(crhm_idx), crhm_swe, '-dk')
legend ('stake 190', 'HRU 3','stake 11','HRU3')

%% For hru 4
figure
a = find(Stake(:, 2) == 190);
for i =1:length(a)
crhm_idx(i,1) = find(timeCRHM == Stake_bwdate(a(i)));
end
crhm_swe = SWE(crhm_idx, HRU_corrStake(2,1))
plot(Stake_bwdate(a), Stake(a, 6)*10, '-xb'); hold on
plot(timeCRHM(crhm_idx), crhm_swe, 'xk')

clear crhm_swe a crhm_idx
a = find(Stake(:, 2) == 11);
for i =1:length(a)
crhm_idx(i,1) = find(timeCRHM == Stake_bwdate(a(i)));
end
crhm_swe = SWE(crhm_idx, HRU_corrStake(2,1))
plot(Stake_bwdate(a), Stake(a, 6)*10, 'db'); hold on
plot(timeCRHM(crhm_idx), crhm_swe, 'dk')

clear crhm_swe a crhm_idx
a = find(Stake(:, 2) == 180);
for i =1:length(a)
crhm_idx(i,1) = find(timeCRHM == Stake_bwdate(a(i)));
end
crhm_swe = SWE(crhm_idx, HRU_corrStake(2,1))
plot(Stake_bwdate(a), Stake(a, 6)*10, 'ob'); hold on
legend ('stake 190', 'HRU 4','stake 11','HRU4', 'stake 180', 'hru 4')

%%
a = find(Stake(:, 2) == 171);
plot(Stake_bwdate(a), Stake(a, 6), '-xr'); hold on
a = find(Stake(:, 2) == 175);
plot(Stake_bwdate(a), Stake(a, 6), '-xc'); hold on












for i = 1:length(HRU_corrStake)
    
  
for ii = 1:length(stk)
% find allt the stakes fat correpon to that HRU
plot(Stake_bwdate(a),Stake(a, 6), '-xk'); hold on
for iii = 1:length(a)
b(iii) = find(timeCRHM== Stake_bwdate(a(iii)));
end
b= sort(b);
sc3 = plot(timeCRHM(b), SWE(b, HRU_corrStake(i,1))./10,'-xg');
clear b
end 
end 

%% Per year
% for 2003




%
%%
fig = figure('units','inches','outerposition',[0 0 8 6]);
yr = [2003:2009, 2011:2017];
hru = [1:16, 18:21, 34];
 
for i = 1:length(yr)
    subplot(5,4,i)
% GSC
a = find(Stake(:, 1) == yr(i));
sc1 = scatter(Stake(a, 5), Stake(a, 6)/100,10, '^', 'k'); hold on; 

% WGMS
a = find(WGMS_winter(:, 1) == yr(i));
sc1 = scatter(WGMS_winter(a, 2), WGMS_winter(a, 3)/1000,10, '^', 'k'); hold on; 

 % HRU CRHM
 b = find(timeCRHM== Stake_bwdate(a(1));
            sc3 = scatter(hruelev(hru), SWE(b, hru)/1000,10, 'd', 'filled', 'b'); hold on
xlim([2100 3000]);
ylim([0 2.2]);
yticks([0.5:0.5:2.2]);
xticks ([2100:250:3050]);

text(2130, 2, num2str(yr(i)))
 a = get(gca,'XTickLabel');  
 set(gca,'XTickLabel',a,'fontsize',10);
  a = get(gca,'YTickLabel');  
  set(gca,'YTickLabel',a,'fontsize',10);
grid on
box on
set(gca, 'FontSize', 10');

end 
