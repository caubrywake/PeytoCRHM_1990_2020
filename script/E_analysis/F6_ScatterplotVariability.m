%% F6 - Scatter plot

 % Correlation between the annual and sesonal variables and streamflow
addpath('D:\4_PeytoCRHM_1990_2020')
clear all
close all

%% Import CRHM results
load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\HRUelev_area.mat')
load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\output_20210323.mat', 'SWE', 'hru_rain', 'hru_snow', 'hru_t', 'time')
load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\output_20210323.mat', 'basinflow', 'basingw')
load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\output_20210323.mat', 'SWE', 'SWEmelt', 'cumSWE_in', 'cumSWE_out', 'firnmelt', 'hru_actet', 'hru_drift', 'hru_rain', 'hru_snow', 'hru_subl', 'icemelt', 'infil', 'meltrunoff', 'runoff', 'snowinfil', 'time')

savedir = 'D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\'
figdir = 'D:\4_PeytoCRHM_1990_2020\fig\analysis\'

swemelt = SWEmelt;
hru_elev = hruelev;
hru_area = hruarea;
snow = hru_snow; 
rain = hru_rain;
rainfallrunoff = (runoff + infil) -icemelt/24  -firnmelt/24;
rainfallrunoff(rainfallrunoff<0)=0;
swemelt = SWEmelt;
Ta = hru_t;
snow = hru_snow;
rain = hru_rain;
precip = snow+rain;
timeCRHM = time;

yr = 1988:2020;
hru = 1:37;
for i = 1:length (yr)
    % Annual values
t1 = strcat ('01-Oct-', num2str(yr(i)-1),{' '}, '1:00');
t2 = strcat ('30-Sep-', num2str(yr(i)), {' '}, '1:00');  
a = find(timeCRHM==datetime(t1));
b = find(timeCRHM==datetime(t2));

x = precip;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
P(i, 1) = xc;
x = rain;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
R(i, 1) = xc;
x = snow;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
S(i, 1) = xc;
x = Ta+273.15;
xc  = sum(mean(x(a:b,hru)).*(hru_area(hru)/sum(hru_area))); 
T(i, 1) = xc-273.15;


% Fall
 t1 = strcat ('01-Sep-', num2str(yr(i)-1),{' '}, '1:00');
 t2 = strcat ('30-Nov-', num2str(yr(i)-1),{' '}, '1:00');
a = find(timeCRHM==datetime(t1));
b = find(timeCRHM==datetime(t2));
x = precip;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
P(i, 2) = xc;
x = rain;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
R(i, 2) = xc;
x = snow;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
S(i, 2) = xc;
x = Ta+273.15;
xc  = sum(mean(x(a:b,hru)).*(hru_area(hru)/sum(hru_area))); 
T(i, 2) = xc-273.15;

% Winter
 t1 = strcat ('01-Dec-', num2str(yr(i)-1),{' '}, '1:00');
 t2 = strcat ('28-Feb-', num2str(yr(i)),{' '}, '1:00');
a = find(timeCRHM==datetime(t1));
b = find(timeCRHM==datetime(t2));
x = precip;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
P(i, 3) = xc;
x = rain;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
R(i, 3) = xc;
x = snow;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
S(i, 3) = xc;
x = Ta+273.15;
xc  = sum(mean(x(a:b,hru)).*(hru_area(hru)/sum(hru_area))); 
T(i, 3) = xc-273.15;

% Spring
 t1 = strcat ('01-Mar-', num2str(yr(i)),{' '}, '1:00');
 t2 = strcat ('30-May-', num2str(yr(i)),{' '}, '1:00');
a = find(timeCRHM==datetime(t1));
b = find(timeCRHM==datetime(t2));
x = precip;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
P(i, 4) = xc;
x = rain;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
R(i, 4) = xc;
x = snow;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
S(i, 4) = xc;
x = Ta+273.15;
xc  = sum(mean(x(a:b,hru)).*(hru_area(hru)/sum(hru_area))); 
T(i, 4) = xc-273.15;

% Summer
 t1 = strcat ('01-Jun-', num2str(yr(i)),{' '}, '1:00');
 t2 = strcat ('30-Aug-', num2str(yr(i)),{' '}, '1:00');
a = find(timeCRHM==datetime(t1));
b = find(timeCRHM==datetime(t2));
x = precip;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
P(i, 5) = xc;
x = rain;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
R(i, 5) = xc;
x = snow;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
S(i, 5) = xc;
x = Ta+273.15;
xc  = sum(mean(x(a:b,hru)).*(hru_area(hru)/sum(hru_area))); 
T(i, 5) = xc-273.15;

% Components
 t1 = strcat ('01-Oct-', num2str(yr(i)-1),{' '}, '1:00');
 t2 = strcat ('30-Sep-', num2str(yr(i)),{' '}, '1:00');
a = find(timeCRHM==datetime(t1));
b = find(timeCRHM==datetime(t2));

x = swemelt/24;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
C(i, 1) = xc;

x = icemelt/24;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ; xc = sum(xb)/(sum(hru_area*10^6)); 
C(i, 2) = xc;

x = firnmelt/24;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ; xc = sum(xb)/(sum(hru_area*10^6)); 
C(i, 3) = xc;

x = rainfallrunoff;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ; xc = sum(xb)/(sum(hru_area*10^6)); 
C(i, 4) = xc;

% peak swe
 t1 = strcat ('01-Oct-', num2str(yr(i)-1),{' '}, '1:00');
 t2 = strcat ('30-Sep-', num2str(yr(i)),{' '}, '1:00');
a = find(timeCRHM==datetime(t1));
b = find(timeCRHM==datetime(t2));
C(i,5) = max(SWE(a:b));

% Annual Streamflo
Sf(i, 1) = sum(basinflow(a:b))/(sum(hru_area*10^6));% annual flow
Sf(i, 2) = max(basinflow(a:b)/3600);% peka flow
x =  find(basinflow(a:b)/3600 == max(basinflow(a:b)/3600));
y =timeCRHM(a:b);
Sf(i, 3) =day(y(x), 'dayofyear');% timing of peak
fiftyper = 0.5*Sf(i, 1);
x = cumsum((basinflow(a:b))/(sum(hru_area*10^6)));
y = find(x >= fiftyper, 1);
z = timeCRHM(a:b);
Sf(i, 4)= day(z(y), 'dayofyear');% triming of 50%
end

% Find fit and significance
 for i = 1:5
         
     [rhoa, pvala] = corrcoef(Sf(:, 1), P(:, i), 'rows', 'pairwise');
      rhoP(i) = rhoa(2); pvalP(i) = pvala(2);

       [rhoa, pvala] = corrcoef(Sf(:, 1), R(:, i), 'rows', 'pairwise');
      rhoR(i) = rhoa(2); pvalR(i) = pvala(2);
      
    [rhoa, pvala] = corrcoef(Sf(:, 1), S(:, i), 'rows', 'pairwise');
      rhoS(i) = rhoa(2); pvalS(i) = pvala(2);
      
           [rhoa, pvala] = corrcoef(Sf(:, 1), T(:, i), 'rows', 'pairwise');
      rhoT(i) = rhoa(2); pvalT(i) = pvala(2);
      
      [rhoa, pvala] = corrcoef(Sf(:, 1), C(:, i), 'rows', 'pairwise');
      rhoC(i) = rhoa(2); pvalC(i) = pvala(2);
     end 

 % Compiled MET in a table
 A = [rhoT;pvalT;  rhoP;pvalP;  rhoR;pvalR; rhoS; pvalS;]
%;rhoC; pvalC
 A = round(A, 4);
 rowname = {'Ta_r';'Ta_pval';'Precip_r';'Precip_pval';'Rainfall_r';'Rainfall_pval';'Snowfall_r';'Snowfall_pval'};
 colname = {'Variable', 'Annual','Fall','Winter','Spring','Summer'};
Correlation_significance_MET_Flow = table(rowname, A(:, 1),A(:, 2),A(:, 3), A(:, 4), A(:, 5));
Correlation_significance_MET_Flow.Properties.VariableNames = colname;
writetable( Correlation_significance_MET_Flow, strcat(savedir, 'Correlation_significance_MET_Flow.txt'))
save(strcat(savedir, 'Correlation_significance_MET_Flow.mat'), 'Correlation_significance_MET_Flow'); 

% Tab;e Flow component
 A2 = [rhoC; pvalC]
 A2 = round(A2, 4);
 rowname = {'Component_r';'Component_pval'};
 colname = {'Variable', 'Snowmelt','Icemelt','Firnmelt','RainfallRunoff','PeakSWE'}; 
Correlation_significance_Comp_Flow = table(rowname, A2(:, 1),A2(:, 2),A2(:, 3), A2(:, 4), A2(:, 5));
Correlation_significance_Comp_Flow.Properties.VariableNames = colname;
writetable( Correlation_significance_Comp_Flow, strcat(savedir, 'Correlation_significance_Comp_Flow.txt'))
save(strcat(savedir, 'Correlation_significance_Comp_Flow.mat'), 'Correlation_significance_Comp_Flow'); 


 %% Flow component
 
 %% make scatter plot of these
 % 3 rows: Temp, Preicp, Component + peak SWE
 % Annual, fall, winter, spring, summer, or snowmele,t icemel,t firn,
 % rainfall runoff, and peak swe
 % 
 close all
sz = 5;
fig = figure('units','inches','outerposition',[0 0 8 10]);
colormap(parula)
fs = 8;
bf = 0.2
ylab = {'Annual';'Fall';'Winter';'Spring';'Summer'}
ylabC = {'Snow melt';'Ice melt';'Firn melt';'Rain. Run.';'Peak SWE'}



x = Sf(:, 1); % streamflow
% Temperature, top row
% Temp
txtloc = [-1.5, -0.5, -6, -1, 8]
for num = 1:5
subplot (5,5,num)
y = T(:, num); 
c = linspace(yr(1), yr(end),length(x)); 
scatter(x,y, sz, c, 'filled' );lsline
 xlim ([1.0 2.2]); xticks([1.5:.5:2.5]); xticklabels ([])
%ylim ([y1 y2]); ylabel(ylab {num});
if num == 1
    ylabel ('Ta(^{\circ}C)')
else 
    ylabel ([])
end

if pvalT(num)<= 0.05
title (strcat(ylab(num),',', num2str(round(rhoT(num), 2))), 'FontSize', fs, 'FontWeight', 'bold');
else
title (strcat(ylab(num),',', num2str(round(rhoT(num), 2))), 'FontSize', fs, 'Fontweight','normal');
end 
grid on; box on
end 

% Precip
txtloc = [1.9, .9, .45,.55, .9]
for num = 1:5
subplot (5,5,5+num)
y = P(:, num); 
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );lsline
 xlim ([1.0 2.2]); xticks([1.5:.5:2.5]); xticklabels ([])
%ylim ([y1 y2]); ylabel(ylab {num});
ylabel(ylab{num});
if pvalP(num)<= 0.05
text (1.3, txtloc(num), strcat(num2str(round(rhoP(num), 2))), 'FontSize', fs, 'FontWeight', 'bold');
else
text (1.3,  txtloc(num), strcat(num2str(round(rhoP(num), 2))), 'FontSize', fs);
end 
grid on; box on
end 

% Rain
txtloc = [.9, .18,.0009, .09, .45]
for num = 1:5
subplot (5,5,10+num)
y = R(:, num); 
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );lsline
 xlim ([1.0 2.2]); xticks([1.5:.5:2.5]); 
%xticklabels ([])
%ylim ([y1 y2]); ylabel(ylab {num});
ylabel(ylab{num});
if pvalR(num)<= 0.05
text (1.3, txtloc(num), strcat(num2str(round(rhoR(num), 2))), 'FontSize', fs, 'FontWeight', 'bold');
else
text (1.3,  txtloc(num), strcat(num2str(round(rhoR(num), 2))), 'FontSize', fs);
end 
grid on; box on
end 

% Snow
txtloc = [1.8, 0.9,0.45, 0.55, 0.35]
for num = 1:5
subplot (5,5,15+num)
y = S(:, num); 
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );lsline
 xlim ([1.0 2.2]); xticks([1.5:.5:2.5]); xticklabels ([])
% ylim ([y1 y2]); ylabel(ylab {num});
ylabel(ylab{num});
if pvalS(num)<= 0.05
text (1.3, txtloc(num), strcat(num2str(round(rhoS(num), 2))), 'FontSize', fs, 'FontWeight', 'bold');
else
text (1.3,  txtloc(num), strcat(num2str(round(rhoS(num), 2))), 'FontSize', fs);
end 
grid on; box on
end 

% Component
txtloc = [1.4, 0.9,0.35, 0.21, 1900]
for num = 1:5
subplot (5,5,20+num)
y = C(:, num); 
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );lsline
 xlim ([1.0 2.2]); xticks([1.5:.5:2.5]); xticklabels ([])
%ylim ([y1 y2]); ylabel(ylab {num});
ylabel(ylabC{num});
if pvalC(num)<= 0.05
text (1.3, txtloc(num), strcat(num2str(round(rhoC(num), 2))), 'FontSize', fs, 'FontWeight', 'bold');
else
text (1.3,  txtloc(num), strcat(num2str(round(rhoC(num), 2))), 'FontSize', fs);
end 
grid on; box on
end 
cb = colorbar;
pos = cb.Position;
cb.Position = pos + [pos(1) pos(2) pos(3)+0.001 pos(4)];


%% Only show the significant ones
fig = figure('units','inches','outerposition',[0 0 8 6]);
fig.Renderer='Painters';

fs = 10
sz = 20
% row 1: annual ta, spring ta, Summer Ta
subplot (3,3,1)
y = T(:, 1); 
x = Sf(:, 1);
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );ls = lsline; ls.Color = 'k';
% xlim ([1 2.4]); xticks([1.0:.5:2.5]);
% ylim ([-5 -1]);
ylabel('Annual Ta (^{\circ}C)', 'FontSize', fs);
text (1.05, -2.3, strcat('r = ',  num2str(round(rhoT(1), 2))), 'FontSize', fs);
grid on; box on

subplot (3,3,2)
y = T(:, 3); 
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );ls = lsline; ls.Color = 'k';
% xlim ([1.2 2.4]); xticks([1.5:.5:2.5]); 
% ylim ([-17 -7]);
ylabel('Winter Ta (^{\circ}C)', 'FontSize', fs);
text (1.05, -6.5, strcat('r = ', num2str(round(rhoT(4), 2))), 'FontSize', fs);
grid on; box on

subplot (3,3,3)
y = T(:, 4); 
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );ls = lsline; ls.Color = 'k';
% xlim ([1 2.4]); xticks([1.0:.5:2.5]);
% ylim ([-10 -2]);yticks ([-10:2:0])
ylabel('Spring Ta (^{\circ}C)', 'FontSize', fs);
text (1.05, -1, strcat('r = ', num2str(round(rhoT(4), 2))), 'FontSize', fs);
grid on; box on
cb = colorbar;
cb.Position = cb.Position + ([0.08 0 0 0])

% Row 2: Summer Ta, Summer Rain, Winter snow
subplot (3,3,4)
y = T(:, 5); 
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );ls = lsline; ls.Color = 'k';
% xlim ([1 2.4]); xticks([1.0:.5:2.5]);
% ylim ([3 8]);
ylabel('Summer Ta (^{\circ}C)', 'FontSize', fs);
text (1.05, 7.5, strcat('r = ', num2str(round(rhoT(5), 2))), 'FontSize', fs);
grid on; box on

subplot (3,3,5)
y = R(:, 5); 
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );ls = lsline; ls.Color = 'k';
% xlim ([1 2.4]); xticks([1.0:.5:2.5]);
% ylim ([0 0.5]); yticks([0:0.25:0.5])
ylabel('Summer Rain (m)', 'FontSize', fs);
text (1.05, 0.54, strcat('r = ', num2str(round(rhoR(5), 2))), 'FontSize', fs);
grid on; box on

subplot (3,3,6)
y = S(:,3); 
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );ls = lsline; ls.Color = 'k';
% xlim ([1 2.4]); xticks([1.0:.5:2.5]);
% ylim ([0 1.2]);yticks ([0:0.5:1]);
ylabel('Winter Snow (m)', 'FontSize', fs);
text (1.05, 1.35, strcat('r = ', num2str(round(rhoS(3), 2))), 'FontSize', fs);
grid on; box on

% row 3 icemelt, firnmelt, rainfall runoff
subplot (3,3,7)
y = C(:, 2); 
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );ls = lsline; ls.Color = 'k';
% xlim ([1 2.4]); xticks([1.0:.5:2.5]);
% ylim ([0 1.5]);
ylabel('Icemelt (m w.e.)', 'FontSize', fs);
text (1.05, 1.38, strcat('r = ', num2str(round(rhoC(2), 2))), 'FontSize', fs);
grid on; box on

subplot (3,3,8)
y = C(:, 3); 
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );ls = lsline; ls.Color = 'k';
% xlim ([1 2.4]); xticks([1.0:.5:2.5]);
% ylim ([0 0.5]);yticks ([0:0.25:0.5])
ylabel('Firnmelt (m w.e.)', 'FontSize', fs);
xlabel ('Streamflow (m)', 'FontSize', fs);
text (1.05,0.55, strcat('r = ', num2str(round(rhoC(3), 2))), 'FontSize', fs);
grid on; box on

subplot (3,3,9)
y = C(:, 4); 
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );ls = lsline; ls.Color = 'k';
% xlim ([1 2.4]); xticks([1.0:.5:2.5]);
% ylim ([0 0.4]);
ylabel('Rainfall Runoff (m)', 'FontSize', fs);
text (1.05, 0.37, strcat('r = ', num2str(round(rhoC(4), 2))), 'FontSize', fs);
grid on; box on

figname ='SignficantMET_Flow_Scatterplot';
saveas (gcf, strcat(figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname,  '.png'))
savefig(gcf, strcat(figdir, figname))

%%%%%%%%%%%%%%%%%%%%%%
%%














%%

num = 2;
y1= -5.5; y2 = -1; t1 = -1.5;
subplot (3,5,2)
y = T(:, num); %annual air t
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );lsline
xlim ([1.2 2.7]); xticks([1.5:.5:2.5]); ylim ([y1 y2]); ylabel(ylab {num});
if rhoT(1) >= 0.005
text (1.3, t1, strcat(num2str(round(rhoT(num), 2))), 'FontSize', fs, 'FontWeight', 'bold');
else
text (1.3, t1, strcat(num2str(round(rhoT(num), 2))), 'FontSize', fs);
end 
grid on; box on

num = 3;
y1= -15; y2 = -9; t1 = -9.5;
subplot (3,5,3)
y = T(:, num); %annual air t
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );lsline
xlim ([1.2 2.7]); xticks([1.5:.5:2.5]); ylim ([y1 y2]); ylabel(ylab {num});
if rhoT(1) >= 0.005
text (1.3, t1, strcat(num2str(round(rhoT(num), 2))), 'FontSize', fs, 'FontWeight', 'bold');
else
text (1.3, t1, strcat(num2str(round(rhoT(num), 2))), 'FontSize', fs);
end 
grid on; box on

num = 4;
y1= -15; y2 = -9; t1 = -9.5;
subplot (3,5,3)
y = T(:, num); %annual air t
c = linspace(yr(1), yr(end),length(x)); scatter(x,y, sz, c, 'filled' );lsline
xlim ([1.2 2.7]); xticks([1.5:.5:2.5]); ylim ([y1 y2]); ylabel(ylab {num});
if rhoT(1) >= 0.005
text (1.3, t1, strcat(num2str(round(rhoT(num), 2))), 'FontSize', fs, 'FontWeight', 'bold');
else
text (1.3, t1, strcat(num2str(round(rhoT(num), 2))), 'FontSize', fs);
end 
grid on; box on




%%

num = 2;
subplot (3,5,2)
x = Sf; % streamflow
y = T(:, num); %annual air t
c = linspace(yr(1), yr(end),length(x));
scatter(x,y, sz, c, 'filled' );lsline
xlim ([1.2 2.7]); xticks([1.5:.5:2.5])
ylim ([min(y)-bf max(y)+bf]); 
ylabel('Fall T_a')
if pvalT(num) <= 0.005
text (1.3, max(y), strcat(num2str(round(rhoT(num), 2))), 'FontSize', fs, 'FontWeight', 'bold');
else
text (1.3, max(y), strcat(num2str(round(rhoT(num), 2))), 'FontSize', fs);
end 
grid on
box on

num = 3;
subplot (3,5,3)
x = Sf; % streamflow
y = T(:, num); %annual air t
c = linspace(yr(1), yr(end),length(x));
scatter(x,y, sz, c, 'filled' );lsline
xlim ([1.2 2.7]); xticks([1.5:.5:2.5])
ylim ([min(y)-bf max(y)+bf]); 
ylabel('Winter T_a')
if pvalT(num) <= 0.005
text (1.3, max(y), strcat(num2str(round(rhoT(num), 2))), 'FontSize', fs, 'FontWeight', 'bold');
else
text (1.3, max(y), strcat(num2str(round(rhoT(num), 2))), 'FontSize', fs);
end 
grid on
box on

num = 4;
subplot (3,5,4)
x = Sf; % streamflow
y = T(:, num); %annual air t
c = linspace(yr(1), yr(end),length(x));
scatter(x,y, sz, c, 'filled' );lsline
xlim ([1.2 2.7]); xticks([1.5:.5:2.5])
ylim ([min(y)-bf max(y)+bf]); 
ylabel('Spring T_a')
if pvalT(num) <= 0.005
text (1.3, max(y), strcat(num2str(round(rhoT(num), 2))), 'FontSize', fs, 'FontWeight', 'bold');
else
text (1.3, max(y), strcat(num2str(round(rhoT(num), 2))), 'FontSize', fs);
end 
grid on
box on

num = 5;
subplot (3,5,5)
x = Sf; % streamflow
y = T(:, num); %annual air t
c = linspace(yr(1), yr(end),length(x));
scatter(x,y, sz, c, 'filled' );lsline
xlim ([1.2 2.7]); xticks([1.5:.5:2.5])
ylim ([min(y)-bf max(y)+bf]); 
ylabel('Summer T_a')
if pvalT(num) <= 0.005
text (1.3, max(y), strcat(num2str(round(rhoT(num), 2))), 'FontSize', fs, 'FontWeight', 'bold');
else
text (1.3, max(y), strcat(num2str(round(rhoT(num), 2))), 'FontSize', fs);
end 
grid on
box on

%%
num = 3;
x = A2(:,5); % streamflow
y = A4(:,6);
c = linspace(yr(1), yr(end),length(x));
scatter(x,y, sz, c, 'filled' );lsline
xlim ([1.2 2.7])
xlabel ('Streamflow (mm)')
ylim ([-16 -8])
ylabel('Winter Ta')
text (1.3, -9, strcat('r = ', num2str(round(rho(6), 2))), 'FontSize', 10);
text (1.3,-9.8, strcat('p = ',num2str(round(pval(6), 3))), 'FontSize', 10);
grid on
box on

 subplot (3,4,3)
x = A2(:,5); % streamflow
y = A4(:,8);
c = linspace(yr(1), yr(end),length(x));
scatter(x,y, sz, c, 'filled' );lsline
xlim ([1.2 2.7])
xlabel ('Streamflow (mm)')
ylim ([-9 -2])
ylabel('Spring Ta')
text (1.3, -8, strcat('r = ', num2str(round(rho(8), 2))), 'FontSize', 10);
text (1.3,-7.2, strcat('p = ',num2str(round(pval(8), 3))), 'FontSize', 10);
grid on
box on

 subplot (3,4,4)
x = A2(:,5); % streamflow
y = A4(:,10);
c = linspace(yr(1), yr(end),length(x));
scatter(x,y, sz, c, 'filled' );lsline
xlim ([1.2 2.7])
xlabel ('Streamflow (mm)')
ylim ([4.2 8])
ylabel('Summer Ta')
text (1.3, 7.5, strcat('r = ', num2str(round(rho(10), 2))), 'FontSize', 10);
text (1.3,7, strcat('p = ', num2str(round(pval(10), 5))), 'FontSize', 10);
grid on
box on

% Second Row: Precip
subplot (3,4,5)
x = A2(:,5); % streamflow
y = A2(:,1);
c = linspace(yr(1), yr(end),length(x));
scatter(x,y, sz, c, 'filled' );lsline
xlim ([1.2 2.7])
xlabel ('Streamflow (mm)')
ylim ([.8 1.55])
ylabel('Snow melt')
text (1.3, 1.4, strcat('r = ', num2str(round(rho2(1), 2))), 'FontSize', 10);
text (1.3,1.3, strcat('p = ', num2str(pval2(1), 2)), 'FontSize', 10);
grid on
box on


 subplot (3,3,)
x = A2(:,5); % streamflow
y = A2(:,2);
c = linspace(yr(1), yr(end),length(x));
scatter(x,y, sz, c, 'filled' );lsline
xlim ([1.2 2.7])
xlabel ('Streamflow (mm)')
ylim ([0.2 1])
ylabel('Ice melt')
text (1.3, .9, strcat('r = ', num2str(round(rho2(2), 2))), 'FontSize', 10);
text (1.3,.82, strcat('p = ',num2str(round(pval2(2), 3))), 'FontSize', 10);
grid on
box on

 subplot (2,4,7)
x = A2(:,5); % streamflow
y = A2(:,3);
c = linspace(yr(1), yr(end),length(x));
scatter(x,y, sz, c, 'filled' );lsline
xlim ([1.2 2.7])
xlabel ('Streamflow (mm)')
ylim ([-0.0250 0.3])
ylabel('Firn melt')
text (1.3, .25, strcat('r = ', num2str(round(rho2(3), 2))), 'FontSize', 10);
text (1.3,.2, strcat('p = ',num2str(round(pval2(3), 3))), 'FontSize', 10);
grid on
box on

 subplot (2,4,8)
x = A2(:,5); % streamflow
y = A2(:,4);
c = linspace(yr(1), yr(end),length(x));
scatter(x,y, sz, c, 'filled' );lsline
xlim ([1.2 2.7])
xlabel ('Streamflow (mm)')
ylim ([0 0.25])
ylabel('Rainfall Runoff')
text (1.3, 0.24, strcat('r = ', num2str(round(rho2(4), 2))), 'FontSize', 10);
text (1.3,0.2, strcat('p = ', num2str(round(pval2(4), 5))), 'FontSize', 10);
grid on
box on

filename = 'D:\PeytoVariabilityPaper\FIG\Scatterplot_streamflow_Mar1'
%%
savefig (filename);
saveas(gcf,strcat([filename, '.png']))
saveas(gcf,strcat([filename, '.eps']))
 