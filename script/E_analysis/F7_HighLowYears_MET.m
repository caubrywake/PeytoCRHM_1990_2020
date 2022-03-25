%% Figure 7_ High and low flow years analysisHigh vs low flow years

clear all
addpath('D:\4_PeytoCRHM_1990_2020')
savedir = 'D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\'
figdir = 'D:\4_PeytoCRHM_1990_2020\fig\analysis\'
%% Import CRHM results
load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\HRUelev_area.mat')
load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\HRUelev_area.mat')
load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\output_20210323.mat', 'SWE', 'hru_rain', 'hru_snow', 'hru_t', 'time')
load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\output_20210323.mat', 'basinflow')
load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\output_20210323.mat', 'SWEmelt', 'basinflow', 'basingw', 'firnmelt', 'icemelt', 'infil', 'runoff')
rainfallrunoff = (runoff + infil) -icemelt/24  -firnmelt/24;
rainfallrunoff(rainfallrunoff<0)=0;
swemelt = SWEmelt;

 snow = hru_snow;
rain = hru_rain;
precip = snow+rain;
Ta=hru_t;
timeCRHM = time;
hru_area=hruarea;
hru_elev = hruelev;
%% Compile annual values for annual and seasonal hydrometeorological variables

basinflow = basinflow+basingw;
yr = 1989:2020
for i = 1:length (yr)
   t1 = strcat ('01-Oct-', num2str(yr(i)-1), {' '}, '01:00');
   t2 = strcat ('30-Sep-', num2str(yr(i)), {' '}, '01:00');
 
a = find(timeCRHM==datetime(t1));
b = find(timeCRHM==datetime(t2));

x = basinflow;
xa  = sum(x(a:b))/sum(hru_area*10^6); %m w.e per hru
OUT(i, 1) = xa;
end

BFstd = nanstd(OUT(2:end))
BFmean = nanmean(OUT(2:end))
figure
plot (yr, OUT, '-xk')
refline(0, BFmean)
refline(0, BFmean+BFstd)
refline(0, BFmean-BFstd)

LowFlow = [1995, 1996, 1997, 2000, 2003, 2008];
HighFlow= [1992, 1994, 2006, 2013, 2015, 2016];

    for i = 1:length(LowFlow)
LFyr(i) = find(yr == LowFlow(i))
HFyr(i)= find(yr == HighFlow(i))
    end
%% Calculate low and high flow years
% average conditions for thos 7 years?
 % swemelt
 hru = 1:37
 x = Ta;
   clear Xl Xh
 for i = 1:length(LFyr)
t1 = strcat ('01-Oct-', num2str(LowFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(LowFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, 4); %m w.e per hru
xc = mean(xa, 2); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'daily', 'mean');
var1 = table2array(TT);
t1 = TT.Time;

if length(var1) == 366
    var1= [var1(1:59); var1(61:end)]
end 
Xl(:, i) = var1;

t1 = strcat ('01-Oct-', num2str(HighFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(HighFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, 4); %m w.e per hru
xc = mean(xa, 2); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'daily', 'mean');
var1 = table2array(TT);
t1 = TT.Time;

if length(var1) == 366
    var1= [var1(1:59); var1(61:end)]
end 
 
Xh(:, i) = var1;
 end 
 Xlm(:, 1) = mean(Xl, 2); % daily temp
 Xhm(:, 1) = mean(Xh, 2);
 
  x = precip; 
  clear Xl Xh
    for i = 1:length(LFyr)
t1 = strcat ('01-Oct-', num2str(LowFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(LowFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'monthly', 'sum');
var1 = table2array(TT);
t1 = TT.Time;

Xl(:, i) = var1;

t1 = strcat ('01-Oct-', num2str(HighFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(HighFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'monthly', 'sum');
var1 = table2array(TT);
t1 = TT.Time;


Xh(:, i) = var1;
    end 
  Xlprecip(:, 1) = mean(Xl, 2); 
 Xhprecip(:, 1) = mean(Xh, 2);
 
   x = rain;
  clear Xl Xh
    for i = 1:length(LFyr)
t1 = strcat ('01-Oct-', num2str(LowFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(LowFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'monthly', 'sum');
var1 = table2array(TT);
t1 = TT.Time;

Xl(:, i) = var1;

t1 = strcat ('01-Oct-', num2str(HighFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(HighFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'monthly', 'sum');
var1 = table2array(TT);
t1 = TT.Time;


Xh(:, i) = var1;
    end 
  Xlprecip(:, 2) = mean(Xl, 2)
 Xhprecip(:, 2) = mean(Xh, 2)
 
   x = snow;
  clear Xl Xh
    for i = 1:length(LFyr)
t1 = strcat ('01-Oct-', num2str(LowFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(LowFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'monthly', 'sum');
var1 = table2array(TT);
t1 = TT.Time;

Xl(:, i) = var1;

t1 = strcat ('01-Oct-', num2str(HighFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(HighFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'monthly', 'sum');
var1 = table2array(TT);
t1 = TT.Time;


Xh(:, i) = var1;
    end 
  Xlprecip(:, 3) = mean(Xl, 2);
 Xhprecip(:, 3) = mean(Xh, 2);
 
    x =Ta;
  clear Xl Xh
    for i = 1:length(LFyr)
t1 = strcat ('01-Oct-', num2str(LowFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(LowFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, 4); %m w.e per hru
xc = mean(xa, 2); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'monthly', 'mean');
var1 = table2array(TT);
t1 = TT.Time;

Xl(:, i) = var1;

t1 = strcat ('01-Oct-', num2str(HighFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(HighFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, 4); %m w.e per hru
xc = mean(xa, 2); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'monthly', 'mean');
var1 = table2array(TT);
t1 = TT.Time;


Xh(:, i) = var1;
    end 
  Xlprecip(:, 4) = mean(Xl, 2);
 Xhprecip(:, 4) = mean(Xh, 2);

   x = SWE;
     clear Xl Xh
     hru = 1:36
    for i = 1:length(LFyr)
t1 = strcat ('01-Oct-', num2str(LowFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(LowFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'daily', 'mean');
var1 = table2array(TT);
t1 = TT.Time;

if length(var1) == 366
    var1= [var1(1:59); var1(61:end)]
end 
Xl(:, i) = var1;

t1 = strcat ('01-Oct-', num2str(HighFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(HighFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'daily', 'mean');
var1 = table2array(TT);
t1 = TT.Time;

if length(var1) == 366
    var1= [var1(1:59); var1(61:end)]
end 
 
Xh(:, i) = var1;
    end 
  Xlm(:,2) = mean(Xl, 2);
 Xhm(:, 2) = mean(Xh, 2);
 
   x = swemelt/24;
      clear Xl Xh
    for i = 1:length(LFyr)
t1 = strcat ('01-Oct-', num2str(LowFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(LowFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'daily', 'sum');
var1 = table2array(TT);
t1 = TT.Time;

if length(var1) == 366
    var1= [var1(1:59); var1(61:end)]
end 
Xl(:, i) = var1;

t1 = strcat ('01-Oct-', num2str(HighFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(HighFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'daily', 'sum');
var1 = table2array(TT);
t1 = TT.Time;

if length(var1) == 366
    var1= [var1(1:59); var1(61:end)]
end 
 
Xh(:, i) = var1;
    end 
  Xlm(:, 3) = mean(Xl, 2);
 Xhm(:, 3) = mean(Xh, 2);
 
   x = icemelt/24;
      clear Xl Xh
    for i = 1:length(LFyr)
t1 = strcat ('01-Oct-', num2str(LowFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(LowFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'daily', 'sum');
var1 = table2array(TT);
t1 = TT.Time;

if length(var1) == 366
    var1= [var1(1:59); var1(61:end)]
end 
Xl(:, i) = var1;

t1 = strcat ('01-Oct-', num2str(HighFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(HighFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'daily', 'sum');
var1 = table2array(TT);
t1 = TT.Time;

if length(var1) == 366
    var1= [var1(1:59); var1(61:end)]
end 
 
Xh(:, i) = var1;
    end 
  Xlm(:, 4) = mean(Xl, 2);
 Xhm(:, 4) = mean(Xh, 2);
 
    x = firnmelt/24;
       clear Xl Xh
    for i = 1:length(LFyr)
t1 = strcat ('01-Oct-', num2str(LowFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(LowFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'daily', 'sum');
var1 = table2array(TT);
t1 = TT.Time;

if length(var1) == 366
    var1= [var1(1:59); var1(61:end)]
end 
Xl(:, i) = var1;

t1 = strcat ('01-Oct-', num2str(HighFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(HighFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'daily', 'sum');
var1 = table2array(TT);
t1 = TT.Time;

if length(var1) == 366
    var1= [var1(1:59); var1(61:end)]
end 
 
Xh(:, i) = var1;
    end 
  Xlm(:,5) = mean(Xl, 2);
 Xhm(:, 5) = mean(Xh, 2);
 
    x = rainfallrunoff;
       clear Xl Xh
    for i = 1:length(LFyr)
t1 = strcat ('01-Oct-', num2str(LowFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(LowFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'daily', 'sum');
var1 = table2array(TT);
t1 = TT.Time;

if length(var1) == 366
    var1= [var1(1:59); var1(61:end)]
end 
Xl(:, i) = var1;

t1 = strcat ('01-Oct-', num2str(HighFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(HighFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b, hru)./1000; %m w.e per hru
xb = xa.*hru_area(hru)*10^6;% snow per hru in m3
xc = sum(xb, 2)/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'daily', 'sum');
var1 = table2array(TT);
t1 = TT.Time;

if length(var1) == 366
    var1= [var1(1:59); var1(61:end)]
end 
 
Xh(:, i) = var1;
    end 
  Xlm(:, 6) = mean(Xl, 2);
 Xhm(:, 6) = mean(Xh, 2);

 
 x = basinflow;
       clear Xl Xh
    for i = 1:length(LFyr)
t1 = strcat ('01-Oct-', num2str(LowFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(LowFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b) %m w.e per hru
xc = xa/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'daily', 'sum');
var1 = table2array(TT);
t1 = TT.Time;

if length(var1) == 366
    var1= [var1(1:59); var1(61:end)]
end 
Xl(:, i) = var1;

t1 = strcat ('01-Oct-', num2str(HighFlow(i)-1), {' '}, '01:00'); a = find(timeCRHM==datetime(t1));
t2 = strcat ('30-Sep-', num2str(HighFlow(i)), {' '}, '01:00');b = find(timeCRHM==datetime(t2));
xa  = x(a:b) %m w.e per hru
xc = xa/(sum(hru_area(hru)*10^6)); % m3 w.e. 

T = timetable (timeCRHM(a:b), xc);
TT = retime(T, 'daily', 'sum');
var1 = table2array(TT);
t1 = TT.Time;

if length(var1) == 366
    var1= [var1(1:59); var1(61:end)]
end 
 
Xh(:, i) = var1;
    end 
    
  Xlm(:, 7) = mean(Xl, 2);
 Xhm(:, 7) = mean(Xh, 2);
 %% Figure of that
 close all
clow = [0.2 .2 .2]
chigh = [0   102   204]/255;

fig = figure('units','inches','outerposition',[0 0 8 6]);
sp1 = subplot (4,1,1)
 plot (Xlm(:, 1), 'Color', clow, 'Linewidth', 1) ;
 hold on ; plot(Xhm(:,1),'Color', chigh, 'Linewidth', 1) 
mean(Xhm(:, 1) - Xlm(:, 1))
ylabel ({'Temperature'; '(^{\circ}C)'})
xlim ([0 365])
xticks([15:30.5:365])
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1) possb(2) possb(3) possb(4)+0.03]);
lg = legend ( 'Low', 'High', 'Location', 'Northwest', 'Orientation', 'Horizontal')
pos = get(lg,'Position');
set(lg,'Position',[pos(1)+0.01 pos(2)+0.01 pos(3) pos(4)]);
xticklabels ([])
text (2, 10, 'a)')
ylim([-20 13])
grid on
box on


sp1 = subplot(4,1,2)
h = bar ([Xlprecip(:, 2), Xhprecip(:,2)]*1000, 1)
set(h(1),'facecolor',clow)
set(h(2),'facecolor',chigh)
xlim ([0.5 12.49])
ylim ([0 125])
xticklabels ([])
ylabel (' Rain (mm)')
text (0.6, 115, 'b)')
grid on
box on
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1) possb(2) possb(3) possb(4)+0.03]);
lg = legend ( 'Low', 'High', 'Location', 'Northwest', 'Orientation', 'Horizontal')
pos = get(lg,'Position');
set(lg,'Position',[pos(1)+0.01 pos(2)+0.01 pos(3) pos(4)]);
xticklabels ([])

sp1 =  subplot(4,1,3)
h = bar ([Xlprecip(:, 3), Xhprecip(:, 3)]*1000, 1)
set(h(1),'facecolor',clow)
set(h(2),'facecolor',chigh)
xlim ([0.5 12.49])
ylim ([0 450])
ylabel (' Snow (mm)')
text (.6, 405, 'c)')
grid on
box on
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1) possb(2) possb(3) possb(4)+0.03]);
lg = legend ( 'Low', 'High', 'Location', 'Northwest', 'Orientation', 'Horizontal')
pos = get(lg,'Position');
set(lg,'Position',[pos(1)+0.22 pos(2)+0.01 pos(3) pos(4)]);
xticklabels ([])

sp1 = subplot (4,1,4)
 plot (Xlm(:,2)*1000, 'Color',clow, 'Linewidth', 1) ;
 hold on ; plot(Xhm(:,2)*1000,'Color', chigh, 'Linewidth', 1) 
xlim ([0 365])
xticks([15:30.5:365])
ylim ([0 1500])
xticklabels ({'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep'})
ylabel ('SWE (mm)')
grid on
box on
text (3, 1350, 'd)')
box on
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1) possb(2) possb(3) possb(4)+0.03]);
lg = legend ( 'Low', 'High', 'Location', 'Northwest', 'Orientation', 'Horizontal')
pos = get(lg,'Position');
set(lg,'Position',[pos(1)+0.02 pos(2)+0.01 pos(3) pos(4)]);

%
tightfig(fig)
filename = 'High_lowFlow_Met'
saveas (gcf,strcat(figdir, filename), 'png')
saveas (gcf,strcat(figdir, filename), 'pdf')
savefig (gcf,strcat(figdir, filename))
%%
DifferenceTemp = Xhm(:, 1) - Xlm(:, 1)
MeanDifference = mean(DifferenceTemp)
t = datetime ('01-Oct-2000'):days(1):datetime ('30-Sep-2001')
T = timetable(t', DifferenceTemp)
TT=retime(T, 'monthly','mean')
DifferenceTemp_monthly = table2array(TT);
DifferenceTemp_monthly_time = TT.Time;

figure;
plot(DifferenceTemp_monthly_time, DifferenceTemp_monthly, '-x')

% difference in total precip
TotalRain= sum([Xlprecip(:, 2), Xhprecip(:,2)])*1000 % rain
TotalRain(2)*100/TotalRain(1)
TotalRain(2)-TotalRain(1)
TotalSnow= sum([Xlprecip(:, 3), Xhprecip(:,3)])*1000 % snow
TotalSnow(2)*100/TotalSnow(1)
TotalSnow(2)-TotalSnow(1)
%% Check statistical difference
% montlhy 
% check if normally distruibuted

%ta % significant
[h, p] = signrank(Xhprecip(:,4), Xlprecip(:, 4))

% rain % not significant
[h, p] = signrank(Xhprecip(:,2), Xlprecip(:, 2))

% snow %  significant
[h, p] = signrank(Xhprecip(:,3), Xlprecip(:, 3))

% swe % significant
[h, p] = signrank(Xhprecip(:,2), Xlprecip(:, 2))
% for the melt components

[r, p] = corrcoef(Xhm(:,1), Xlm(:,1) )
figure
scatter(Xhm(:,1), Xlm(:,1))
ref = refline (1,0)
lsline




 %% Low and Hig flow year composiiton
 % 
clow = [0.2 .2 .2]
chigh = [0   102   204]/255;
% COLOR PALETTE
cavy = [102 170 255]/255%pale blue
csnow= [0 102 204]/255%  mid blue
crain = [0 51 102]/255% dark blue
crunoff = [0 153 153]/255 % light green
cdrift = [0 102 102]/255% dark green
cET = [0 204 204]/255 % black
csnowm = [200 200 200]/255% light grey
cicem =  [140 140 140]/255% mid grey 
cfirnm = [50 50 50]/255% dark grey

fig = figure('units','inches','outerposition',[0 0 8 6]);
lw = 1
sp1=  subplot (3,2,1)
plot (smooth(Xlm(:, 3)*1000,7), 'Color', clow, 'Linewidth', lw); hold on
plot (smooth(Xhm(:, 3)*1000,7), 'Color', chigh, 'Linewidth', lw); hold on
xticks([15:30.5:365])
xticklabels ({'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep'})
xlim ([190 365])
ylim ([0 30])
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1)-0.02 possb(2) possb(3)+0.02 possb(4)]);
lg = legend ( 'Low', 'High', 'Location', 'Northwest' )
pos = get(lg,'Position');
set(lg,'Position',[pos(1)+0.01 pos(2)+0.015 pos(3) pos(4)]);
ylabel ('Snowmelt (mm d^{-1})')
grid on; box on
text(193, 27, 'a)')

sp1 =  subplot (3,2,2)
plot (smooth(Xlm(:, 4)*1000,7), 'Color', clow,  'Linewidth', lw); hold on
plot (smooth(Xhm(:, 4)*1000,7),'Color', chigh,  'Linewidth', lw); hold on
xticks([15:30.5:365])
xticklabels ({'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep'})
xlim ([190 365])
ylim ([0 20])
ylabel ('Icemelt (mm d^{-1})')
grid on; box on
text(193, 18, 'b)')
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1)-0.03 possb(2) possb(3)+0.02 possb(4)]);


sp1 =  subplot (3,2,3)
plot (smooth(Xlm(:, 5)*1000,7), 'Color', clow,  'Linewidth', lw); hold on
  plot (smooth(Xhm(:, 5)*1000,7),  'Color', chigh,  'Linewidth', lw); hold on
xticks([15:30.5:365])
xticklabels ({'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep'})
xlim ([190 365])
ylim ([0 3])
ylabel ('Firnmelt (mm d^{-1})')
grid on; box on
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1)-0.03 possb(2) possb(3)+0.02 possb(4)]);
lg = legend ( 'Low', 'High', 'Location', 'Northwest' )
pos = get(lg,'Position');
set(lg,'Position',[pos(1)+0.01 pos(2)+0.015 pos(3) pos(4)]);
grid on; box on
text(193, 2.7, 'c)')

sp1 = subplot (3,2,4)
plot (smooth(Xlm(:, 6)*1000,7), 'Color', clow,  'Linewidth', lw); hold on
  plot (smooth(Xhm(:, 6)*1000,7), 'Color', chigh,  'Linewidth', lw); hold on
xticks([15:30.5:365])
xticklabels ({'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep'})
xlim ([190 365])
ylim ([0 5.6])
ylabel ('Rainfall Runoff (mm d^{-1})')
grid on; box on
text(193, 5.0, 'd)')
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1)-0.03 possb(2) possb(3)+0.02 possb(4)]);

sp1 =  subplot (3,2,5)
h =  area([smooth(Xlm(:, 6)*1000,7) smooth(Xlm(:, 5)*1000,7) smooth(Xlm(:, 4)*1000,7) smooth(Xlm(:, 3)*1000, 7)])
set(h(1),'facecolor', crunoff)
set(h(3),'facecolor', cicem)
set(h(2),'facecolor', cfirnm)
set(h(4),'facecolor', csnowm)
xlim ([190 365])
ylim ([0 50])
xticks([15:30.5:365])
xticklabels ({'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep'})
ylabel ('Low Flow (mm d^{-1})')
grid on; box on
neworder = [4 3 2 1];
labels= {'Rainfall Runoff','Firnmelt', 'Icemelt','Snowmelt'};
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1)-0.03 possb(2) possb(3)+0.02 possb(4)]);
lg = legend(h(neworder), labels(neworder), 'Location', 'northwest', 'Color', [1 1 1], 'Orientation', 'Vertical', 'Fontsize', 8);
pos = get(lg,'Position');
set(lg,'Position',[pos(1)-0.015 pos(2)-0.015 pos(3) pos(4)]);
text(193, 45, 'e)')

sp1 = subplot (3,2,6)
h= area([smooth(Xhm(:, 6)*1000,7) smooth(Xhm(:, 5)*1000,7) smooth(Xhm(:, 4)*1000,7) smooth(Xhm(:, 3)*1000, 7)])
set(h(1),'facecolor', crunoff)
set(h(3),'facecolor', cicem)
set(h(2),'facecolor', cfirnm)
set(h(4),'facecolor', csnowm)
xlim ([190 365])
ylim ([0 50])
xticks([15:30.5:365])
xticklabels ({'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep'})
ylabel ('High Flow (mm d^{-1})')
grid on; box on
text(193, 45, 'f)')
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1)-0.03 possb(2) possb(3)+0.02 possb(4)]);
tightfig(fig)

filename = 'High_lowFlow_Component'
saveas (gcf,strcat(figdir, filename), 'png')
saveas (gcf,strcat(figdir, filename), 'pdf')
savefig (gcf,strcat(figdir, filename))

%% Gwt the number for these;
IncreaseinComponent(1, :) = [sum(Xhm(:, 3)) sum(Xlm(:, 3)) sum(Xhm(:, 3))*100/sum(Xlm(:, 3))-100];
IncreaseinComponent(2, :) = [sum(Xhm(:, 4)) sum(Xlm(:, 4)) sum(Xhm(:, 4))*100/sum(Xlm(:, 4))-100];
IncreaseinComponent(3, :) = [sum(Xhm(:, 5)) sum(Xlm(:, 5)) sum(Xhm(:, 5))*100/sum(Xlm(:, 5))-100]
IncreaseinComponent(4, :) = [sum(Xhm(:, 6)) sum(Xlm(:, 6)) sum(Xhm(:, 6))*100/sum(Xlm(:, 6))-100]
IncreaseinComponent(5, :) = [sum(Xhm(:, 7)) sum(Xlm(:, 7)) sum(Xhm(:, 7))*100/sum(Xlm(:, 7))-100]
IncreaseinComponent(6, :) = [mean(Xhm(:, 1)) mean(Xlm(:, 1)) mean(Xhm(:, 1))-mean(Xlm(:, 1))]

lab = {'snowmelt_m','icemelt_m','firnmelt_m','rainfallrunoff_m','basnflow_M','ta_C'}'
colname = {'Variable' 'HighFlow','LowFlow','Difference'}
IncreaseinComponent_HighLowFlow = table(lab, IncreaseinComponent(:, 1), IncreaseinComponent(:,2), IncreaseinComponent(:, 3))
IncreaseinComponent_HighLowFlow.Properties.VariableNames = colname;
writetable(IncreaseinComponent_HighLowFlow, strcat(savedir, 'IncreaseinComponent_HighLowFlow.txt'))

% for precip)
IncreaseinMonthlyValues(1, :) = [sum(Xhprecip(:, 1)) sum(Xlprecip(:, 1)) sum(Xhprecip(:, 1))*100/sum(Xlprecip(:, 1))-100];
IncreaseinMonthlyValues(2, :) = [sum(Xhprecip(:, 2)) sum(Xlprecip(:, 2)) sum(Xhprecip(:, 2))*100/sum(Xlprecip(:, 2))-100];
IncreaseinMonthlyValues(3, :) = [sum(Xhprecip(:, 3)) sum(Xlprecip(:, 3)) sum(Xhprecip(:, 3))*100/sum(Xlprecip(:, 3))-100];


lab = {'precip_m',  'rain_m', 'snow_m'}'
colname = {'Variable' 'HighFlow','LowFlow','Difference'}
IncreaseinPrecip_HighLowFlow = table(lab, IncreaseinMonthlyValues(:, 1), IncreaseinMonthlyValues(:,2),IncreaseinMonthlyValues(:, 3))
IncreaseinPrecip_HighLowFlow.Properties.VariableNames = colname;
writetable(IncreaseinPrecip_HighLowFlow, strcat(savedir, 'IncreaseinPrecip_HighLowFlow.txt'))
