function [VarBasinTot,VarMperHRU, TimeM] = MonthlyHruWeigthedBasinMean(Time, Var, hru, hru_area)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
T = timetable (Time,Var(:, hru));
TT = retime(T, 'monthly', 'mean');
VarM = table2array(TT);
TimeM= TT.Time;

for i = 1:length(hru)
VarMperHRU(:, i) = VarM(:, i); % m w.e.
VarMperHRUmulArea(:, i) = VarMperHRU(:, i)*[(hru_area(i)*10^6)/sum(hru_area(hru)*10^6)];% m3 w e.
end 
VarBasinTot = mean(VarMperHRUmulArea, 2);
end

