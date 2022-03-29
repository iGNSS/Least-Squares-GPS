function [m, n] = FindInt(t, tb, Interval, m)

%find vector interval (indexes m,n) that corresponds to the fit interval 
% t - vector of epochs (time given in JD)
% tb - time of beginning
% Interval - duration of the interval in JD
% m - index of vector t that corresponds to the beginning of the interval

n = -1;
te = tb + Interval;  % time of end of interval
for i = m:length(t)
    if (t(i) - tb) <= 0
        m = i;
    end
    if (t(i) - te) >= 0
        n = i;
        return;
    end
end