function  en = FindEn(vt, t)

%FindEn finds interval number of epoch t in the vector vt
% vt - vector of epochs (time given in [DateTime])
% t - given epoch in [DateTime]
% en - index of vt, where t falls 

%Written by Milan Horemuz, last modified 2005-02-03


en = -1;
lvt = length(vt);
if lvt < 1
    return;
end
if lvt == 1
    en = 1;
    return;
end
interv = vt(2) - vt(1);
last = length(vt);
vt(last+1) = vt(last) + interv;
for i = 2:(last+1)
    if  (vt(i) - t) >= 0
        en = i-1;
        return;
    end
end