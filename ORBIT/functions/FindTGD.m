function Tgd = FindTGD(stTgd, ep)

%Finds value of group delay (Tgd) closest to epoch ep

ind = FindEn(stTgd.vDateTime, ep);
if ind < 0
    ep
    error('Could not find any Tgd for given epoch');
end
Tgd = stTgd.cTGD(:,ind);  %vector of Tgd for all matrixes
