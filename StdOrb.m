classdef StdOrb
properties (SetAccess=public)
   
          dTbeg DateTime
          dTend DateTime
          dInterval,iPolyOrder,dXCoef,dXmu,dYCoef,dYmu,dZCoef,dZmu,dFitInt
          dTimeSet DateTime
          iPRN
          dtClTime DateTime
          dClErr,dTGD
                 
end

    methods
        function this = StdOrb() %Constructor
        %StdOrb - Standard orbits in form of polynoms
        %dTbeg -  time of beginning [DateTime]
        %dTend -  time of end [DateTime]
        %dInterval - time interval of one set of coefficients [seconds]
        %iPolyOrder - Number of coefficients in polynoms
        %dXCoef, dYCoef, dZCoef  - matrix of coefficients for X, Y, Z coordinate, size (s,m,n).s-number of satellites, m - number
        %of intevals, n - number of coefficients
        %FitInt - interval used to fit the data; FitInt >= dInterval [seconds]
        %dTimeSet - vector containing time of start of interval
        %iPRN - satellite numbers; if iPRN == 0, no data for the satellite
        %dtClTime - vector of DateTime, containing time of clock corrections stored in vector dClErr
        %dClErr - vector of satellite clock corrections; clock correction for a
        %given epoch is linearly interpolated
        
        %NofCoef = 12;
        iPolyOrder = 16;  %
        sz = iPolyOrder + 1;  %number of polynomial coefficients
        dtTime = DateTime;
        this.dTbeg= dtTime; this.dTend = dtTime; this.dInterval=2*3600;
         this.iPolyOrder=iPolyOrder;this.dXCoef = zeros(1,1,sz);this.dXmu=zeros(1,1,2);
         this.dYCoef= zeros(1,1,sz);this.dYmu=zeros(1,1,2); this.dZCoef= zeros(1,1,sz);
         this.dZmu=zeros(1,1,2);this.dFitInt=4*3600;this.dTimeSet=dtTime;
         this.iPRN=zeros(1,1);this.dtClTime=dtTime;this.dClErr=0;this.dTGD=0;
         this.dtClTime=DateTime;this.dTimeSet=DateTime;this.dTend=DateTime;this.dTbeg=DateTime;
        end
        function so = CompStdOrb(so, pe, stTgd)
        %CompStdOrb computes polynomial coefficients for all satellites available in precise
        %ephemeris peEph
        %so = CompStdOrb(so, pe)
        %so - Standard orbit (class StdOrb)
        %pe - Coordinate ephemeris (class CoordEph)
        %stTgd = struct('vDateTime', dtTime, 'cTGD',0) - this structure is created by getTgd() function
        %Coordinates for an epoch are computed by method GetSatCoord
        x = get(pe,'dX'); %extract x-coordinate for all satellites
        y = get(pe,'dY'); 
        z = get(pe,'dZ'); 
        t = get(pe,'dTime');  %epochs 
        cl = get(pe,'dDts');   %clock corrections
        prn = get(pe,'iPRN'); %satellie numbers
        ne = length(t);     %number of epochs
        tb = t(1);          %first epoch, for which the coordinates are given
        te = t(ne);  %last epoch
        tm = (so.dFitInt - so.dInterval)/2; %lead-in and lead-out interval
        Nint = floor((te - tb - 2*tm)/so.dInterval); %number of fit intervals
        m = 1;
        for i = 1:Nint
            so.dTimeSet(i) = tb + tm + (i-1)*so.dInterval; %time of start of interval
            tg = FindTGD(stTgd, so.dTimeSet(i));  %find value of Tgd closest to time so.dTimeSet(i) for all satellites
            [m, n] = FindInt(t, so.dTimeSet(i)-tm, so.dFitInt, m); %find vector interval (indexes m,n) that corresponds to the fit interval 
            if n < 0     %it was not possible to find data for current interval
                so.dTimeSet(i) = [];
                break;
            end
            tt = (t(m:n) - so.dTimeSet(i));  %reduce time to get better numerical stability
            for j = 1:size(x,1) %loop over satellites
                if prn(j, n) < 1  %no data for satellite j
                    coef = zeros(1,size(so.dXCoef,3));
                    so.dXCoef(j,i,:) =  coef;   %placeholder for satellite
                    so.dYCoef(j,i,:) =  coef;
                    so.dZCoef(j,i,:) =  coef;
                    continue;
                end
                [coef, S, mu] = polyfit(tt, x(j,m:n), so.iPolyOrder);
                %yy=polyval(coef,tt,[],mu);
                %dd= yy-x(j,m:n);
                %plot(dd);
                so.dXCoef(j,i,:) =  coef;
                so.dXmu(j,i,:) =  mu';
                [coef, S, mu] = polyfit(tt, y(j,m:n), so.iPolyOrder);
                so.dYmu(j,i,:) =  mu';
                so.dYCoef(j,i,:) =  coef;
                [coef, S, mu] = polyfit(tt, z(j,m:n), so.iPolyOrder);
                so.dZCoef(j,i,:) =  coef;
                so.dZmu(j,i,:) =  mu';
                so.dTGD(j,i) = tg(j);  
            end
        end
        %interval, for which the std. orbit is valid
        
        so.dTbeg = so.dTimeSet(1); %begin
        so.dTend = so.dTimeSet(length(so.dTimeSet)) + so.dInterval;  %end
        so.iPRN = prn;
        so.dtClTime = t;
        so.dClErr = cl;
        end
        function a = get(b,par)

            a=b.(par);
        end
        function display(b)
        
        fprintf('Standard orbit interval:\n');
        fprintf('Beginning:  %d %d %d %d %d %.1f\n', get(b.dTbeg,'year'), get(b.dTbeg,'month'), get(b.dTbeg,'day'), get(b.dTbeg,'hour'), get(b.dTbeg,'min'), get(b.dTbeg,'sec'));
        fprintf('End:  %d %d %d %d %d %.1f\n', get(b.dTend,'year'), get(b.dTend,'month'), get(b.dTend,'day'), get(b.dTend,'hour'), get(b.dTend,'min'), get(b.dTend,'sec'));
        fprintf('Polynomial order:  %d\n', b.iPolyOrder);
        fprintf('Fit interval:  %.1f hours\n', b.dFitInt/3600);
        fprintf('Validity interval:  %.1f hours\n', b.dInterval/3600);
        end
        function ret = GetSatCoord(std, prn, t)
        
        %GetSatCoord - computes satellite's coordinates and clock correction for a
        %given epoch, using standard ephemeris
        %ret = GetSatCoord(std, prn, t)
        %std - standard orbit in polynom form (class StdOrb) created by CompStdOrb
        %prn - satellite number
        %t - eopch in [DateTime]
        
        if (t - std.dTbeg) < 0  | (t - std.dTend) > 0
            %if DateTime.difference(t,std.dTbeg) | DateTime.difference(t,std.dTend)
            
            smaller=t-std.dTbeg
            larger=t-std.dTend
            result_date = std.dTend-std.dTbeg
        
            time_not_included=t
            
            error('Given epoch is outside of the standard ephemeris');
        end
        en = FindEn(std.dTimeSet, t);  %Finds epoch number en
        if en < 1
            t
            error('Should not be here: something is wrong in findEn');
        end
        ent = FindEn(std.dtClTime, t);
        if ent < 1
            t
            error('Should not be here: something is wrong in findEn finding dtClTime');
        end
        
        ret = struct('x',0,'y',0,'z',0,'dt',0, 'Tgd',0);
        ret.Tgd = std.dTGD(prn, en);
        coefX = zeros(1,size(std.dXCoef,3));
        coefY = coefX;
        coefZ = coefX;
        coefDT =coefX;
        Xmu = zeros(1,2);
        Ymu = Xmu;
        Zmu = Xmu;
        %Tmu = Xmu;
        coefX(:) = std.dXCoef(prn, en, :);
        if abs(coefX(length(coefX))) < 1e-12
            prn
            error('No data for given prn');
        end
        coefY(:) = std.dYCoef(prn, en, :);
        coefZ(:) = std.dZCoef(prn, en, :);
        %coefDT(:) = std.dDtsCoef(prn, en, :);
        Xmu(:) = std.dXmu(prn, en, :);
        Ymu(:) = std.dYmu(prn, en, :);
        Zmu(:) = std.dZmu(prn, en, :);
        %Tmu(:) = std.dTmu(prn, en, :);
        tt = t - std.dTimeSet(en);
        ret.x = polyval(coefX, tt, [], Xmu);
        ret.y = polyval(coefY, tt, [], Ymu);
        ret.z = polyval(coefZ, tt, [], Zmu);
        %ret.dt = polyval(coefDT, tt, [], Tmu);
        %coefVX = polyder(coefX);
        %tv = (tt - Xmu(1))/Xmu(2);
        %ret.vx = polyval(coefVX, tv);
        
        coefVX = coefX;
        coefVY = coefY;
        coefVZ = coefZ;
        n = length(coefX);
        coefVX(n) = []; %delete the last coefficient (first derivative = 0)
        coefVY(n) = [];
        coefVZ(n) = [];
        n = n-1; %length of coefVX
        for i=1:n
            coefVX(i) = coefVX(i)*(n-i+1);
            coefVY(i) = coefVY(i)*(n-i+1);
            coefVZ(i) = coefVZ(i)*(n-i+1);
        end
        tv = (tt - Xmu(1))/Xmu(2);
        ret.vx = polyval(coefVX, tv)/Xmu(2);
        tv = (tt - Ymu(1))/Ymu(2);
        ret.vy = polyval(coefVY, tv)/Ymu(2);
        tv = (tt - Zmu(1))/Zmu(2);
        ret.vz = polyval(coefVZ, tv)/Zmu(2);
        
        targ = std.dtClTime(ent:ent+1) - std.dtClTime(ent);
        [coefDT, S, mu] = polyfit(targ, std.dClErr(prn, ent:ent+1), 1);
        tt = t - std.dtClTime(ent);
        ret.dt = polyval(coefDT, tt, [], mu);
        end
        function b = set(b,par,in)
        b.(par)=in;
        end

    end

end
