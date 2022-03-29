classdef DateTime 
%DateTime -  class that stores time and performes conversion between various format
%A = DateTime() - fills in zeros to all variables
%A - class DateTime
%A = DateTime(y, mo, d, ho, min, sec)
%A = DateTime([y mo d ho min sec])
%y, mo, d, ho, min, sec - date and time
%A = DateTime(gw, ws) 
%gs - GPS week
%ws - seconds of GPS week 
%A = DateTime(MJD) 
%MJD - modified julian date
%dweek - day of GPS week

%Written by Milan Horemuz, last modified 2004-11-01
properties (SetAccess = public)
    MJD
    gweek
    dweek
    wsec
    year
    month
    day
    hour
    min
    sec
    DOY
end
methods
    function this = DateTime(varargin)
        
        nargin=length(varargin);
        if nargin==0
            this.MJD=0;this.gweek=0;this.dweek=0;this.wsec=0;this.year=0;this.month=0;
            this.day=0;this.hour=0;this.min=0;this.sec=0;this.DOY=0;
        end
        if nargin==6
            [y,mo,d,ho,min,sec]=deal(varargin{:});
            this.MJD=0;this.gweek=0;this.dweek=0;this.wsec=0;this.year=y;this.month=mo;
            this.day=d;this.hour=ho;this.min=min;this.sec=sec;this.DOY=0;
            this=date2GPS(this);
        end
        if isvector(varargin) && nargin==1
            varargin=cell2mat(varargin);
            y=varargin(1);mo=varargin(2);d=varargin(3);ho=varargin(4);min=varargin(5);sec=varargin(6);
            this.MJD=0;this.gweek=0;this.dweek=0;this.wsec=0;this.year=y;this.month=mo;
            this.day=d;this.hour=ho;this.min=min;this.sec=sec;this.DOY=0;
            this=date2GPS(this);
        end
        if nargin==2
            [gw,ws]=deal(varargin{:});
            this.MJD=0;this.gweek=gw;this.dweek=0;this.wsec=ws;this.year=0;this.month=0;
            this.day=0;this.hour=0;this.min=0;this.sec=0;this.DOY=0;
            this=GPS2date(this);
        end
        if nargin==1 && not(isvector(varargin))
            [mjd]=deal(varargin{:});
            this.MJD=mjd;this.gweek=0;this.dweek=0;this.wsec=0;this.year=0;this.month=0;
            this.day=0;this.hour=0;this.min=0;this.sec=0;this.DOY=0;
        end
        if nargin>6
            in=varargin{1};
            y = in(1);
            mo = in(2);
            d = in(3);
            ho = in(4);
            min = in(5);
            sec = in(6);
            this.MJD=0;this.gweek=0;this.dweek=0;this.wsec=0;this.year=y;this.month=mo;
            this.day=d;this.hour=ho;this.min=min;this.sec=sec;this.DOY=0;
            this=date2GPS(this);
        end
    end
    function  dtA = date2GPS(dtA)
        %date2GPS converts datum to GPS time and MJD
        %dtA = date2GPS(dtA)
        %dtA - object of type DateTime
        %Written by Milan Horemuz, last modified 2004-11-01
           y = dtA.year; 
           mo = dtA.month;
           if mo <= 2 
              y = y - 1; 
              mo = mo + 12;
          end
           a = 365.25*y;
           b = (mo+1)*30.6001;
           dh = dtA.hour + dtA.min/60 + dtA.sec/3600;  %hours in day
           jd = floor(a) + floor(b) + dtA.day + 1720981.5;  %+ dh/24  
           dtA.MJD = jd-2400000.5 + dh/24; 
           a = (jd - 2444244.5)/7;
           dtA.gweek = floor(a);
           wsec = (a - dtA.gweek)*7.*86400.;         % seconds of the week - not sufficient precision
           dtA.dweek = round(wsec/86400.);
           dtA.wsec = dtA.dweek*86400 + dh*3600;     % seconds of the week -  sufficient precision
    end
    function A = GPS2date(A)

        %GPS2date converts datum to GPS time and MJD
        %A = GPS2date(A)
        %A - object of type DateTime

        %Written by Milan Horemuz, last modified 2005-02-09


           jd = A.gweek*7 + A.wsec/86400 + 2444244.5;
           A.MJD = jd-2400000.5;
           a = floor(jd+0.5); 
           b = a + 1537; 
           c = floor((b-122.1)/365.25);
           d = floor(365.25*c); 
           e = floor((b-d)/30.6001);
           f = jd+0.5;
           A.day = b - d - floor(30.6001*e); % + (int) modf(jd+0.5,&pom);
           A.month = e-1-12* floor(e/14);
           A.year = c - 4715 - floor((7+ A.month)/10);
           A.dweek = floor(A.wsec/86400.);
           pom = A.wsec/3600 - A.dweek*24;
           A.hour = floor(pom);
           pom = (pom - A.hour)*60;
           %A.min = floor(pom);
        %   A.sec = (pom - A.min)*60;
           A.sec = A.wsec - A.dweek*86400 - A.hour*3600;  % - A.min*60;
           A.min = floor(A.sec/60);
           A.sec = A.sec - A.min*60;
    end
    function [doy, ep] = dayofyear(ep)

        %computes day of year
        %Written by Milan Horemuz 2005-03-04

        %y = get(ep,'year');
        ep0 = DateTime(ep.year,1,1,0,0,0);
        %doy = floor(get(ep,'MJD') - get(ep0,'MJD'))+1;
        ep.DOY = floor(ep.MJD - ep0.MJD) + 1;
        doy = ep.DOY;
    end
    function display(A)

        fprintf('Year                     : %4d\n', A.year);
        fprintf('Month                    : %4d\n', A.month);
        fprintf('Day                      : %4d\n', A.day);
        fprintf('Hour                     : %4d\n', A.hour);
        fprintf('Minute                   : %4d\n', A.min);
        fprintf('Seconds                  : %.5f\n', A.sec);
        fprintf('GPS week                 : %4d\n', A.gweek);
        fprintf('Seconds of GPS week      : %.5f\n', A.wsec);
        fprintf('Day of week              : %4d\n', A.dweek);
        fprintf('MJD                      : %.5f\n', A.MJD);
        fprintf('Day of Year              : %.5f\n', A.DOY);
        %struct('Year', A.year, 'Month', A.month, 'Day', A.day, 'Hour', A.hour, 'Minute', A.minute, 'Sec', A.sec, 'GPSweek', A.gweek, 'GPSsec', A.wsec);
    end
    function a = get(b,par)
        a=b.(par);
    end
    function r = minus(a,b)
        %Computes differences between two DateTime a-b
        %or DateTime a minus b seconds
        %Written by Milan Horemuz, last modified 2004-11-01
        if isa(a, 'DateTime') & isa(b, 'DateTime')
            [k,l] = size(a);
            [m,n] = size(b);
            if (k==m & l==n)
                dw = a.gweek - b.gweek;
                ds = a.wsec - b.wsec;
                r = dw*86400*7 + ds;
            elseif (k==1 & l>1) & (m ==1 & n==1)
                for j=1:l
                    dw(j) = a(j).gweek - b.gweek;
                    ds(j) = a(j).wsec - b.wsec;
                    r(j) = dw(j)*86400*7 + ds(j);
                end
            end


        elseif isa(a, 'DateTime') & isa(b, 'double')
            pom = a.wsec - b;
            r = DateTime(a.gweek, pom);
        else
            err = sprintf('Operator minus in DateTime does not allow arguments %s %s', class(a), class(b));
            error(err);
        end
    end
    function A = MJD2date(A)

        %MJD2date converts MJD to GPS time
        %A = MJD2date(A)
        %A is DateTime object

        %Written by Milan Horemuz, last modified 2004-11-01

        jd = A.MJD + 2400000.5;
        week = (jd - 2444244.5)/7;
        A.gweek = floor(week);
        A.wsec = (week - A.gweek)*86400*7;
        A = GPS2Date(A);
    end
    function r = plus(a,b)
        if isa(a, 'DateTime') & isa(b, 'double')
            pom = a.wsec + b;
            r = DateTime(a.gweek, pom);
        elseif isa(b, 'DateTime') & isa(a, 'double')
            pom = a + b.wsec;
            r = DateTime(b.gweek, pom);
        else
            err = sprintf('Operator plus in DateTime does not allow arguments %s %s', class(a), class(b));
            error(err);
        end
    end
    function  b= set(b,par,val)
        b.(par) = val;
    end
    


    

    
end
methods (Static)
    function r = difference(a,b)
        %Computes differences between two DateTime a-b
        %or DateTime a minus b seconds
        %Written by Milan Horemuz, last modified 2004-11-01
        if isa(a, 'DateTime') & isa(b, 'DateTime')
            [k,l] = size(a);
            [m,n] = size(b);
            if (k==m & l==n)
                dw = a.gweek - b.gweek
                ds = a.wsec - b.wsec
                r = dw*86400*7 + ds
            elseif (k==1 & l>1) & (m ==1 & n==1)
                for j=1:l
                    dw(j) = a(j).gweek - b.gweek;
                    ds(j) = a(j).wsec - b.wsec;
                    r(j) = dw(j)*86400*7 + ds(j);
                end
            end


        elseif isa(a, 'DateTime') & isa(b, 'double')
            pom = a.wsec - b;
            r = DateTime(a.gweek, pom);
        else
            err = sprintf('Operator minus in DateTime does not allow arguments %s %s', class(a), class(b));
            error(err);
        end
    end
end

end
