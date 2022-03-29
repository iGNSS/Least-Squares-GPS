classdef EphData
properties
dtStart DateTime,dtEnd DateTime,dTransTime,dSpare1,dSpare2,dSpare3,dGpsWeek,dPDataFlag
dSVaccur,dSVhealth,dTGD,dIODC,di0,dCrc,dOmega,dOMEGADot,dIdot,dCodeOnL2,dEccent,dCus
dSqrtA,dToe,dCic,dOMEGA,dCis,iEpochSecond,dClockBias,dClockDrift,dClockDriftRate,dIDOE
dCrs,dDeltaN,dM0,dCuc,iPRN,iEpochYear,iEpochMonth,iEpochDay,iEpochHour,iEpochMinute

end

methods
    function this = EphData()

        % Broadcast Emepheris from Rinex file 
        % PRN / EPOCH / SV CLK
        % Satellite PRN number           : iPRN
        % Epoch year                     : iEpochYear
        % Epoch month                    : iEpochMonth
        % Epoch day                      : iEpochDay
        % Epoch hour                     : iEpochHour
        % Epoch minute                   : iEpochMinute
        % Epoch second                   : iEpochSecond
        % Sat Clock Bias (sec)           : dClockBias
        % Sat Clock drift(sec/sec)       : dClockDrift
        % Sat Clock drift rate (sec/sec2): dClockDriftRate
        %
        % BROADCAST ORBIT - 1
        % IDOE Issue of Data, Ephemeris  : dIDOE
        % Crs (meters)                   : dCrs
        % Delta n (radians/sec)          : dDeltaN
        % M0 (radians)                   : dM0
        %
        % BROADCAST ORBIT - 2
        % Cuc (radians)                  : dCuc
        % e Eccenricity                  : dEccent
        % Cus (radians)                  : dCus
        % sqrt(A) (sqrt(m))              : dSqrtA
        %
        % BROADCAST ORBIT - 3
        % Toe Time of Ephemeris (sec of GPS week): dToe
        % Cic (radians)                  : dCic
        % OMEGA (radians)                : dOMEGA
        % Cis (radians)                  : dCis
        %
        % BROADCAST ORBIT - 4
        % i0                             : di0
        % Crc (radians)                  : dCrc
        % omega (radians)                : dOmega
        % OMEGA Dot (radians)            : dOMEGADot
        %
        % BROADCAST ORBIT - 5
        % Idot                           : dIdot
        % Codes on L2 channel            : dCodeOnL2
        % GPS Week # (to go with TOE)    : dGpsWeek
        % L2 P data flag                 : dPDataFlag
        %
        % BROADCAST ORBIT - 6
        % SV Accuracy (meters)           : dSVaccur
        % SV health   (MSB only)         : dSVhealth
        % TGD                            : dTGD
        % IODC Issue of Data, Clock      : dIODC
        %
        % BROADCAST ORBIT - 7
        % Transmission time of message   : dTransTime
        % Spare1                         : dSpare1
        % Spare2                         : dSpare2
        % Spare3                         : dSpare3
        
        ep = DateTime;
        this.iPRN=0;
        this.iEpochYear=0;
        this.iEpochMonth=0;
        this.iEpochDay=0;
        this.iEpochHour=0;
        this.iEpochMinute=0;
                     
                   this.iEpochSecond=0;
                   this.dClockBias=0;
                   this.dClockDrift=0;
                   this.dClockDriftRate=0;
                   this.dIDOE=0;
                   this.dCrs=0;
                   this.dDeltaN=0;
                   this.dM0=0;
                   this.dCuc=0;
                   this.dEccent=0;
                   this.dCus=0;
                   this.dSqrtA=0;
                   this.dToe=0;
                   this.dCic=0;
                   this.dOMEGA=0;
                   this.dCis=0;
                   this.di0=0;
                   this.dCrc=0;
                   this.dOmega=0;
                   this.dOMEGADot=0;
                   this.dIdot=0;
                   this.dCodeOnL2=0;
                   this.dGpsWeek=0;
                   this.dPDataFlag=0;
                   this.dSVaccur=0;
                   this.dSVhealth=0;
                   this.dTGD=0;
                   this.dIODC=0;
                   this.dGpsWeek=0;
                   this.dPDataFlag=0;
                   this.dSVaccur=0;
                   this.dSVhealth=0;
                   this.dTGD=0;
                   this.dIODC=0;
                   this.dTransTime=0;
                   this.dSpare1=0;
                   this.dSpare2=0;
                   this.dSpare3=0;
                   this.dtStart=ep;
                   this.dtEnd=ep;
           
        

    end
        function a = get(b,par)

        a=b.(par);
        end

        function b = set(b,par,in)
        b.(par)=in;
        end
        
        function display(b)
        
        % PRN / EPOCH / SV CLK
        fprintf('Satellite PRN number               iPRN            : %4i\n',b.iPRN)
        fprintf('Epoch year                         iEpochYear      : %4i\n',b.iEpochYear)
        fprintf('Epoch month                        iEpochMonth     : %4i\n',b.iEpochMonth)
        fprintf('Epoch day                          iEpochDay       : %4i\n',b.iEpochDay)
        fprintf('Epoch hour                         iEpochHour      : %4i\n',b.iEpochHour)
        fprintf('Epoch minute                       iEpochMinute    : %4i\n',b.iEpochMinute)
        fprintf('Epoch second                       iEpochSecond    : %4i\n',b.iEpochSecond)
        fprintf('Sat Clock Bias (sec)               dClockBias      : %6.8f\n',b.dClockBias)
        fprintf('Sat Clock drift(sec/sec)           dClockDrift     : %6.8f\n',b.dClockDrift)
        fprintf('Sat Clock drift rate (sec/sec2)    dClockDriftRate : %6.8f\n',b.dClockDriftRate)
        
        % BROADCAST ORBIT - 1
        fprintf('IDOE Issue of Data, Ephemeris      dIDOE   : %6.8f\n',b.dIDOE)
        fprintf('Crs (meters)                       dCrs    : %6.8f\n',b.dCrs)
        fprintf('Sat Clock drift rate (sec/sec2)    dM0     : %6.8f\n',b.dDeltaN)
        fprintf('M0 (radians)                       dDeltaN : %6.8f\n',b.dM0)
        
        % BROADCAST ORBIT - 2
        fprintf('Cuc (radians)                      dCuc    : %6.8f\n',b.dCuc)
        fprintf('e Eccenricity                      dEccent : %6.8f\n',b.dEccent)
        fprintf('Cus (radians)                      dCus    : %6.8f\n',b.dCus)
        fprintf('sqrt(A) (sqrt(m)                   dSqrtA  : %6.8f\n',b.dSqrtA)
        
        % BROADCAST ORBIT - 3
        fprintf('Toe Time of Ephemeris(sec of Week) dToe    : %6.8f\n',b.dToe)
        fprintf('Cic (radians)                      dCic    : %6.8f\n',b.dCic)
        fprintf('OMEGA (radians)                    dOMEGA  : %6.8f\n',b.dOMEGA)
        fprintf('Cis (radians)                      dCis    : %6.8f\n',b.dCis)
        
        % BROADCAST ORBIT - 4
        fprintf('i0                                 di0      : %6.8f\n',b.di0)
        fprintf('Crc (radians)                      dCrc     : %6.8f\n',b.dCrc)
        fprintf('omega (radians)                    dOmega   : %6.8f\n',b.dOmega)
        fprintf('OMEGA Dot (radians)                dOMEGADot: %6.8f\n',b.dOMEGADot)
        
        % BROADCAST ORBIT - 5
        fprintf('Idot                               dIdot       : %6.8f\n',b.dIdot)
        fprintf('Codes on L2 channel                dCodeOnL2   : %6.8f\n',b.dCodeOnL2)
        fprintf('GPS Week # (to go with TOE)        dGpsWeek    : %6.8f\n',b.dGpsWeek)
        fprintf('L2 P data flag                     dPDataFlag  : %6.8f\n',b.dPDataFlag)
        
        % BROADCAST ORBIT - 6   
        fprintf('SV Accuracy (meters)               dSVaccur    : %6.8f\n',b.dSVaccur)
        fprintf('SV health   (MSB only)             dSVhealth   : %6.8f\n',b.dSVhealth)
        fprintf('TGD (seconds)                      dTGD        : %6.8f\n',b.dTGD)
        fprintf('IODC Issue of Data, Clock          dIODC       : %6.8f\n',b.dIODC)
        
        % BROADCAST ORBIT - 7
        fprintf('Transmission time of message       dTransTime  : %6.8f\n',b.dTransTime)
        fprintf('Spare 1                            dSpare1     : %6.8f\n',b.dSpare1)
        fprintf('Spare 2                            dSpare2     : %6.8f\n',b.dSpare2)
        fprintf('Spare 3                            dSpare3     : %6.8f\n',b.dSpare3)
        
        b.dtStart
        b.dtEnd
        
        
        end

end

end
