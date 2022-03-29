function [x,y,z,ds,TGD]=brEphCoordCalcOneDataSet(vObs,cObsTime)

% This function calculates coordinates of one satellite at observation time
% "cObsTime" from observations in the observation vector "vObs". 
% 
% To call the function, use 
%
%    [x,y,z,ds,TGD]=brEphCoordCalcOneDataSet(vObs,cObsTime)
%
% where 
%   (x,y,z)  - are satellite coordinates in WGS84
%   ds       - delta SV PRN code phase time offset in seconds
%   TGD      - Differential Group delay in seconds
%   vObs     - observation vector
%   cObsTime - The current observation time in reciver time
%


% Constatants *******************************************************

    my=3.986005e14;             % Gravitional constant for WGS84

    OmegaDotE=7.2921151467e-5;  % Earth's rotation rate for WGS84  
    
    C=299792458;                % The speed of light in vacuum
    
    F=-2*sqrt(my)/C^2;
   
   
% Data conversion ***************************************************
    
% PRN/EPOCH/SC CLK

    SatNr=get(vObs,'iPRN');

    EpochYear =get(vObs,'iEpochYear');

    EpochMonth=get(vObs,'iEpochMonth');
    
    EpochDay=get(vObs,'iEpochDay');
    
    EpochHour=get(vObs,'iEpochHour');
    
    EpochMinute=get(vObs,'iEpochMinute');
    
    EpochSecond=get(vObs,'iEpochSecond');
    dtTocl = DateTime(EpochYear, EpochMonth, EpochDay, EpochHour, EpochMinute, EpochSecond);
    
    a0=get(vObs,'dClockBias');
    
    a1=get(vObs,'dClockDrift');
    
    a2=get(vObs,'dClockDriftRate');

    
% BroadcastOrbit - 1

    IODE = get(vObs,'dIDOE');        % Issue of data (ephemeris)
    
    Crs = get(vObs,'dCrs');          % Amplitude of second-order harmonic pertubations
    
    Delta_n = get(vObs,'dDeltaN');   % Mean motion difference from computed value
    
    M0 = get(vObs,'dM0');            % Mean anomaly at reference time
    
    
% BroadcastOrbit - 2

    Cuc=get(vObs,'dCuc');           % Amplitude of second-order harmonic pertubations
    
    e=get(vObs,'dEccent');          % Eccentricity
    
    Cus=get(vObs,'dCus');           % Amplitude of second-order harmonic pertubations
    
    a=(get(vObs,'dSqrtA'))^2;       % (Square root of the semi major axis)^2
    
    
% BroadcastOrbit - 3

    Toe=get(vObs,'dToe');           % Ephemeris reference time
    
    Cic=get(vObs,'dCic');           % Amplitude of second-order harmonic pertubations
    
    OMEGA=get(vObs,'dOMEGA');       % Longitude of ascending node of orbit plane at beginning of week
    
    Cis=get(vObs,'dCis');           % Amplitude of second-order harmonic pertubations
    
    
% BroadcastOrbit - 4

    i0=get(vObs,'di0');             % Inclination angle at reference time
    
    Crc=get(vObs,'dCrc');           % Amplitude of second-order harmonic pertubations
    
    omega=get(vObs,'dOmega');       % Argument of perigee
    
    OMEGA_DOT=get(vObs,'dOMEGADot');% Rate of right ascension
    
    
% BroadcastOrbit - 5

    IDOT=get(vObs,'dIdot');                     % Rate of incliniation angle
    
    Codes_On_L2_channel=get(vObs,'dCodeOnL2');
    
    GPS_WEEK=get(vObs,'dGpsWeek');              % To go with TOE
    dtToe = DateTime(GPS_WEEK, Toe);  %convert to DateTime
    
    L2_data_flag=get(vObs,'dPDataFlag'); 

    
% BroadcastOrbit - 6

    SV_accuracy=get(vObs,'dSVaccur');       % (meters)
    
    SV_health=get(vObs,'dSVhealth');        % (MSB only)
    
    TGD=get(vObs,'dTGD');                   % (seconds)
    
    IODC_Issue_of_data=get(vObs,'dIODC'); ; % (Clock)
    
    
% BroadcastOrbit - 7
    
    tr_time_of_message=get(vObs,'dTransTime'); % Transmission time of message (sec of GPS week))
                                               % derived e. g. from Z-count in Hand Over Word
                                               % (HOW)
    
    n0=sqrt(my/a^3); % Computed mean motion - rad/sec

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SV PRN code phase time offset
    
      diffTcl = cObsTime - dtTocl;
    
      delta_ts = a0 + a1*diffTcl + a2*diffTcl^2; % Delta SV PRN code phase time offset in seconds
                                                 % determined the first time without relativistic 
                                                 % effects.
      
      
    % Time from ephemeris reference epoch
    
      tk = cObsTime - dtToe - delta_ts;
      
      
    % Corrected mean motion
      
      n=n0+Delta_n; 
      
        
    % Mean anaomaly
    
      Mk=M0+n*tk;
      
    
    % Itteration to determine Kepler's equation for eccentric anomaly
    % initial values
        
      Ek=Mk; 
        
      diff=1;
        
    % While loop
      
      while abs(diff)>1.0e-13

            diff=Mk-Ek+e*sin(Ek);
       
            Ek=Ek+diff;
      end

%%%% Run a second time to delta_ts where the relativistic effects are
%%%% applied

   % Relativistic effects
   
      dTr = F * e * sqrt(a) * sin(Ek);
    
   
   % delta_ts with relativistic effects
   
      delta_ts=a0+a1*diffTcl + a2*diffTcl^2 + dTr; % Delta SV PRN code phase time offset in seconds
                                                   % determined the first time without relativistic 
                                                   % effects.
      
   % Time from ephemeris reference epoch
    
     tk = cObsTime - dtToe - delta_ts;

      
    % Itteration to determine Kepler's equation for eccentric anomaly
    % initial values

      diff=1;

    % While loop
      
      %while abs(diff)>1.0e-13

       %     diff=Mk-Ek+e*sin(Ek);
       
        %    Ek=Ek+diff;
      %end
   
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % True anomaly
    
      Cvk=(cos(Ek)-e)/(1-e*cos(Ek));
    
      Svk=(sqrt(1-e^2)*sin(Ek))/(1-e*cos(Ek));
      
      fk=atan2(Svk,Cvk);
      
        if fk<0
        
            fk=fk+2*pi;
            
        end

    % Argument of latituide
        fi_k=fk+omega;
    
    % Second harmonic perturbations
    
        ra_uk=Cus*sin(2*fi_k)+Cuc*cos(2*fi_k);
        
        ra_rk=Crc*cos(2*fi_k)+Crs*sin(2*fi_k);
        
        ra_ik=Cic*cos(2*fi_k)+Cis*sin(2*fi_k);
    
        
    % Corrected argument of latitude
    
        uk=fi_k+ra_uk;
    
        
    % Corrected radius
    
        rk=a*(1-e*cos(Ek))+ra_rk;
     
        
    % Corrected inclination
    
        ik=i0+ra_ik+IDOT*tk;
    
        
    % Position in orbital plane
    
        xk=rk*cos(uk);
        
        yk=rk*sin(uk);

        
    % Corrected longitude of ascending node
    
        OMEGA_k=OMEGA+(OMEGA_DOT-OmegaDotE)*tk-OmegaDotE*Toe;
    
        
    % Earth fixed coordinates
  

        x=xk*cos(OMEGA_k)-yk*cos(ik)*sin(OMEGA_k);
        
        y=xk*sin(OMEGA_k)+yk*cos(ik)*cos(OMEGA_k);
        
        z=yk*sin(ik);
        
        ds = delta_ts;              