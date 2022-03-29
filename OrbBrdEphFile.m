classdef OrbBrdEphFile < handle
    
    properties (Access = 'private')
        brdEphData
        current_t = 0
        settings
    end
        
    methods
        
        function this=OrbBrdEphFile(settings)
        % FUNCTION
        %   Function that reads GPS broadcast ephemerides file
        % INPUTS
        %   filePath: path to file including GPS broadcast ephemerides file
        %----------------------------------------------------------------    
        % Add path to auxiliary toolbox

          this.settings=settings;
         
        end
        
        function [GPS_pos,GPSclc_corr,GPS_Tgd]=getGpsSatParam(this,t,x,sv_no,recClcBias)
        %-----------------------------------------------------------------                                           
        % FUNCTION
        %   Computes the GPS position and GPS satellite clock error at
        %   signal transmission time via  light time iteration
        % INPUTS
        %   epoch : epoch number
        %   t     : receiver date and time, (GPS time)
        %   x     : receiver position vector
        %   sv_no : satellite vehicle number
        %   recClcBias : receiver clock bias
        % OUTPUTS
        %   GPS_pos     : GPS positions for given satellites, [m]
        %   GPSclc_corr : GPS clock cprrections for given satellites, [m]
        %   GPS_Tgd     : Intrumental delays for given satellites, [m]
        %-----------------------------------------------------------------
       
        % Set Constants
          Omegae = 7.292115147e-5;  % Earth rotation rate (rad/sec)
          c=299792458;              % speed of light 
        % Check the file
          epo = DateTime(t) - recClcBias/c;
          if isnumeric(this.current_t)
              this.loadObsFile(epo)
          elseif (floor(get(this.current_t,'MJD')) ~= floor(get(epo,'MJD')))             
            this.loadObsFile(epo)
          end          
        % Compute GPS satellite coordinates using light time iteration
          
          tau0=0.072*ones(length(sv_no),1);
          
          for i=1:length(sv_no)
             d_p=1;
             while d_p>0.0001
                % Approximate transmition time
                  epo_trans=minus(epo,tau0(i));  
                % Get GPS satellite coordinates
                  coordsBr = GetSatCoord(this.brdEphData,sv_no(i),epo_trans);
                  GPS_pos(1,i)=coordsBr.x;
                  GPS_pos(2,i)=coordsBr.y;
                  GPS_pos(3,i)=coordsBr.z;
                  GPSclc_corr(i,1)=coordsBr.dt;
                  GPS_Tgd(i,1)=coordsBr.Tgd;
                % Earth rotation correction
                  omegae_dot_dt = Omegae*tau0(i);  
                  R3 = [  cos(omegae_dot_dt)   sin(omegae_dot_dt) 0;
                          -sin(omegae_dot_dt)  cos(omegae_dot_dt) 0;
                          0                    0                  1];
                  GPS_pos(:,i)=R3*GPS_pos(:,i);
                % Compute Geometric range
                  p=sqrt((GPS_pos(1,i)-x(1))^2+(GPS_pos(2,i)-x(2))^2+(GPS_pos(3,i)-x(3))^2);
                % Compute new time delay
                  tau=p/c;
                  d_p=abs(tau-tau0(i));
                  tau0(i)=tau;
             end
          end
        % Compute the return values
          GPS_pos=GPS_pos';           % GPS Satellite position    
          GPSclc_corr=GPSclc_corr.*c; % GPS Satellite clock correction 
          GPS_Tgd=c*GPS_Tgd*0;          % GPS instrumental delay
        end
        

        function loadObsFile(this,t)
        % FUNCTION
        %   Loads the next broadcast observation file 
        %---------------------------------------------------------------- 
        % file name
          doy = dayofyear(t);
          year = num2str(get(t,'year'));
          fileName = strcat('brdc',num2str(doy),'0.',year(3:4),'n');
          filePath = strcat(this.settings.ws,'inputs/brdc/', fileName);        
        % Open the observation file
          % Create classes
            BEph = CoordEph(); 
            brd_eph_data = StdOrb; 
          % Read Broadcast Ephemerides
            [BEph, OrbPar] = ReadEphBroadcast(BEph,filePath,15*60);
          % Read group delay 
            stTgd = getTgd(filePath);
          % Compute standard orbits  
            this.brdEphData = CompStdOrb(brd_eph_data,BEph,stTgd); 
            
            this.current_t = t;
          
        end        
    end
    
end

