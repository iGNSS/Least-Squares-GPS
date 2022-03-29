classdef LSSolver < handle
    
    properties
        settings
        obsFileObj
        brdEphFileObj    
        measModelObj
        dynModelObj
        timeRefSysObj
        outFileObj
    end
                                                  
  methods
      function this=LSSolver()
          
      
      end

 function [x]=kinematicPositioning(this,t,z,sv_no)
        %-----------------------------------------------------------------            
        % FUNCTION
        %   Computes position, velocity and receiver clock bias
        % INPUTS
        %   t    : date and time, such as [year month day hour minute second] 
        %   z    : pseudo code observation
        %   sv_no: vehicle number of GPS satellites
        % OUTPUTS 
        %  x     : vector including estimated position and receiver clock 
        %          correction, such ad;    x=[x y z clc_rec]
        %--------------------------------------------------------------                 
            
            %c=299792458; %speed of light in meters
            x=zeros(4,1);
            for i=1:50
              % Compute satellite position
                [GPS_pos,GPSclc_corr]=this.brdEphFileObj.getGpsSatParam( ...
                                      t,x(1:3),sv_no,x(4));  
                                  
              % Compute geometric distance,d
               zg = sqrt((power(GPS_pos(:,1)-x(1),2)+ ...
                           power(GPS_pos(:,2)-x(2),2)+ ...
                           power(GPS_pos(:,3)-x(3),2)));
              % Measurement partials      
                H = zeros(length(z),length(x));
                for j = 1:length(z)
                    H(j,1) =-(GPS_pos(j,1)-x(1))/zg(j);
                    H(j,2) =-(GPS_pos(j,2)-x(2))/zg(j);
                    H(j,3) =-(GPS_pos(j,3)-x(3))/zg(j);
                    H(j,4) = 1;
                end 
              % Approximate range
                zc=zg+ones(size(zg)).*(x(4)-GPSclc_corr);% computed range
              % Least Square Estimation
                dz=(z-zc);
                dx=pinv(H'*H)*H'*dz;
                x=x+dx;
                
            end
           
 end
  function solverObj=getInstance(this,varargin)
        %-----------------------------------------------------------------            
        % FUNCTION
        %   Creates the solver object as an instance
        % INPUTS
        %   varargin : cell structure including input parameters. The class
        %              can be initiated as follows;
        %
        %              For Code and GRAPHIC measurements
        %                 obj=LSSolver(settings,obsFileObj,outFileObj,brdEphObj)
        %              For navigation solution measurements
        %                 obj=LSSolver(settings,obsFileObj,outFileObj)    
        %              
        %              where
        %                settings   : data structure including filter settings
        %                obsFileObj: Handle to instance of 'OrbObsFile' class
        %                brdEphObj: Handle to instance of 'BrdEphObj' class  
        %                outFileObj: Handle to instance of 'OrbFilterOutputFile'
        %                            to store filter outputs
        % OUTPUTS 
        %   ModelObj : Handle to child class instance
        %--------------------------------------------------------------  
      
          solverObj=this;
          settings=varargin{1};
        % Store solver settings
          solverObj.settings=settings;
        % Create Measurement Model Objects
          baseMeasModelObj=OrbMeasurementModel();
          if strcmp(settings.obsType,'NavSol')
             solverObj.measModelObj=baseMeasModelObj.getInstance(...
                                           settings,varargin{2});
          elseif strcmp(settings.obsType,'Graphic') || ...
                 strcmp(settings.obsType,'Code')
                 solverObj.measModelObj=baseMeasModelObj.getInstance( ...
                                         settings,varargin{2},varargin{4});
          else
             error ('Input measurement model type can not be recognised')
          end
       
        % Set data files
          solverObj.obsFileObj=varargin{2};     
          solverObj.outFileObj=varargin{3};
          n=length(varargin);  
          if n==4
            solverObj.brdEphFileObj=varargin{4};
          end  
       
        end

        function run(this)
      
            if strcmp(this.settings.obsType,'Code') || ...
                        strcmp(this.settings.obsType,'Graphic')
            % Read the observation file that corresponds with the beginning
            % time of the settings
              init_flag=true;
              while init_flag
                % Read observation data
                  if strcmp(this.settings.obsType,'Code')
                     [orbEDobj,health_flag]= ...
                                       this.obsFileObj.getNextObs('Code');
                  elseif strcmp(this.settings.obsType,'Graphic')
                     [orbEDobj,health_flag]= ...
                                       this.obsFileObj.getNextObs('Graphic');  
                  end

                  orbepo=orbEDobj.epochDateTime;
                  if get(DateTime(orbepo(1),orbepo(2),orbepo(3),orbepo(4),orbepo(5),orbepo(6)),'MJD')...
                          <get(this.settings.beginDateTime,'MJD')
                      disp(strcat('Search the observation file for the solver beginning,obs time:', ...
                                   num2str(orbEDobj.epochDateTime)))
                     
                      continue;
                  else
                      init_flag=false;
                      disp(orbepo)
             
                  end
              end
            end
        % Run the least squares loop
          while 1
                while 1
                    
                  % Get the observations
                    [nextEDobj,eof_flag]= ...
                          this.obsFileObj.getNextObs(this.settings.obsType);
                    % Control end of file
                      if eof_flag; break ; end                        
                  % Date and time of the current observation epoch
                    nextDTobsObj=DateTime(nextEDobj.epochDateTime);
                    %Debug
                    %nextDTobsObj
                    % compute the position
                    x=this.kinematicPositioning(nextEDobj.epochDateTime, ...
                                                      nextEDobj.obs.C1,...
                                                    nextEDobj.svn);
                    %plot(sqrt(x(1)^2+x(2)^2+x(3)^2))
                    plot(x(1:3)) % in meters
                    % Check the solver end time
                    if get(nextDTobsObj,'MJD')>=get(this.settings.endDateTime,'MJD')
                      disp('END OF LS SOLVER')
%                       close('all')
                      return;
                    end
          
                end
              % Control end of file
                if eof_flag; fclose('all');break ; end 
          end
          
        
        end
  end
end

