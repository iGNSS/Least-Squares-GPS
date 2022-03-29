classdef OrbCodeMeasModel < OrbMeasurementModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function this=OrbCodeMeasModel()
            
        end
        
        function [updSVobj,updEDobj,minObsFlag]=dataEditing(this,...
                                            currEDobj,predSVobj)
        %-----------------------------------------------------------------                                        
        % FUNCTION
        %   Function search for the observation, elevation angles of which
        %   are below the threshold and outliers than removes them from the data 
        %   set. In addition, function controls the minumum number of observations. 
        % INPUTS
        %   currEDobj : instance of 'OrbEpochData' class including observations
        %               GPS ephemerides and other auxiliary data at current
        %               time   
        %   predSVobj : instance of 'OrbStateVector' class including predicted 
        %               state vector parameters  
        % OUTPUTS
        %   updSVobj  : instance of 'OrbStateVector' class including updated 
        %               state vector parameters and covariance matrix           
        %   updEDobj : instance of 'OrbEpochData' class including updated 
        %              observations, GPS ephemerides and other auxiliary data
        %   minObsFlag  : If number of minimum observations is under the 
        %                 threshold after editing process, minObsFlag takes 
        %                 value "1" else value "0". If minObsFlag takes the
        %                 value "0", updSVobj and updEDobj variables are
        %                 set to empty "[]"
        %-----------------------------------------------------------------      
        
        % Get the GPS ephemerides and update the epoch data
          [gpsPos,gpsClcBias]=this.brdEphObj.getGpsSatParam(...
                                          currEDobj.epochDateTime, ...
                                          predSVobj.position, ...
                                          currEDobj.svn);
          % Update epoch data to include GPS ephemerides data
            currEDobj.gpsPos=gpsPos; currEDobj.gpsClcCorr=gpsClcBias;        
          
        % Remove the observations, elevation angles of which under threshold
          [updEDobj,maskFlag,minObsFlag]=this.removeObsBelowElevThresh(...
                                               predSVobj,currEDobj);
          % Control the minumum number of observation                                 
            if minObsFlag==1
               updSVobj=[];updEDobj=[];return;
            end                                           
        % Outlier Detection    
          [updSVobj,updEDobj,outlierFlag,minObsFlag]=this.removeOutliers( ...
                                                     predSVobj,updEDobj);
          % Control the minumum number of observation
            if minObsFlag==1
               updSVobj=[];updEDobj=[];return;

            end
            
        
        end
        
        function [updEDobj,maskFlag,minObsFlag]=removeObsBelowElevThresh( ...
                                                this,SVobj,EDobj)
        %-----------------------------------------------------------------                                        
        % FUNCTION
        %   Function serch for the observation, elevation angles of which
        %   are below the threshold, than removes them from the data set
        % INPUTS
        %   EDobj : instance of 'OrbEpochData' class including observations
        %           GPS ephemerides and other auxiliary data
        %   SVobj : instance of 'OrbStateVector' class including state
        %           vector variables  
        % OUTPUTS
        %   updEDobj : instance of 'OrbEpochData' class including updated 
        %              observations, GPS ephemerides and other auxiliary data
        %   maskFlag : If elevation angle of any observations are below the 
        %              threshold, maskFlag takes value '1', else value '0'
        %   minObsFlag  : If number of minimum observations is under the 
        %                 threshold, minObsFlag takes value '1' else value '0'
        %----------------------------------------------------------------- 
        
        % Set Parameters
          minNumObs=this.filtSet.dataEditing.minNumObs;
          elvThreshold=this.filtSet.dataEditing.elevationThreshold; % in degree
        % Rotation matrix from ECEF system to body system
          ecef2body=(this.timeRefSysObj.orb2geo(...
                                        SVobj.position(1),...
                                        SVobj.position(2),...
                                        SVobj.position(3),...
                                        SVobj.velocity(1),...
                                        SVobj.velocity(2),...
                                        SVobj.velocity(3)))';
        % Find Elevation angel for each receiver and remove the
        % observations that below the threshold
          maskFlag=0; remInd=zeros(length(EDobj.svn));
          for i=1:length(EDobj.svn)
            % Position difference vector between satellite and receiver
              pos_diff = [EDobj.gpsPos(i,1)-SVobj.position(1)
                          EDobj.gpsPos(i,2)-SVobj.position(2)
                          EDobj.gpsPos(i,3)-SVobj.position(3)];
            % Compute the unit vector
              e_GPS=pos_diff/norm(pos_diff);
            % Convert the GPS position unit vector given in ECEF system to 
            % body system
              e_body=ecef2body*e_GPS;
            % elevation angle
              elevationangle=atan2(-e_body(1),sqrt(e_body(3)^2+e_body(2)^2));
              elevationangle=elevationangle*180/pi;
            % Compare to threshold
              if elevationangle<elvThreshold
                 remInd(i)=1;
                 maskFlag=1;
              end
          end
          remInd=logical(remInd);
        % Update epoch data object, state vector and covariance matrix
          if maskFlag==1
             remSatId=EDobj.svn(remInd);  
             [updEDobj]=this.removeObsFromDataSet(EDobj,remSatId); 
          else
             updEDobj=OrbEpochData(EDobj.obs,EDobj.svn,EDobj.gpsPos,...
                                   EDobj.gpsClcCorr,EDobj.epochDateTime,...
                                   EDobj.hFlag);                   
          end
        % Control the minimum number of observation 
          if (length(updEDobj.svn))<minNumObs
             minObsFlag=1;
          else
             minObsFlag=0; 
          end            
        end
        
        function [updSVobj,updEDobj,outlierFlag,minObsFlag]=removeOutliers( ...
                                                   this,SVobj,EDobj)
        %-----------------------------------------------------------------             
        % FUNCTION
        %   Function removes the outlier from observation set and updates
        %   both epoch data and state vector. Besides, function make a priory 
        %   estimation for the receiver clock bias and its statistics, then 
        %   update the both state vector and state cvariance matrix
        % INPUTS
        %   EDobj : instance of 'OrbEpochData' class including observations
        %           GPS ephemerides and other auxiliary data
        %   SVobj : instance of 'OrbStateVector' class including state
        %           vector variables  
        % OUTPUTS
        %   updEDobj : instance of 'OrbEpochData' class including updated 
        %              observations, GPS ephemerides and other auxiliary data
        %   updSVobj : instance of 'OrbStateVector' class including updated 
        %              state vector variables    
        %   outlierFlag : If any outlier detected, outlierFlag takes value 
        %                 '1', else value '0'
        %   minObsFlag  : If number of minimum observations is under the 
        %                 threshold, minObsFlag takes value '1' else value '0'
        %----------------------------------------------------------------- 

        % Set parameters 
          obsSISRE=this.filtSet.stat.std.obsSISRE;
          stdObsNoise=this.filtSet.stat.std.measurementNoise;
          minNumObs=this.filtSet.dataEditing.minNumObs;
          stdFactor=this.filtSet.dataEditing.outlierFactor;
        % Get apriory observations
          z=EDobj.obs.C1;
        % Get GPS observations
          gpsPos=EDobj.gpsPos; gpsClcBias=EDobj.gpsClcCorr;
        % Get the receiver positions
          xp=SVobj.position;
        % Get a priory geomatric range
          geoRange = sqrt((gpsPos(:,1)-xp(1)).^2+ ...
                          (gpsPos(:,2)-xp(2)).^2+ ...
                          (gpsPos(:,3)-xp(3)).^2);
        % Compute expected uncertanity of the biased measurements (modelled 
        % measurements without receiver clock bias )
          % Expected uncertanity of each biased observation due to only 
          % position and ambiguity biases projected on line of sight vector 
            [H]=this.getDesignMatrix (EDobj,SVobj);
            [rH,cH]=size(H); 
            Pp=SVobj.stateCov;
            H(:,4:end-1)=[];  Pp(:,4:end-1)=[]; Pp(4:end-1,:)=[];
            for i=1:length(EDobj.svn)
                sigmaResWoRcb(i,1)=sqrt(H(i,:)*Pp*H(i,:)'+ ...
                                         stdObsNoise^2+obsSISRE^2);  
            end
        % Search for the bad measurement
          outlierFlag=0;  remObs=[];
          while 1
            % Computed observable without receiver clock correction
              [zcWoRcb]=geoRange-gpsClcBias; 
            % Estimated standart deviation of receiver clock bias
              sigmaRcb=sqrt(1/(sum(1./(sigmaResWoRcb.^2))));                     
            % Priory estimation of receiver clock bias 
              recClcBias=(sigmaRcb^2)* ...
                            sum((z-zcWoRcb)./(sigmaResWoRcb.^2));
            % Residuals
              res=z-recClcBias-zcWoRcb;  
            % Predicted uncertanity of measurment residuals
              sigmaRes(:,1)=sqrt(sigmaResWoRcb.^2+sigmaRcb^2);       
            % Normalized residuals
              resNorm=res./sigmaRes;
            % RMS of normalized residuals
              rmsResNorm=std(resNorm);   
            % Compare to the threshold
              if rmsResNorm <= stdFactor
                 break
              else
                % Set the outlier flag
                  outlierFlag=1;                  
                % Compute RMS of each observation subset
                  for i=1:length(z);
                      % Create new observation set by excluding one 
                      % measurement from set
                        subZ=z; subZ(i)=[];
                        subSigmaResWoRcb=sigmaResWoRcb;subSigmaResWoRcb(i)=[];
                        subZcWoRcb=zcWoRcb; subZcWoRcb(i)=[];
                      % Compute mean receiver clock bias                              
                        % Estimated standart deviation of receiver clock bias
                          subSigmaRcb(i,1)=sqrt(1/ ...
                                     (sum(1./(subSigmaResWoRcb.^2))));                     
                        % Priory estimation of receiver clock bias 
                          subRecClcBias(i,1)=(subSigmaRcb(i)^2)* ...
                                          sum((subZ-subZcWoRcb)./...
                                          (subSigmaResWoRcb.^2));
                      % Residuals
                        res=subZ-subRecClcBias(i)-subZcWoRcb;       
                      % Predicted uncertanity of measurment residuals
                        sigmaRes=sqrt(subSigmaResWoRcb.^2+ ...
                                       subSigmaRcb(i)^2);                        
                      % Normalized residuals
                        resNorm=res./sigmaRes;
                      % RMS of normalized residuals
                        rmsResNorm(i)=sqrt(sum(resNorm.^2)/(length(resNorm)-1));                        
                  end
                % Find the observation subset with best RMS to remove bad
                % measurement
                  [min_val,min_ind]=min(rmsResNorm);
                  recClcBias=subRecClcBias(min_ind);
                  sigmaRcb=subSigmaRcb(min_ind);
                % Update internal processing parameters
                  sigmaResWoRcb(min_ind)=[]; subRecClcBias=[min_ind];
                  subSigmaRcb(min_ind)=[];
                % Update the variables for the next loop
                  remObs=[remObs,z(min_ind)]; % observations that will be removed
                  geoRange(min_ind)=[]; z(min_ind)=[]; 
                  gpsClcBias(min_ind)=[];
                % If bad quality measurements removed, stop the search              
                  if rmsResNorm(min_ind)<stdFactor
                     break 
                  else
                      rmsResNorm=[];
                  end
              end
          end
        % Set the updated output epoch data and state vector objects
          if ~isempty(remObs)
             % Find the satellite vehicle number of observations that will 
             % be removed
               remInd=logical(0);
               for i=1:length(remObs)
                   remInd=logical(remInd+(remObs(i)==EDobj.obs.C1));
               end
               remSatId=EDobj.svn(remInd);
             % Update the epoch data, state vector and covariance
               [updEDobj]=this.removeObsFromDataSet( EDobj,remSatId);  
          else
             updEDobj=OrbEpochData(EDobj.obs,EDobj.svn,EDobj.gpsPos, ...
                                   EDobj.gpsClcCorr,EDobj.epochDateTime,...
                                   EDobj.hFlag);               
          end
          % Initialize output State Vector
            updSVobj=OrbStateVector(SVobj.stateVector, ...
                                    this.filtSet.enableDynModParam, ...
                                    this.filtSet.obsType,...
                                    SVobj.stateDateTime);    
            updSVobj.stateCov=SVobj.stateCov;                   
        % Set the a priory estimated receiver clock bias parameter
          % Find the index of receiver clock bias
            ind=updSVobj.stateVector==updSVobj.recClcBias;
          % Set the new receiver clock bias 
            updSVobj.recClcBias=recClcBias;     
            updSVobj.stateVector(ind)=recClcBias;
            updSVobj.stateCov(ind,:)=zeros;updSVobj.stateCov(:,ind)=zeros;
            updSVobj.stateCov(ind,ind)=sigmaRcb^2;
        % Control the observation count
          if length(updEDobj.svn)<minNumObs
             minObsFlag=1;
          else
             minObsFlag=0;
          end             
            
        end
        
        function [H]=getDesignMatrix(this,EDobj,SVobj)
        %-----------------------------------------------------------------             
        % FUNCTION
        %   Function computes the design matrix for the GRAPHIC measurement 
        %   model
        % INPUTS
        %   EDobj : instance of 'OrbEpochData' class including observations
        %           GPS ephemerides and other auxiliary data
        %   SVobj : instance of 'OrbStateVector' class including state
        %           vector variables       
        % OUTPUTS
        %   H     : Computed design matrix
        %            
        %----------------------------------------------------------------- 
        
        % Partial derivatives of measurement model with respect 
        % to state vector parameters
          % Compute geometric range
            GPS_pos=EDobj.gpsPos; xp=SVobj.position;
            geoRange = sqrt((GPS_pos(:,1)-xp(1)).^2+ ...
                            (GPS_pos(:,2)-xp(2)).^2+ ...
                            (GPS_pos(:,3)-xp(3)).^2);
          % Measurement partials with respect to position, velocity, 
          % receiver clock bias  and ambiguity biases components    
             m=length(SVobj.stateVector);
             n=length(EDobj.obs.C1);
             H = zeros(n,m);
            for i = 1:n
                H(i,1) =-(GPS_pos(i,1)-xp(1))/geoRange(i);
                H(i,2) =-(GPS_pos(i,2)-xp(2))/geoRange(i);
                H(i,3) =-(GPS_pos(i,3)-xp(3))/geoRange(i);
                H(i,m-n) = 1;
            end             
        end
        
        function [zc]=getCompObs(this,EDobj,SVobj)
        % FUNCTION
        %   Function computes the observations based on the measurement model 
        % INPUTS
        %   EDobj : instance of 'OrbEpochData' class including observations
        %           GPS ephemerides and other auxiliary data
        %   SVobj : instance of 'OrbStateVector' class including state
        %           vector variables       
        % OUTPUTS
        %   zc    : Computed observations
        %            
        %----------------------------------------------------------------- 
        
        % Compute geometric range
            GPS_pos=EDobj.gpsPos; xp=SVobj.position;
            geoRange = sqrt((GPS_pos(:,1)-xp(1)).^2+ ...
                            (GPS_pos(:,2)-xp(2)).^2+ ...
                            (GPS_pos(:,3)-xp(3)).^2);
        % Computed GRAPHIC observable
          RECclc_corr=SVobj.recClcBias; % Receiver clock bias
          GPSclc_corr=EDobj.gpsClcCorr; % Satellite clock bias
          zc=geoRange+RECclc_corr-GPSclc_corr;              
            
        end
        
        function [updEDobj]=removeObsFromDataSet(this,EDobj,svn)
        % FUNCTION
        %   Function removes the observation given by satellite vehicle
        %   number from the data set
        % INPUTS
        %   EDobj : instance of 'OrbEpochData' class including observations
        %           GPS ephemerides and other auxiliary data
        % OUTPUTS
        %   updEDobj : instance of 'OrbEpochData' class including updated 
        %              observations, GPS ephemerides and other auxiliary data
        %-----------------------------------------------------------------             

        % Find the index of observations that will be removed
          remInd=logical(0);
          for i=1:length(svn)
              remInd=logical(remInd+(svn(i)==EDobj.svn));
          end
        % Update the epoch data
          updEDobj=OrbEpochData(EDobj.obs,EDobj.svn,EDobj.gpsPos, ...
                                EDobj.gpsClcCorr,EDobj.epochDateTime,...
                                EDobj.hFlag); 
          updEDobj.obs.C1(remInd)=[];
          updEDobj.svn(remInd)=[];
          updEDobj.gpsPos(remInd,:)=[];
          updEDobj.gpsClcCorr(remInd)=[];
                      
        end        
    end
    
end

