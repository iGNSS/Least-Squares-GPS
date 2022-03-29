classdef OrbGraphicMeasModel < OrbMeasurementModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        bodyToGps % vector from body center to antenna center in ECEF
    end
    
    methods
        function this=OrbGraphicMeasModel()
            
        end
        
        function [updSVobj,updEDobj,minObsFlag]=dataEditing(this,...
                                            prevEDobj,currEDobj,predSVobj)
        %-----------------------------------------------------------------                                        
        % FUNCTION
        %   Function search for the observation, elevation angles of which
        %   are below the threshold, outliers than removes them from the data 
        %   set. Besides updates ambiguity bias parameters in predicted covariance 
        %   matrix and state vector for newly allocated and untracked satellites. 
        %   In addition, function controls the minumum number of observation.  
        % INPUTS
        %   prevEDobj : instance of 'OrbEpochData' class including observations
        %               GPS ephemerides and other auxiliary data at previous             
        %               time of epoch
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
        
        % Set initial parameters
          minObsFlag=0;
        % Get the GPS ephemerides and update the epoch data
          [gpsPos,gpsClcBias,gpsTgd]=this.brdEphObj.getGpsSatParam(...
                                          currEDobj.epochDateTime, ...
                                          predSVobj.position, ...
                                          currEDobj.svn,...
                                          predSVobj.recClcBias);
          % Update epoch data to include GPS ephemerides data
            currEDobj.gpsPos=gpsPos; currEDobj.gpsClcCorr=gpsClcBias;  
            currEDobj.gpsTgd=gpsTgd;
        % Remove the ambiguity bias parameters which belongs to the satellites
        % that no longer tracked from state vector and covariance matrix,
        % and update the ambiguity bias parameters that are newly allocated
          [updSVobj]=this.stateAmbBiasUpd(prevEDobj,currEDobj,predSVobj);
        % Transformation vector from body center to antenna center
          orb2geo=this.timeRefSysObj.orb2geo(predSVobj.position(1),...
                                             predSVobj.position(2),...
                                             predSVobj.position(3),...
                                             predSVobj.velocity(1),...
                                             predSVobj.velocity(2),...
                                             predSVobj.velocity(3));
          this.bodyToGps=orb2geo*this.filtSet.dataEditing.AntOffFromSatBaseLoc(:);          
          
        % Outlier detection
          if this.filtSet.dataEditing.mode==1 % Recursive outlier detection mode
             % Remove the observations, elevation angles of which under threshold
              [updSVobj,updEDobj,maskFlag,minObsFlag]=this. ...
                              removeObsBelowElevThresh(updSVobj,currEDobj);
              
             % Control the minumum number of observation                                 
               if minObsFlag==1
                  updSVobj=[];updEDobj=[];return;
               end                                                                                      
             % Outlier Detection   

                       %[updSVobj,updEDobj,outlierFlag,minObsFlag]=this.removeOutliers( ...
              %                                      updSVobj,updEDobj);
              %                                      % original
              %Debug
              outlierFlag=0;
              minObsFlag=0;
              %DataEditUpdSV=updSVobj
              %flag=minObsFlag

          elseif this.filtSet.dataEditing.mode==2 % Robust filtering mode
             % If the filter will executed in robust filtering mode, do
             % nothing, only return inputs as outputs
              %[updSVobj,updEDobj,minObsFlag]=this.priorRecClcEst( updSVobj,currEDobj);
             % Set the outputs    
               updEDobj=OrbEpochData(currEDobj.obs,currEDobj.svn,currEDobj.gpsPos, ...
                                currEDobj.gpsClcCorr,currEDobj.epochDateTime,...
                                currEDobj.hFlag,currEDobj.gpsTgd);    
               Ptemp=updSVobj.stateCov;                    
               updSVobj=OrbStateVector(updSVobj.position,updSVobj.velocity, ...
                                  updSVobj.atmDragCoef,updSVobj.solarRadCoef, ...
                                  updSVobj.empAccel,updSVobj.corelTime,...
                                  updSVobj.recClcBias,updSVobj.ambBias,...
                                  updSVobj.stateDateTime);
               updSVobj.stateCov=Ptemp;           
          else
              error ('Wrong input for filtSet.dataEditing.mode')
          end
        % Control the minumum number of observation
          if minObsFlag==1
             updSVobj=[];updEDobj=[];return;
          end          
        end
        
        function [updSVobj,updEDobj,maskFlag,minObsFlag]=removeObsBelowElevThresh( ...
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
        %   updSVobj  : instance of 'OrbStateVector' class including updated 
        %               state vector parameters and covariance matrix           
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
            %Debug all obs omitted
            %SVobjposX=SVobj.position(1)
            %gpsposX=EDobj.gpsPos(i,1)
              pos_diff = [EDobj.gpsPos(i,1)-SVobj.position(1)
                          EDobj.gpsPos(i,2)-SVobj.position(2)
                          EDobj.gpsPos(i,3)-SVobj.position(3)];
            % Compute the unit vector
              e_GPS=pos_diff/norm(pos_diff);
            % Convert the GPS position unit vector given in ECEF system to 
            % body system
            %Debug
            %ECEFunitvector=ecef2body
            %gpsUnitVector=e_GPS
              e_body=ecef2body*e_GPS;
            % elevation angle
              elevationangle=atan2(-e_body(1),sqrt(e_body(3)^2+e_body(2)^2));
              elevationangle=elevationangle*180/pi;
              %elvThreshold=elvThreshold
              %display(elevationangle<elvThreshold)
            % Compare to threshold
              if elevationangle<elvThreshold
                  %display('within elvAngle<elvThreshold loop')
                 %Debug
                 %elevasyon=elevationangle
                 %sinir=elvThreshold
                 remInd(i)=1;
                 %Debug
                 %remInd=remInd
                 maskFlag=1;
              end
          end
          remInd=logical(remInd);
        % Update epoch data object, state vector and covariance matrix
          if maskFlag==1
             remSatId=EDobj.svn(remInd);  
             %Debug Empty SVN
             %beforeRemoveObsSVN=EDobj.svn
             %beforeRemoveObsEphemeris=EDobj.obs
             [updSVobj,updEDobj]=this.removeObsFromDataSet(SVobj,EDobj,remSatId);
             %afterREmoveObsSVN=updEDobj.svn
             %afterREmoveObsEphemeris=updEDobj.obs
          else
             updEDobj=OrbEpochData(EDobj.obs,EDobj.svn,EDobj.gpsPos,...
                                   EDobj.gpsClcCorr,EDobj.epochDateTime,...
                                   EDobj.hFlag,EDobj.gpsTgd);  
             
             updSVobj=OrbStateVector(SVobj.position,SVobj.velocity,...
                                     SVobj.atmDragCoef,SVobj.solarRadCoef, ...
                                     SVobj.empAccel,SVobj.corelTime,...
                                     SVobj.recClcBias,SVobj.ambBias,...
                                     SVobj.stateDateTime);
             updSVobj.stateCov=SVobj.stateCov;                     
          end
        % Control the minimum number of observation 
        %Debug minObs problem
        %NoSvn=length(updEDobj.svn)
          if (length(updEDobj.svn))<minNumObs
             minObsFlag=1;
          else
             minObsFlag=0; 
          end            
        end
        
        function [updSVobj,ambUpdFlag]=stateAmbBiasUpd(this,prevEDobj, ... 
                                                    currEDobj,predSVobj)
        %-----------------------------------------------------------------             
        % FUNCTION
        %   Function removes the ambiguity bias data from predicted state 
        %   vector and covariance belongs the satellite that no longer in 
        %   track at current epoch and update the data that are newly allocated.
        % INPUTS
        %   prevEDobj : instance of 'OrbEpochData' class including observations
        %               GPS ephemerides and other auxiliary data at previous             
        %               time of epoch
        %   currEDobj : instance of 'OrbEpochData' class including observations
        %               GPS ephemerides and other auxiliary data at current
        %               time   
        %   predSVobj : instance of 'OrbStateVector' class including predicted 
        %               state vector parameters  
        % OUTPUTS
        %   updSVobj  : instance of 'OrbStateVector' class including updated 
        %               state vector parameters and covariance matrix         
        %   ambUpdFlag: ambUpdFlag takes value '0', if there exist no update,
        %               else '1'
        %-----------------------------------------------------------------               
            
        %-- Set constants
          c=299792458; % speed of light in vacum (m/sec)
          lamda_L1=(c/1575.42e6); % wavelength of L1 carrier
        %-- Compute the updated ambiguity bias of the state vector
          ambUpdFlag=0;
          if ~isequal(prevEDobj.svn,currEDobj.svn)
              ambUpdFlag=1;

              % Update Ambiguity
                % Clear the ambiguity bias parameters from state that are  
                % not exist in observations at current epoch and approximate  
                % ambiguity biases for the new observations 
                  for i=1:length(currEDobj.svn)
                      match_sv=find(prevEDobj.svn==currEDobj.svn(i));
                      if ~isempty(match_sv) % Keep the exist one
                          updAmbBias(i,1)=predSVobj.ambBias(match_sv);
                      else                  % Initialize newly allocated
                          updAmbBias(i,1)=(currEDobj.obs.C1(i)- ...
                                           currEDobj.obs.L1(i)*lamda_L1);
                      end
                  end
          end
        %-- Compute the updated ambiguity bias of the state covariance matrix    
          if ambUpdFlag==1
             % Initialize parameters
               n_SV=length(predSVobj.stateVector); % length of old state vecctor
               n_Amb=length(predSVobj.ambBias);    % length of old ambiguity bias
               updP=predSVobj.stateCov;
             % Remove the components of satellites that are no longer observed 
               % Find the index of components 
                 rem_ind=[];
                 for i=1:length(prevEDobj.svn)
                     coln=find(currEDobj.svn==prevEDobj.svn(i));
                     if isempty(coln)
                        rem_ind=[rem_ind,i];
                     end
                 end
               % Remove the components  
                 if ~isempty(rem_ind)
                   rem_ind=(n_SV-n_Amb)+rem_ind;
                   updP(rem_ind,:)=[];updP(:,rem_ind)=[];
                 end              
             % Extend covariance matrix for newly allocated satellites
               % Find the index of components that will be initialized             
                 add_ind=[];
                 for i=1:length(currEDobj.svn)
                     coln=find(prevEDobj.svn==currEDobj.svn(i));
                     if isempty(coln)
                          add_ind=[add_ind,i];
                      end
                 end  
               % Add the initial variance for newly allocated satellites
                 add_ind=(n_SV-n_Amb)+add_ind;
                 for i=1:length(add_ind)
                     ind=add_ind(i);
                     sz=length(updP);
                     updP=[updP;zeros(1,sz)]; updP=[updP,zeros(sz+1,1)];
                     updP(ind+1:end,:)=updP(ind:end-1,:);
                     updP(:,ind+1:end)=updP(:,ind:end-1);
                     updP(ind,:)=0; updP(:,ind)=0;
                     updP(ind,ind)=this.filtSet.stat.std.init.ambiguityBias^2;
                 end                 
          end
        %-- Set the output parameters
            if ambUpdFlag==1
               updSVobj=OrbStateVector(predSVobj.position,predSVobj.velocity, ...
                                       predSVobj.atmDragCoef,predSVobj.solarRadCoef, ...
                                       predSVobj.empAccel,predSVobj.corelTime,...
                                       predSVobj.recClcBias,updAmbBias,...
                                       predSVobj.stateDateTime);
               updSVobj.stateCov=updP;
            else
               updSVobj=OrbStateVector(predSVobj.position,predSVobj.velocity, ...
                                       predSVobj.atmDragCoef,predSVobj.solarRadCoef, ...
                                       predSVobj.empAccel,predSVobj.corelTime,...
                                       predSVobj.recClcBias,predSVobj.ambBias,...
                                       predSVobj.stateDateTime);
               updSVobj.stateCov=predSVobj.stateCov;
            end
        end
        
        function [updSVobj,updEDobj,minObsFlag]=priorRecClcEst(...
                                                       this,SVobj,EDobj)
        %-----------------------------------------------------------------             
        % FUNCTION 
        %   Function checks the minimum number of observation using a prior 
        %   outlier detection and make a prior estimation for the receiver 
        %   clock bias and its statistics, then update the both state vector 
        %   and state covariance matrix
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
        %-----------------------------------------------------------------
        
        % Set parameters
          obsSISRE=this.filtSet.stat.std.obsSISRE;
          stdObsNoise=this.filtSet.stat.std.measurementNoise;
        % Get apriory observations
          z=EDobj.obs.Graphic;
        % Get GPS observations
          gpsPos=EDobj.gpsPos; gpsClcBias=EDobj.gpsClcCorr;
        % Get a priory ambiguity bias parameters
          ambBias=SVobj.ambBias;
        % Get the receiver positions
          xp=SVobj.position;
          xp=xp+this.bodyToGps; %body system to GPS antenna system
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
            %Debug
            %mantik=isempty(SVobj.stateCov)
            %goster=SVobj.stateCov
            if isempty(SVobj.stateCov)
                Pp=zeros(rH,cH);
            else
                Pp=SVobj.stateCov;
            end
            %display(Pp)
            H(:,4:cH-rH-1)=[];  Pp(:,4:cH-rH-1)=[]; Pp(4:cH-rH-1,:)=[];
            %Debug
            %display(H)
            for i=1:length(EDobj.svn)
                sigmaResWoRcb(i,1)=sqrt(H(i,:)*Pp*H(i,:)'+ ...
                                         stdObsNoise^2+obsSISRE^2);  
            end  
          % Estimated standart deviation of receiver clock bias
            sigmaRcb=sqrt(1/(sum(1./(sigmaResWoRcb.^2))));  
        % Compute a priory clock estimation
          % Computed observable without receiver clock correction
            [zcWoRcb]=geoRange-gpsClcBias-ambBias./2;      
          % Priory estimation of receiver clock bias based on median
          % statistics
            resWoRcb=z-zcWoRcb; % residuals without receiver clock bias 
            recClcBias=median(resWoRcb); 
        % Compare the receiver clock bias with the expected uncertanity 
          n_obs=length(z);
          minNumObs=this.filtSet.dataEditing.minNumObs;
          % Find the observations, standart deviation of which exceeds the
          % 3 sigma threshold
            res=z-recClcBias-zcWoRcb; 
            out_len=length(find(abs(res)> 3*sigmaRcb));
 'Outlier control'
 res
 3*sigmaRcb
 n_obs
 out_len
 
%             if (n_obs-out_len)<minNumObs
%                 minObsFlag=1;
%             else
%                 minObsFlag=0;
%             end     
minObsFlag=0;
            
        % Set the outputs    
          updEDobj=OrbEpochData(EDobj.obs,EDobj.svn,EDobj.gpsPos, ...
                                EDobj.gpsClcCorr,EDobj.epochDateTime,...
                                EDobj.hFlag);               
          updSVobj=OrbStateVector(SVobj.position,SVobj.velocity, ...
                                  SVobj.atmDragCoef,SVobj.solarRadCoef, ...
                                  SVobj.empAccel,SVobj.corelTime,...
                                  SVobj.recClcBias,SVobj.ambBias,...
                                  SVobj.stateDateTime);
          updSVobj.stateCov=SVobj.stateCov; 
        % Set the a priory estimated receiver clock bias parameter
          % Find the index of receiver clock bias
            ind=updSVobj.stateVector==updSVobj.recClcBias;
          % Set the new receiver clock bias 
            tempRecCov=updSVobj.stateCov(ind,ind);
%            updSVobj.recClcBias=recClcBias;     
%            updSVobj.stateVector(ind)=recClcBias;
%            updSVobj.stateCov(ind,:)=zeros;updSVobj.stateCov(:,ind)=zeros;
%            updSVobj.stateCov(ind,ind)=tempRecCov;                     
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
          z=EDobj.obs.Graphic;
        % Get GPS observations
          gpsPos=EDobj.gpsPos; gpsClcBias=EDobj.gpsClcCorr;gpsTgd=EDobj.gpsTgd;
        % Get a priory ambiguity bias parameters
          ambBias=SVobj.ambBias;
        % Get the receiver positions
        %Debug
          xp=SVobj.position;
          xp=xp+this.bodyToGps; %body system to GPS antenna system
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
            %Debug
            %rH=rH
            %cH=cH
            %ambBiasVect=SVobj.ambBias
            %StateVector=SVobj.stateVector
            Pp=SVobj.stateCov;
            %sizePp=size(Pp)
            %Ppinit=Pp
            %H(:,4:cH-rH)=[];  Pp(:,4:cH-rH)=[]; Pp(4:cH-rH,:)=[]; %original
            %Hinit=H
            H(:,4:cH)=0;  Pp(:,4:cH)=0; Pp(4:cH,:)=0;
            Hfinal=H
            Ppfinal=Pp
            %Debug
            %minOBStest=minNumObs
            %test=length(EDobj.svn)
            for i=1:length(EDobj.svn)
                sigmaResWoRcb(i,1)=sqrt(H(i,:)*Pp*H(i,:)'+ ...
                                         stdObsNoise^2+obsSISRE^2);  
            end
            %Debug
            sigmares=sigmaResWoRcb
            %Hfinal=H
            %matrisTest=H*Pp*H'
            %for i=1:length(EDobj.svn)
            %    testsigmaResWoRcb(i,1)=sqrt(H(i,:)*Pp*H(i,:)');  
            %end
            %    testsigmaResWoRcb
            %sigmaResWoRcb=sigmaResWoRcb
        % Search for the bad measurement
          outlierFlag=0;  remObs=[];
          while 1
            % Computed observable without receiver clock correction
            gpsClcBias;
            gpsTgd;
              [zcWoRcb]=geoRange-gpsClcBias+gpsTgd./2-ambBias./2 ;
            % Estimated standart deviation of receiver clock bias
            %Debug sigmaResWoRcb unrecognised function or variable
            %test=sigmaResWoRcb
            %test2=sigmaResWoRcb.^2
            %test3=sqrt(1/(sum(1./(sigmaResWoRcb.^2))))
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
            %Debug
            %display(stdFactor)
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
                  subSigmaRcb(min_ind)=[]; ambBias(min_ind)=[];
                % Update the variables for the next loop
                  remObs=[remObs,z(min_ind)]; % observations that will be removed
                  geoRange(min_ind)=[]; z(min_ind)=[]; 
                  gpsClcBias(min_ind)=[];
                  gpsTgd(min_ind)=[];
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
                   remInd=logical(remInd+(remObs(i)==EDobj.obs.Graphic));
               end
               remSatId=EDobj.svn(remInd);
               %Debug
               display(remSatId)
             % Update the epoch data, state vector and covariance
               [updSVobj,updEDobj]=this.removeObsFromDataSet( ...
                                                     SVobj,EDobj,remSatId);  
               %fivesixsixSVN=updEDobj.svn
          else
              %Debug
              %fivesixsevenSVN=EDobj.svn
             updEDobj=OrbEpochData(EDobj.obs,EDobj.svn,EDobj.gpsPos, ...
                                   EDobj.gpsClcCorr,EDobj.epochDateTime,...
                                   EDobj.hFlag,EDobj.gpsTgd);               
             updSVobj=OrbStateVector(SVobj.position,SVobj.velocity, ...
                                     SVobj.atmDragCoef,SVobj.solarRadCoef, ...
                                     SVobj.empAccel,SVobj.corelTime,...
                                     SVobj.recClcBias,SVobj.ambBias,...
                                     SVobj.stateDateTime);
             updSVobj.stateCov=SVobj.stateCov; 
          end
        % Set the a priory estimated receiver clock bias parameter
          % Find the index of receiver clock bias
            ind=updSVobj.stateVector==updSVobj.recClcBias;
          % Set the new receiver clock bias 
            updSVobj.recClcBias=recClcBias;     
            updSVobj.stateVector(ind)=recClcBias;
            updSVobj.stateCov(ind,:)=zeros;updSVobj.stateCov(:,ind)=zeros;
            updSVobj.stateCov(ind,ind)=sigmaRcb^2;
        % Control the observation count
        %Debug
        %fiveeightsevenEndSVN=length(updEDobj.svn)
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
             %Debug
             m=length(SVobj.stateVector);
             n=length(SVobj.ambBias);
             H = zeros(n,m);
             %positiare=SVobj.position 
             %velocitiare=SVobj.velocity
             %ATMDragCoeff=SVobj.atmDragCoef
             %SolarRadCoeff=SVobj.solarRadCoef
             %EmpAccel=SVobj.empAccel
             %CorelTime=SVobj.corelTime
             %ClcBias=SVobj.recClcBias
             %ambBiass=SVobj.ambBias
            for i = 1:n
                H(i,1) =-(GPS_pos(i,1)-xp(1))/geoRange(i);
                H(i,2) =-(GPS_pos(i,2)-xp(2))/geoRange(i);
                H(i,3) =-(GPS_pos(i,3)-xp(3))/geoRange(i);
                H(i,m-n) = 1;
                H(i,m-n+i) = -0.5;
            end 
            %display(H)
        end
        
        function [zc]=getCompObs(this,EDobj,SVobj)
        % FUNCTION
        %   Function computes the observations based on the measurement model 
        % INPUTS
        %   EDobj : instance of 'OrbEpochData' class including observations
        %           GPS ephemerides and other auxiliary data
        %   SVobj : instance of 'OrbStateVector' class including state
        %           vector variables  
        %   useBody2GPSTrans :
        % OUTPUTS
        %   zc    : Computed observations
        %            
        %----------------------------------------------------------------- 
        
        % Compute geometric range
            GPS_pos=EDobj.gpsPos; 
            GPS_tgd=EDobj.gpsTgd;
            xp=SVobj.position;
            xp=xp+this.bodyToGps; %body system to GPS antenna system
            geoRange = sqrt((GPS_pos(:,1)-xp(1)).^2+ ...
                            (GPS_pos(:,2)-xp(2)).^2+ ...
                            (GPS_pos(:,3)-xp(3)).^2);
        % Computed GRAPHIC observable
          RECclc_corr=SVobj.recClcBias; % Receiver clock bias
          GPSclc_corr=EDobj.gpsClcCorr; % Satellite clock bias
          B=SVobj.ambBias;
          zc=geoRange+RECclc_corr-GPSclc_corr+GPS_tgd./2-B./2;              
        end
        
        function [updSVobj,updEDobj]=removeObsFromDataSet(this,SVobj,EDobj,svn)
        % FUNCTION
        %   Function removes the observation given by satellite vehicle
        %   number from the data set
        % INPUTS
        %   EDobj : instance of 'OrbEpochData' class including observations
        %           GPS ephemerides and other auxiliary data
        %   SVobj : instance of 'OrbStateVector' class including state
        %           vector variables       
        % OUTPUTS
        %   updSVobj  : instance of 'OrbStateVector' class including updated 
        %               state vector parameters and covariance matrix           
        %   updEDobj : instance of 'OrbEpochData' class including updated 
        %              observations, GPS ephemerides and other auxiliary data
        %-----------------------------------------------------------------             

        % Find the index of observations that will be removed
          remInd=logical(0);
          for i=1:length(svn)
              remInd=logical(remInd+(svn(i)==EDobj.svn));
          end
          %Debug
          display(remInd)
        % Update the epoch data
          updEDobj=OrbEpochData(EDobj.obs,EDobj.svn,EDobj.gpsPos, ...
                                EDobj.gpsClcCorr,EDobj.epochDateTime,...
                                EDobj.hFlag,EDobj.gpsTgd); 
          updEDobj.obs.C1(remInd)=[];
          updEDobj.obs.L1(remInd)=[];
          updEDobj.obs.Graphic(remInd)=[];
          updEDobj.svn(remInd)=[];
          updEDobj.gpsPos(remInd,:)=[];
          updEDobj.gpsClcCorr(remInd)=[];
          updEDobj.gpsTgd(remInd)=[];
        % Update the state vector parameters
          % Size of parameters
            n_SV=length(SVobj.stateVector); % length of state vecctor
            n_Amb=length(SVobj.ambBias);    % length of ambiguity bias
          % Update state vector
            updAmbBias=SVobj.ambBias; 
            updAmbBias(remInd,:)=[];
          % Update state covariance matrix
            remId=find(remInd);  
            remId=(n_SV-n_Amb)+remId;
            Pupd=SVobj.stateCov; Pupd(remId,:)=[];Pupd(:,remId)=[];
        % Set output
          updSVobj=OrbStateVector(SVobj.position,SVobj.velocity, ...
                                  SVobj.atmDragCoef,SVobj.solarRadCoef, ...
                                  SVobj.empAccel,SVobj.corelTime,...
                                  SVobj.recClcBias,updAmbBias,...
                                  SVobj.stateDateTime);
          updSVobj.stateCov=Pupd;                       
        end
        
    
        
       
    end
    
end

