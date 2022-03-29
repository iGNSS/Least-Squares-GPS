classdef OrbEpochData < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        obs
        svn
        epochDateTime
        gpsPos
        gpsClcCorr
        hFlag
        gpsTgd
        
    end
    
    methods
        function this=OrbEpochData(obs,svn,gpsPos,gpsClcCorr,epochDateTime,hFlag,gpsTgd)
        % FUNCTION
        %   Function stores the observation data and related auxiliary data 
        %   acquired at a specific time of epoch 
        % INPUTS
        %   obs           : struct including observations, such as;
        %                      obs.C1 , obs.L1 , obs.Graphic, obs.NavSol
        %   svn           : Array including vehicle number of the GPS satellites  
        %   gpsPos        : array including position of tracked GPS
        %                   satellites
        %   gpsClcCorr    : array including clock correction of tracked GPS
        %                   satellites
        %   epochDateTime : Array including date and time, such as;
        %                    EpochDateTime=[year month day hour minute second]
        %   hFlag         : Takes the value "0", if the observations are
        %                   healthy
        %---------------------------------------------------------------- 
        
        % Store the data
          this.obs=obs;
          this.svn=svn;
          this.gpsPos=gpsPos;
          this.gpsClcCorr=gpsClcCorr;
          this.epochDateTime=epochDateTime;
          this.hFlag=hFlag;
          this.gpsTgd=gpsTgd;
            
        end
    end
    
end

