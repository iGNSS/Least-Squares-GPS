classdef OrbMeasurementModel < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        settings
        obsFileObj
        brdEphObj
        timeRefSysObj
    end
    
    methods
        function this=OrbMeasurementModel()
        %-----------------------------------------------------------------            
        % FUNCTION
        %   Contructure function for the class
        %--------------------------------------------------------------  

        end
                    
        function ModelObj=getInstance(this,varargin)
        %-----------------------------------------------------------------            
        % FUNCTION
        %   Creates object related with mesurement models; Navigation
        %   Solution, Code Pseudorange or GRAPHIC. 
        % INPUTS
        %   varargin : cell structure including input parameters. The class
        %              can be initiated as follows;
        %
        %              For Code and GRAPHIC measurements
        %                 obj=OrbMeasurementModel(filtSet,obsFileObj,brdEphObj)
        %              For navigation solution measurements
        %                 obj=OrbMeasurementModel(filtSet,obsFileObj)    
        %              
        %              where
        %                filtSet: data structure including filter settings
        %                obsFileObj: Handle to instance of 'OrbObsFile' class
        %                brdEphObj: Handle to instance of 'BrdEphObj' class   
        % OUTPUTS 
        %   ModelObj : Handle to child class instance
        %--------------------------------------------------------------   
        
        % Get the filter setings
          settings=varargin{1};
        % Create Measurement Objects
          if strcmp(settings.obsType,'NavSol')
             ModelObj=OrbNavSolMeasModel();
          elseif strcmp(settings.obsType,'Graphic')
             ModelObj=OrbGraphicMeasModel();
          elseif strcmp(settings.obsType,'Code')
             ModelObj=OrbCodeMeasModel();
          else
             error ('Input measurement model type can not be recognised')
          end
        % Set the parent class variables
          n=length(varargin);
          ModelObj.settings=varargin{1};
          ModelObj.obsFileObj=varargin{2};          
          if n==3
            ModelObj.brdEphObj=varargin{3};
          end  
            
        end
        
    end
    
end

