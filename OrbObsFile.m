classdef OrbObsFile < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        obsFileReaderObj
        fileType
        filtSet
        current_t
    end
    
    methods
        
        function this=OrbObsFile(filtSet,fileType)
        % FUNCTION
        %   Creates an object from 'OrbObsFile' class for a given specific 
        %   observation file, and reads the relevant information from the file
        % INPUTS
        %   filtSet : filter setting
        %   fileType: type of the observation file. It takes 'rinex' or 'nav_file' 
        %----------------------------------------------------------------      
        this.fileType = fileType;
        this.filtSet  = filtSet;
        % Check the file type 
          if ~strcmp(fileType,'rinex') && ~strcmp(fileType,'nav_file')
               error('File type is not supported')
          end
        % open the initial observation file
          this.loadObsFile(filtSet.beginDateTime)
          
          %this.obsFileReaderObj=rinexObsReader2(filePath);
        end
        
        function [orbEpochDataObj,eof_flag]= getNextObs(this,obsType)
        % FUNCTION
        %   Reads the next observation from the file and returns an object
        %   of 'OrbEpochData' class including all relevant data 
        % INPUTS
        %   obsType: Type of observation. It can take 'Code' or 'Graphic' or 'NavSol'
        % OUTPUTS
        %   orbEpochDataObj: Instance of "orbEpochData" class
        %----------------------------------------------------------------                
        
        % Initialize the variables
          eof_flag=0;
          c = 299792458;          % speed of light in vacum (m/sec)
          lambda_L1=(c/1575.42e6); % wavelength of L1 carrier

        % Get the Data
          while eof_flag==0
            if strcmp(obsType,'Code')
              [health_flag,date_time,GPSobs,sat_vn,sat_id,eof_flag ]= ...
                                this.obsFileReaderObj.getNextObs({'C1'},{'G'});
            elseif strcmp(obsType,'Graphic')
              [health_flag,date_time,GPSobs,sat_vn,sat_id,eof_flag ]= ...
                                this.obsFileReaderObj.getNextObs({'C1','L1'},{'G'}); 
            elseif strcmp(obsType,'NavSol')
              sat_vn=[];
              [health_flag,date_time,GPSobs,eof_flag ]= ...
                                this.obsFileReaderObj.getNextObs();
            end
          % Control the end of file
            if eof_flag==0
               break;
            elseif eof_flag==1
               % Read the next file
                 [file_exist_flag]=this.loadObsFile(this.current_t + 86400);
                 if file_exist_flag == 1
                    eof_flag = 0;
                 else
                    orbEpochDataObj=[];   
                    return;
                 end
            end
          end              
        % Edit the Data
          if strcmp(obsType,'C1') || strcmp(obsType,'Graphic')
             % Sort the data
               [sat_vn,ind]=sort(sat_vn);
                GPSobs=GPSobs(ind,:);
                sat_id=sat_id(ind);
             % Remove the repeated observations
               if length(unique(sat_vn))~=length(sat_vn)
                  rep_sv_ind=(diff(sat_vn )== 0);
                  sat_vn(rep_sv_ind)=[];
                  GPSobs(rep_sv_ind,:)=[];
                  sat_id(rep_sv_ind)=[];  
               end
          end
        % Set the observations to output
          if strcmp(obsType,'Code') 
             obs.C1=GPSobs(:,1); 
          elseif strcmp(obsType,'Graphic')
              obs.C1=GPSobs(:,1);
              obs.L1=GPSobs(:,2);
              obs.Graphic=(GPSobs(:,1)+GPSobs(:,2).*lambda_L1)./2;
          elseif strcmp(obsType,'NavSol')
              obs.NavSol=GPSobs;
          end
          orbEpochDataObj=OrbEpochData(obs,sat_vn,[],[],date_time,health_flag,[]);              
        end
    
        function [file_exist_flag]=loadObsFile(this,t)
        % FUNCTION
        %   Loads the next observation file 
        %---------------------------------------------------------------- 
        file_exist_flag = 1;

        % file name
          doy = dayofyear(t);
          if strcmp(this.fileType,'rinex')
              file_path = strcat(this.filtSet.ws,'inputs/obs/', ...
                           num2str(get(this.filtSet.beginDateTime,'year')), ...
                           '_',num2str(doy),'.rnx');
             
          elseif strcmp(this.fileType,'nav_file')
              file_path = strcat(this.filtSet.ws,'inputs/obs/', ...
                           num2str(get(this.filtSet.beginDateTime,'year')), ...
                           '_',num2str(doy),'.nav');
           
          else
             error('File type are not supported') 
          end          
        % Open the observation file
          if strcmp(this.fileType,'rinex')
              try
                this.obsFileReaderObj=rinexObsReader2(file_path);
                this.current_t = t;
              catch
                file_exist_flag = 0;
              end
          elseif strcmp(this.fileType,'nav_file')
              try              
                this.obsFileReaderObj=navFileObsReader(file_path);
                this.current_t = t;
              catch
                file_exist_flag = 0;

              end
          else
             error('File type are not supported') 
          end
          
        end
    end
    
end

