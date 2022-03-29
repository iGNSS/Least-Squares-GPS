classdef rinexObsReader2

    
    properties 
        rinex_file_id              % file identifier
        num_obs_type               % number of observation types
        obs_type;                  % observation types in the file
        antenna_height_HEN=[0 0 0] % default value for the antenna height
        app_pos_XYZ   =[0 0 0 ]    % default value for the approximate position 
                                   % of the site 
        interval=30;               % default value for the  file sampling time
        tof;                       % time of first observation
        file_data_pointer          
    end
    properties (Access=public)
 %       obsIndex;
    end
    methods
        % Constructer
        function this=rinexObsReader2(file_path)
        % FUNCTION
        %   Contructure function loads the file and read the metadata in
        %   the header
        % INPUTS
        %  file_path : File path to the rinex file 
        % OUTOUTS
        %
        %----------------------------------------------------------------  
        
        % Load the file
          this.rinex_file_id=fopen(file_path,'r');
          if this.rinex_file_id==-1
             error('File can not be opened')
          end
        % Read header file
          while 1
                strline=fgetl(this.rinex_file_id);
                if ~isempty(strfind (strline,'RINEX VERSION / TYPE'))
                    ver=sscanf(strline,'%f');
                    if (ver== 2.11) || (ver== 2.10) || (ver== 2)
                        % continue
                    else
                        error ('Version of rinex file are not supported')
                    end
                if isempty(strfind (strline,'OBSERVATION DATA'))
                    error ('Only observation files can be red')
                end
                
                elseif ~isempty(strfind (strline,'APPROX POSITION XYZ'))
                    this.app_pos_XYZ=sscanf(strline,'%f %f %f'); 
                    
                    
                elseif ~isempty(strfind (strline,'ANTENNA: DELTA H/E/N'))
                    this.antenna_height_HEN=sscanf(strline,'%f %f %f');                       
                    
                elseif ~isempty(strfind (strline,'# / TYPES OF OBSERV'))
                       % Read the umber of observations
                        remain=strline(1:60);
                        [token,remain]=strtok(remain); 
                        this.num_obs_type=str2num(token);
                       % Read the observation types
                         obs_type={};
                         while 1 
                               % Read observation type line
                                 obs_type_temp=textscan(remain,'%s');
                                 obs_type_temp=obs_type_temp{:};
                                 obs_type={obs_type{:},obs_type_temp{:}};
                               % Control whether exist another line
                                 if length(obs_type)==this.num_obs_type
                                     this.obs_type=obs_type;
                                     break
                                 else
                                     % Get the next line
                                       strline=fgetl(this.rinex_file_id);
                                       remain=strline(1:60);
                                 end
                          end
                elseif ~isempty(strfind (strline,'INTERVAL'))
                        this.interval=sscanf(strline,'%d'); 
                        
                elseif ~isempty(strfind (strline,'TIME OF FIRST OBS'))
                        this.tof=sscanf(strline,'%f %f %f %f %f %f'); 
                        
                elseif ~isempty(strfind (strline,'END OF HEADER'))
                        break       
                end
          end
        % Set the file data pointer
          this.file_data_pointer=ftell(this.rinex_file_id);

        end
    end
    methods
        function [ flag,date_time,obs,sat_vn,sat_id,eof_flag ]= getNextObs(...
                                                             this,obs_name,sat_type)
        %------------------------------------------------------------------            
        %FUNCTION 
        %   return desired observable at current epoch
        %     [ flag,date_time,obs,sat_vn,sat_id,eof_flag ] = getNextObs(this,obsType )
        %   Example
        %     [ flag,date_time,obs,sat_vn,sat_id,eof_flag ] = getNextObs({'L1','L2','C1'})
        %INPUTS
        %   obs_name :cell array for name of required observables
        %               such as; obsType={'L1','L2'}
        %   sat_type :cell array for the satellite type that will be
        %             returned
        %               such as; obsType={'G','R'}
        %            It can take;
        %                G          : GPS
        %                R          : GLONASS
        %                S          : Geostationary signal payload
        %                E          : Galileo              
        %OUTPUTS
        %   flag:  Epoch flag, 0: OK 
        %                      1: power failure between previous and current epoch 
        %                      >1: Event flag    
        %   time: [year,month,day,hour,minute,second]
        %   obs : requested observation
        %   sat_vn : PRN number of satellite
        %   sat_id : cell including PRN code of the satellite
        %               such as sv={id1;id2;id3}
        %            It takes;
        %                G          : GPS
        %                R          : GLONASS
        %                S          : Geostationary signal payload
        %                E          : Galileo        
        %   eof_flag: End of file flag for rinex file
        %------------------------------------------------------------------
        
        % read epoch identifier line
          strline=fgetl(this.rinex_file_id);  
          eof_flag=0;
        % Extract data
          if ~isempty(strfind (strline,'END')) || isequal(feof(this.rinex_file_id),1) 
             eof_flag=1; 
             flag=false; date_time=false; obs=false; sat_vn=false; 
             sat_id=false;
          else
             % Read the time and satellite vehicle numbers

               date_time=sscanf(strline(1:28),'%d %d %d %d %d %f')';
               if date_time(1)<=79
                   date_time(1)=date_time(1)+2000;
               end
             % Read the time health flag
               flag=sscanf(strline(27:29),'%d');
             % Nunmer of satellites
               num_of_sat=str2double(strline(30:32));
             % Read observed satellite vehicle numbers
               if num_of_sat>12
                  str_sat=strline(33:33+12*3-1);
               else
                  str_sat=strline(33:33+num_of_sat*3-1);  
               end
               %debug str_sat
               %if num_of_sat==0
                %    str_sat={}
               %end
               %debug end
                   
               sat_id={};sat_vn=[]; i=1;data_str=[];
               if num_of_sat>0
               while 1
                   sat_id={sat_id{:},str_sat(i)};
                   sat_vn=[sat_vn(:);str2double(str_sat(i+1:i+2))];
                   i=i+3;
                   if length(sat_vn) == num_of_sat
                       break
                   elseif length(sat_vn) == 12
                       strline=fgetl(this.rinex_file_id);
                       str_sat=strline(33:end);
                       i=1;
                   end
               end
               % If satellite identifier number is blank change it with
               % char 'G'y
                 for i=1:num_of_sat
                     if sat_id{i}==' ';
                        sat_id{i}='G';
                     end
                 end
                 sat_id=sat_id';
            % Get Epoch Data       
              for i=1:num_of_sat
                  % Get the data string
                    k=0; data_str=[];
                    while k<this.num_obs_type 
                        line_str=fgetl(this.rinex_file_id);
                        if length(line_str)<80
                           line_str(end:80)=char(0); 
                        end
                        data_str=[data_str,line_str];
                        k=k+5;
                    end
              %end
               %end
                  % Convert to numeric array
                    for kk=1: this.num_obs_type
                        if isempty(data_str)
                            num_val=[];
                        else
                        num_val=str2num(data_str(16*(k-1)+1:16*k-2));
                        end
                        if isempty(num_val)
                           obs(i,kk) = nan;
                        else
                           obs(i,kk) = num_val;
                        end
                         
                    end
                    
              end
               end

              end
            % Select the requested data
              % Find the index of requested data
              %uzunluk=length(obs_name)
              %tip=this.obs_type;
                for i=1:length(obs_name)
                    ind(i,:)=logical(strcmp(obs_name{i},this.obs_type)); 
                end
              % Set the requested data 
                for i=1:length(obs_name)
                    temp_obs(:,i)=obs(:,ind(i,:));   
                end
                obs=temp_obs;
            % Select the requested satellites
              ind=logical(0);
              for i=1:length(sat_type)
                  ind=logical(ind+strcmp(sat_id,sat_type{i}));
              end
              obs=obs(ind,:);
              sat_vn=sat_vn(ind);
              sat_id=sat_id(ind);
              %sat_id or str_sat empty debug
              size(sat_vn);
              size(sat_id);
              
              end
         
       
        function closeFile (this)
            fclose(this.rinex_file_id);
        end
    end
    
end

