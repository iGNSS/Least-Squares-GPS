classdef OrbFilterOutputFile < handle
    
    
    properties
        fileID
        filePathforTextFile
        filePathforMatFile       
        filtSet
        dataFormat
        epochCounter=0;
    end
    
    methods
        function this=OrbFilterOutputFile(filtSet)
            this.filtSet=filtSet;
        end        
        function setOutputTextFile(this,date_time_array)
        % FUNCTION
        %   Set the output file path for text document including filter
        %   navigation outputs
        % INPUTS
        %   out_fil_path : file path to output file
        %---------------------------------------------------------------- 
        outTextFilePath = strcat(this.filtSet.ws,'outputs/', ...
                                 sprintf('%04d%02d%02d',date_time_array(1), ...
                                 date_time_array(2),date_time_array(3)),'.txt');
        if ~isequal(this.filePathforTextFile,outTextFilePath) || isempty(this.fileID)
            % Close the previous file
              if ~isempty(this.fileID)
                fclose(this.fileID);
              end
            % Create the output file
              this.fileID=fopen(outTextFilePath,'wt+');   
              this.filePathforTextFile = outTextFilePath;
            % Write file header
              % Write Filter Settings
                fprintf(this.fileID,'FILTER SETTINGS\n');
                fprintf(this.fileID,'Filter Type     : %s \n',this.filtSet.filterType);
                fprintf(this.fileID,'Observation Type: %s \n',this.filtSet.obsType);
                fprintf(this.fileID,'\n');
              % Write Filter Output Data Format
                fprintf(this.fileID,'FILTER DATA FORMAT\n'); 
                fprintf(this.fileID,'[ year month day hour minute second ] [fFlag]\n');            
                fprintf(this.fileID,'[ x y z vx vy vz ');
                if this.filtSet.enableDynModParam==1 || ...
                   this.filtSet.enableDynModParam==2  
                   fprintf(this.fileID,'Cdrag Csrad aR aT aN '); 
                    if this.filtSet.enableDynModParam==2 
                       fprintf(this.fileID,'Ctime ');    
                    end
                end
                if strcmp(this.filtSet.obsType,'Graphic') || ...
                   strcmp(this.filtSet.obsType,'Code')
                   fprintf(this.fileID,'RecClc '); 
                end
                fprintf(this.fileID,']\n');
                fprintf(this.fileID,'      Where;\n'); 
                fprintf(this.fileID,'         x,y,z   : positions   \n');   
                fprintf(this.fileID,'         vx,vy,vz: velocities   \n');              
                if this.filtSet.enableDynModParam==1 || ...
                   this.filtSet.enableDynModParam==2  
                    fprintf(this.fileID,'         Cdrad   : atmospheric drag coefficient   \n');  
                    fprintf(this.fileID,'         Csrad   : solar radiation coefficient   \n');  
                    fprintf(this.fileID,'         aR,aT,aN: empirical accelerations  \n');
                    fprintf(this.fileID,'                   radial, tangential and normal directions, respectively  \n');  

                   if this.filtSet.enableDynModParam==2 
                    fprintf(this.fileID,'         Ctime   : Markov process corelletion time  \n');  
                   end
                end  
                if strcmp(this.filtSet.obsType,'Graphic') || ...
                   strcmp(this.filtSet.obsType,'Code')
                    fprintf(this.fileID,'         RecClc  : Receiver clock   \n');  
                end
                fprintf(this.fileID,'         fFlag    : Filter prediction-update flag. \n');              
                fprintf(this.fileID,'                    1 for only time update, 0 for measurement update   \n');              
                fprintf(this.fileID,'\n');
              % Write Filter Outputs
                fprintf(this.fileID,'FILTER OUTPUTS\n');
              % Set the filter data format string
                % Set the format string for the position and velocity
                  dFormat=['%15.4f',' ','%15.4f',' ','%15.4f',' ','%15.6f',' ','%15.6f',' ','%15.6f'];
                  len=6;
                % Set the format string for dynamical model parameters
                  if this.filtSet.enableDynModParam==1 || this.filtSet.enableDynModParam==2  
                     dFormat=[dFormat,' ','%15.4f',' ','%15.4f',' ','%.15e',' ','%.15e',' ','%.15e'];
                     len=len+5;
                     if this.filtSet.enableDynModParam==2
                        dFormat=[dFormat,' ','%15.4f'];
                        len=len+1;
                     end
                  end
                % Set the format string for the clock estimation
                  if strcmp(this.filtSet.obsType,'Graphic') || ...
                     strcmp(this.filtSet.obsType,'Code')            
                     dFormat=[dFormat,' ','%15.6f'];
                     len=len+1;
                  end
                % Set the class variable
                  this.dataFormat.formatString=dFormat;
                  this.dataFormat.length=len;
        end
        end
        
        function setOutputMatFile(this,date_time_array)
        % FUNCTION
        %   Contructer function for 'OrbFilterOutputFile' class
        % INPUTS
        %   out_fil_path : file path to output file
        %---------------------------------------------------------------- 
        outMatFilePath = strcat(this.filtSet.ws,'outputs/', ...
                                 sprintf('%04d%02d%02d',date_time_array(1), ...
                                 date_time_array(2),date_time_array(3)),'.mat');
        % Set the output file path                                               
          if ~isequal(this.filePathforMatFile,outMatFilePath) || isempty(this.fileID)
             this.filePathforMatFile=outMatFilePath;
             this.epochCounter=0;
          end
               
        end
           
        function saveStateDataIntoTextFile(this,SVobj)
        % FUNCTION
        %   save the filte navigation outputs into the text file
        % INPUTS
        %   SVobj : instance of 'OrbStateVector' class including state vector 
        %           parameters taht will be stored in an output file  
        %---------------------------------------------------------------- 
        % Check the file
          this.setOutputTextFile(SVobj.stateDateTime)
        % Save the time of filter output
          fprintf(this.fileID,'%d %d %d %d %d %15.13f ',SVobj.stateDateTime'); 
          fprintf(this.fileID,' %d\n',SVobj.updateFlag);      
        % Save the state vector
          % Length of state vector that will be stored
            sLen=this.dataFormat.length;
          % state vector that will be stored
            sVector=SVobj.stateVector(1:sLen)';
          % Store the state vector parameters  
            fprintf(this.fileID,this.dataFormat.formatString,sVector); 
            fprintf(this.fileID,'\n');
        end
        
        function saveFilterDataIntoMatFile(this,SVobj,Edobj)
        % FUNCTION
        %   save the given state and epoch data in binary format into the
        %   mat file.
        % INPUTS
        %   SVobj : instance of 'OrbStateVector' class including state vector 
        %           parameters taht will be stored in an output file  
        %   Edobj : instance of 'OrbEpochData' class including epoch data 
        %           that will be stored in the output mat file  
        %---------------------------------------------------------------- 
        
        % Check the file
         this.setOutputMatFile(SVobj.stateDateTime)
        % Save the data  
         this.epochCounter=this.epochCounter+1;
         data_string_1=['epoch_data_',num2str(this.epochCounter),'.SVobj','=SVobj;'];
         data_string_2=['epoch_data_',num2str(this.epochCounter),'.Edobj','=Edobj;'];
         eval(data_string_1); eval(data_string_2);
         if this.epochCounter==1
            save(this.filePathforMatFile,['epoch_data_',num2str(this.epochCounter)]);
         else
            save(this.filePathforMatFile,['epoch_data_',num2str(this.epochCounter)],'-append');
         end

        end
        

    end
    
end

