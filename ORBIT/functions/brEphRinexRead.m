function [Header,out] = brEphRinexRead(sFileName)

% This function reads the header and ephemeris from RINEX format
%
% The function gives/takes the following output/input parameters 
%
%    [Header,out] = brEphRinexRead(sFileName)
%  
%    sFileName - String with the filename
%    Header  - Return the header of the RINEX file
%    out - Returns the data from the RINEX file
%



fid=fopen(sFileName);

%% ----------------------------------------------------------
rad=0; % Start value Counter
num=0; % Start value Counter
tline= fgetl(fid);
countComments=0;
z=[]; % z is a vector with all the satellites in the file (see below)
A=EphHeader; % Creates a Header object

while num<1 % Loop until END OF HEADER
    
    matches = findstr(tline,'END OF HEADER');
    num = length(matches);
        
    if num==0; % If header is not on the line take next
        
        % RINEX VERSION /Type of observation
        matches = findstr(tline,'RINEX VERSION');
        num = length(matches);
        if num>0
            A = set(A,'iRinexVersion',str2double(tline(6:(6+14))));
            A = set(A,'sFileType',tline(21:40));
         
            num=0;
        end
        
        % PGM / RUN BY / DATE
        matches = findstr(tline,'PGM / RUN BY / DATE');
        num = length(matches);
        if num>0
            A =set(A,'sPGM',tline(1:20));
            A =set(A,'sRunBy',tline(21:40));
            A =set(A,'sDate',tline(41:60));
            num=0;    
        end
        
        % COMMENTS Comment line(s)  
        matches = findstr(tline,'COMMENT');
        num = length(matches);
        if num>0
            countComments=countComments+1;
            if countComments==1;
                sComment=tline(1:60);
            else 
                sComment=[A.aComment;tline(1:60)];
            end
            A =set(A,'sComment',sComment);
            num=0;
        end
        
        
        % ION ALPHA, Ionosphere parameters A0-A3 of almanac
        matches = findstr(tline,'ION ALPHA');

        num = length(matches);
        if num>0
            A=set(A,'dAlfaIon1',str2double(strrep(tline( 3:15),'D','e')));
            A=set(A,'dAlfaIon2',str2double(strrep(tline(16:28),'D','e')));
            A=set(A,'dAlfaIon3',str2double(strrep(tline(29:40),'D','e')));
            A=set(A,'dAlfaIon4',str2double(strrep(tline(41:54),'D','e')));
            num=0;    
        end
        
        % ION BETA, Ionosphere parameters B0-B3 of almanac
        matches = findstr(tline,'ION BETA');
        num = length(matches);
        if num>0
            A=set(A,'dBetaIon1',str2double(strrep(tline( 3:15),'D','e')));
            A=set(A,'dBetaIon2',str2double(strrep(tline(16:28),'D','e')));
            A=set(A,'dBetaIon3',str2double(strrep(tline(29:40),'D','e')));
            A=set(A,'dBetaIon4',str2double(strrep(tline(41:54),'D','e')));
            num=0;    
        end
        
        % DELTA-UTC, Almanac parameters to compute time in UTC
        matches = findstr(tline,'DELTA-UTC');
        num = length(matches);
        if num>0
           A=set(A,'dA0',str2double(strrep(tline( 4:23),'D','e')));
           A=set(A,'dA1',str2double(strrep(tline(24:43),'D','e')));
           A=set(A,'dRefTime',str2double(strrep(tline(44:53),'D','e')));
           A=set(A,'dRefWeek',str2double(strrep(tline(54:60),'D','e')));
           num=0;    
        end
       
        
        % LEAP SECONDS, Delta time due to leap seconds   
        matches = findstr(tline,'LEAP SECONDS');
        num = length(matches);
        if num>0
            
            A=set(A,'iLeapSeconds',str2double(strrep(tline(1:6),'D','e')));
            num=0; 
        end
        
    tline= fgetl(fid);
    rad=rad+1; % Counter

    end
   
end
Header=A;

%% DATA
%% ************************************************************************
   C = EphData;

while feof(fid) == 0; % Loop continues until the end of the file
     B=EphData;
     
    
    % PRN / EPOCH / SV CLK
        tline= fgetl(fid); rad=rad+1; % Take a new line from the file
   
        B=set(B,'iPRN',str2double(strrep(tline(1:2),'D','e')));
    
        year=str2double(tline(3:5));
            % Convert the year to 4 digets
            if year<94;
                year=year+2000;
            else
                year=year+1900;
            end       

        B=set(B,'iEpochYear',year);
        B=set(B,'iEpochMonth',str2double(strrep(tline(6:8),'D','e')));    
        B=set(B,'iEpochDay',str2double(strrep(tline(9:11),'D','e')));
        B=set(B,'iEpochHour',str2double(strrep(tline(12:14),'D','e')));
        B=set(B,'iEpochMinute',str2double(strrep(tline(15:18),'D','e')));        
        B=set(B,'iEpochSecond',str2double(strrep(tline(19:22),'D','e')));    
        B=set(B,'dClockBias',str2double(strrep(tline(23:41),'D','e')));
        B=set(B,'dClockDrift',str2double(strrep(tline(42:60),'D','e')));
        B=set(B,'dClockDriftRate',str2double(strrep(tline(61:79),'D','e')));       
    
    % Check if there are observation data, and if switch the boolean
    % operator to inform that there are data in the object
        
    % BROADCAST ORBIT - 1
        tline= fgetl(fid);rad=rad+1;
    
        B=set(B,'dIDOE',str2double(strrep(tline(4:22),'D','e')));    
        B=set(B,'dCrs',str2double(strrep(tline(23:41),'D','e')));
        B=set(B,'dDeltaN',str2double(strrep(tline(42:60),'D','e')));
        B=set(B,'dM0',str2double(strrep(tline(61:79),'D','e')));  
 
    % BROADCAST ORBIT - 2
        tline= fgetl(fid);rad=rad+1;
        
        B=set(B,'dCuc',str2double(strrep(tline(4:22),'D','e')));    
        B=set(B,'dEccent',str2double(strrep(tline(23:41),'D','e')));
        B=set(B,'dCus',str2double(strrep(tline(42:60),'D','e')));
        B=set(B,'dSqrtA',str2double(strrep(tline(61:79),'D','e')));
    
    % BROADCAST ORBIT - 3
        tline= fgetl(fid);rad=rad+1;
        
        B=set(B,'dToe',str2double(strrep(tline(4:22),'D','e')));    
        B=set(B,'dCic',str2double(strrep(tline(23:41),'D','e')));
        B=set(B,'dOMEGA',str2double(strrep(tline(42:60),'D','e')));
        B=set(B,'dCis',str2double(strrep(tline(61:79),'D','e')));
  
     % BROADCAST ORBIT - 4
        tline= fgetl(fid);rad=rad+1;
        
        B=set(B,'di0',str2double(strrep(tline(4:22),'D','e')));    
        B=set(B,'dCrc',str2double(strrep(tline(23:41),'D','e')));
        B=set(B,'dOmega',str2double(strrep(tline(42:60),'D','e')));
        B=set(B,'dOMEGADot',str2double(strrep(tline(61:79),'D','e')));
      
        
  % BROADCAST ORBIT - 5
        tline= fgetl(fid);rad=rad+1;
        
        B=set(B,'dIdot',str2double(strrep(tline(4:22),'D','e')));    
        B=set(B,'dCodeOnL2',str2double(strrep(tline(23:41),'D','e')));
        B=set(B,'dGpsWeek',str2double(strrep(tline(42:60),'D','e')));
        B=set(B,'dPDataFlag',str2double(strrep(tline(61:79),'D','e')));
 
 % BROADCAST ORBIT - 6
        tline= fgetl(fid);rad=rad+1;
     
        B=set(B,'dSVaccur',str2double(strrep(tline(4:22),'D','e')));    
        B=set(B,'dSVhealth',str2double(strrep(tline(23:41),'D','e')));
        B=set(B,'dTGD',str2double(strrep(tline(42:60),'D','e')));  
        B=set(B,'dIODC',str2double(strrep(tline(61:79),'D','e')));

  % BROADCAST ORBIT - 7
        tline= fgetl(fid);rad=rad+1;
        
        
        
        if length(tline)>78
            B=set(B,'dTransTime',str2double(strrep(tline(4:22),'D','e')));
            B=set(B,'dSpare1',str2double(strrep(tline(23:41),'D','e')));
            B=set(B,'dSpare2',str2double(strrep(tline(42:60),'D','e')));
            B=set(B,'dSpare3',str2double(strrep(tline(61:79),'D','e')));
           
        elseif length(tline)>55
            B=set(B,'dTransTime',str2double(strrep(tline(4:22),'D','e')));
            B=set(B,'dSpare1',str2double(strrep(tline(23:41),'D','e')));
            B=set(B,'dSpare2',str2double(strrep(tline(42:60),'D','e')));
           
        elseif length(tline)>23
            B=set(B,'dTransTime',str2double(strrep(tline(4:22),'D','e')));
            B=set(B,'dSpare1',str2double(strrep(tline(23:41),'D','e')));
           
        elseif length(tline)<23
            B=set(B,'dTransTime',str2double(strrep(tline(4:22),'D','e')));
           
        end

 % Routine: To prevent that data is over written when data from the same
 % satellite is archived twice
 
    % Start values
        xx=0;       % switch
        i=0;        % Reset counter
       
    % The actual satellite
        x=get(B,'iPRN');
       
        if length(z)==0 % First satellite
            z=x;            % z - new vector with the satellite numbers
            C(x,1) = B;       
            zCount=1;
            
        else % otherwise
        
            while ((xx==0) & (i<length(z))) 
            i=i+1; % Counter
              
              % check if the satellite data already is imported 

              if x==z(i)
                  xx=1;
              else
                  xx=0;
              end
              
            end
           
            % switch case - if data does not exist, enter it on the first position
            % otherwise in the second position
            
            switch xx
                case 0
                    z=[z,x];
                    zCount=[zCount,1];
                    C(x,1)=B;
                case 1
                    zCount(i)=zCount(i)+1;
                    C(x,zCount(i))=B;
                    
            end
            
        end
          
        clear B;

end

% Return the the header and all data
out = [C];