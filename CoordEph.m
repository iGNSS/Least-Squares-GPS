classdef CoordEph
%CoordEph - Coordinate ephemeris container - cartesian coordinates and clock corrections
% Can contain coordinates from precise or broadcast ephemeris
%X-coordinate                   dX
%Y-coordinate                   dY
%Z-coordinate                   dZ
%Satellite clock correction     dDts
%Satellite nr                   iPRN
%Time                           dTime
% Group delay                   dTGD
%Written by Milan Horemuz, last modified 2005-01-31 by Johan Vium Andersson
    properties (SetAccess = public)
        dX,dY,dZ,dDts,iPRN,dTime DateTime,dTGD,jd
    end
   methods
      function this = CoordEph()
        this.dX=[];this.dY=[];this.dZ=[];this.dDts=[];this.iPRN=[];
        this.dTime=DateTime;
        this.dTGD=0; this.jd=DateTime;
      end
      function display(b)
        for i = 1 : size(b.iPRN,2)
            fprintf('Epoch \n');
            b.dTime(i)
            for j = 1 : size(b.iPRN,1)
                if b.iPRN(j,i) == 0
                    continue;
                end
                fprintf('Satellite number %d\n', b.iPRN(j,i));
                fprintf('X-coordinate [m]: %11.3f\n', b.dX(j,i));
                fprintf('Y-coordinate [m]: %11.3f\n', b.dY(j,i));
                fprintf('Z-coordinate [m]: %11.3f\n', b.dZ(j,i));
                fprintf('Clock correction: %g\n', b.dDts(j,i));
            end
        end
      end
      
        function a = get(b,par)
                 a=b.(par);
        end
        function [A, data] = ReadEphBroadcast(A,sFileName,dInterval)
  
        % This function reads broadcast ephemeris from a RINEX file
% and calculates the satellit position with the wanted intervall
%
% To call the function one writes:
%
%      [A, data] = ReadEphBro(A,sFileName,dInterval)
%
%   A - CoordEph instance 
%   data - EphData 
%   sFileName - filename of the RINEX file
%   dInterval - the time intervall
%
%
% The function is written by Johan Vium Andersson
% Last time updated 2005-02-01
            A=CoordEph(); % Constructor
             [header,data] = brEphRinexRead(sFileName);
            [n m]=size(data);
            dLastEpoch =[];
            dFirstEpoch = [];
         % Find the first and last Epoch 
            first_time_flag = 0;  
            for i=1:n
                for j=1:m
                    BB(1) = get(data(i,j) , 'iEpochYear');
                    %BBdata=get(data(i,j),'iEpochYear')
                    %BB(1)=BBdata(1)    
                    BB(2) = get(data(i,j) , 'iEpochMonth');                        
                    %BBdata = get(data(i,j) , 'iEpochMonth')                        
                    %BB(2)=BBdata(2)
                    BB(3) = get(data(i,j) , 'iEpochDay');               
                    %BBdata = get(data(i,j) , 'iEpochDay');               
                    %BB(3)=BBdata(3)
                    BB(4) = get(data(i,j) , 'iEpochHour');               
                    %BBdata = get(data(i,j) , 'iEpochHour');               
                    %BB(4)=BBdata(4)
                    BB(5) = get(data(i,j) , 'iEpochMinute');        
                    %BBdata = get(data(i,j) , 'iEpochMinute');        
                    %BB(5)=BBdata(5)
                    BB(6) = get(data(i,j) , 'iEpochSecond');    
                    %BBdata = get(data(i,j) , 'iEpochSecond');    
                    %BB(6)=BBdata(6)
                    if abs(BB(1)) < 1e-13 % no data for this satellite
                        continue;
                    end
                    dJD = DateTime(BB); % Convert the time to JD
                    %if i==1 & j==1  %ORIGINAL
                    if first_time_flag ==0
                        dFirstEpoch = dJD;
                        dLastEpoch = dJD;
                        first_time_flag=1;
                    elseif (dJD - dFirstEpoch) < 0   
                        dFirstEpoch = dJD;
                    elseif (dJD - dLastEpoch) > 0  
                        dLastEpoch = dJD;
                    end
                end
            end    
        % Define start of the interval
            dStart = dFirstEpoch - 3*3600;
        % Define the end of the interval
            dEnd = dLastEpoch + 3*3600;
        %Time span in [s]
            TimeSpan = dEnd - dStart;
            %display(TimeSpan);
        %Number of intervals
            Nint = floor(TimeSpan/dInterval);
            data(1,1) = set(data(1,1), 'dtStart', dStart);
            data(1,1) = set(data(1,1), 'dtEnd', dEnd);
        % Calculate satellite coordinates for each interval
            epo = dStart; % counter
            for i = 1:Nint
                A.dTime(i) = epo;
                % Here should the resulult be placed in a CoordEph class
                vXYZ = beEphCoordCalc(data, epo);
                disp(vXYZ);
                [iRow,iCol] = size(vXYZ);       
                for p=1:iRow
                        prn = p;
                        A.iPRN(prn,i) = prn;
                        A.dX(prn,i)   =  vXYZ(p,1);
                        A.dY(prn,i)   =  vXYZ(p,2);
                        A.dZ(prn,i)   =  vXYZ(p,3);
                        A.dDts(prn,i) =  vXYZ(p,4);
                        A.dTGD(prn,i) =  vXYZ(p,5);
                end
                epo = epo +  dInterval;
            end
        end
        function A = ReadEphPrecise(varargin)
        %function A = ReadEphPrecise(A,sFileName_1 ... sFileName_n)

        %ReadEphPrecise reads precise ephemeris (coordinates and clock corrections) from sp3
        %file
        %A = ReadEphPrecise(A,sFileName) 
        %A instance of class CoordEph, which is container of coordinate ephemeris
        %sFileName 1 ... n - file(s) containing sp3 ephemeris
        %A is to be used in CompStdOrb

        %Written by Milan Horemuz, last modified 2005-02-03

        A = CoordEph; %create an instance
        i = 1; %counter of epochs for which the satellite coordinates are listed in PE
        for npe = 2:nargin
            %read header
            first ='1'; % variable for testing end of  header
            fid=fopen(varargin{npe}); %open file
            if fid < 1
                fprintf('Could not open file %s', sFileName);
                return;
            end
            tline = fgetl(fid);
            tline(1:3) = [];  %delete the first 3 characters 
            epoch = sscanf(tline,'%f');
            Time0 = DateTime(epoch);
            tline = fgetl(fid); %skip 2 tline
            tline = fgetl(fid);
            tline(1) = [];  %delete + character
            MofSat = sscanf(tline, '%i'); %read number of satellites
            while first ~= '*'
                tline = fgetl(fid);
                if length(tline) < 3
                    fprintf('Not a valid sp3 file');
                    return;
                end
                first = tline(1);
            end
            first ='1'; % variable for testing the line with time
            while feof(fid) == 0 & tline(1) ~= 'E';
                tline(1) = [];  %delete * character
                epoch = sscanf(tline,'%f');
                A.dTime(i) = DateTime(epoch(1), epoch(2), epoch(3), epoch(4), epoch(5), epoch(6));  %date2jd(epoch(1), epoch(2), epoch(3), epoch(4), epoch(5), epoch(6));
                tline = fgetl(fid);
                while first ~= '*' 
                    if tline(2) == 'R' %read in only GPS satellites
                        tline = fgetl(fid);
                        first = tline(1);
                        if first == 'E'
                            break;
                        end
                        continue;
                    end
                    tline(1:2) = []; %delete PG
                    prn = sscanf(tline(1:2), '%d'); %satellite number
                    epoch = sscanf(tline,'%f'); 
                    A.iPRN(prn,i) = epoch(1);  
                    A.dX(prn,i) = epoch(2)*1000; %convert to [m]
                    A.dY(prn,i) = epoch(3)*1000; 
                    A.dZ(prn,i) = epoch(4)*1000; 
                    A.dDts(prn,i) = epoch(5)/1e6; %convert to [s]
                    tline = fgetl(fid);
                    if tline(1) == 'E'
                        break;
                    end
                    first = tline(1);
                end
                i = i+1;
                first ='1';
            end
            fclose(fid);
        end
        end
        
        function vpeEph = RemoveEpoch(vpeEphi,tRef)

            %function vpeEph = RemoveEpoch(vpeEphi,tRef)
            %Removes epoch tRef (DateTime) from vpeEphi (CoordEph)
            %Returns vpeEph (CoordEph) with removed epoch

            %Written by Milan Horemuz, last modified 2005-02-04


            vpeEph = vpeEphi;
            t = get(vpeEph,'dTime');  %extract vector of eopchs
            nep = length(t);  %number of epochs in vpeEph
            epnr = -1;
            for i = 1:nep
                if abs(t(i) - tRef) < 1e-5
                    epnr = i;
                    break;
                end
            end
            if epnr < 0
                tRef
                error('Could not find this epoch in vpeEph');
            end

            vpeEph.dX(:,epnr) = [];
            vpeEph.dY(:,epnr) = [];
            vpeEph.dZ(:,epnr) = [];
            vpeEph.dDts(:,epnr) = [];
            vpeEph.iPRN(:,epnr) = [];
            vpeEph.dTime(:,epnr) = [];
        end
        function b = set(b,par,in)
            b.(par)=in;
        end
      end
        


end
