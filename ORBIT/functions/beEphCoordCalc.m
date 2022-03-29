function vCoordinateVector = beEphCoordCalc(mData,vObsTime)

% This function calculates the coordinates from the broadcast ephemeris at
% all aviliable satellits in the the inserted RINEX file
%
% The function takes the following input parameters 
%
%    vCoordinateVector = beEphCoordCalc(sFileName,vObsTime)
%  
%    sFileName - String with the filename
%    vObsTime  - time in [DateTime]
%
%
% This function calls the following functions
%   # brEphRinexRead(sFileName)
%   # brEphCoordCalcOneDataSet(vObs,cObsTime)
%
%
% The function is written by Johan Vium Andersson
% Last time updated 2005-01-31

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START

%DIFLIM = 3*3600; %maximum allowed time difference between toe and vObsTime

% read broadcast ephemeris file and return the header and a data matrix
% with the observations
    
    [n m]=size(mData);
    
% Get all data from one and the same satellite (:,m)

XYZ=[];
vCoordinateVector = zeros(32,5);  %allocate space for 32 satellites
dtStart = get(mData(1,1), 'dtStart');
dtEnd = get(mData(1,1), 'dtEnd');
if (vObsTime - dtStart) < 0 || (vObsTime - dtEnd) > 0
    
    error('No broadcast ephemeris for this epoch');
end
NoData = 1;
for i=1:n  %loop over satellite

    A=[];
    MinDiff = 1e15;  %Time difference in [s]
    MinIdx = 0;  %index to the closest epoch
    for j=1:m  %loop over toe 
        %B=mData(i,j)
        
        B(1) = get(mData(i,j) , 'iEpochYear');
        B(2) = get(mData(i,j) , 'iEpochMonth'); 
        B(3) = get(mData(i,j) , 'iEpochDay');               
        B(4) = get(mData(i,j) , 'iEpochHour');               
        B(5) = get(mData(i,j) , 'iEpochMinute');        
        B(6) = get(mData(i,j) , 'iEpochSecond');    
        if abs(B(1)) < 1e-13  %No data
            %NoData = 1;
            continue;
        end
        
        dJDSat = DateTime(B); % Convert the time  
        Tdif = abs(dJDSat-vObsTime);
        if Tdif < MinDiff
            MinDiff = Tdif;
            MinIdx = j;
            NoData = 0;
        end
    end
%     if MinDiff > DIFLIM  %no eph data for current epoch
%         continue;
%     end
    if NoData
        continue;
    end
    [x,y,z,ds,TGD] =  brEphCoordCalcOneDataSet( mData(i, MinIdx), vObsTime );
    vCoordinateVector(i,:) = [x y z ds TGD];
    NoData = 1;
end

