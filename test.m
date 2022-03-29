clc
clear all
addpath('ORBIT/functions/');
addpath('ORBIT/Classes');
 %Obj=LSSolver();
settings.ws = 'ws/sacc_2008_245_254/';
settings.beginDateTime = DateTime(2008,9,5,1,0,0);
settings.endDateTime   = DateTime(2008,9,5,1,0,10);
%get(settings.beginDateTime,'MJD')
% Observation type
settings.obsType='Graphic';
% Open GPS broadcast ephemerides file
    orbBrdEphFileObj=OrbBrdEphFile(settings);
    % Open GPS observation file  
    orbObsFileObj=OrbObsFile(settings,'rinex');                   
    outFileObj=OrbFilterOutputFile(settings);
    
  % Run the filter
    baseObj=LSSolver();
   Obj=baseObj.getInstance(settings,orbObsFileObj, ...
                                       outFileObj,orbBrdEphFileObj)
   Obj.run()