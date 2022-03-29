classdef EphHeader
    properties (SetAccess=public)
        iRinexVersion,sFileType,sPGM,sRunBy,sDate,sComment,dAlfaIon1,dAlfaIon2,dAlfaIon3
        dAlfaIon4,dBetaIon1,dBetaIon2,dBetaIon3,dBetaIon4,dA0,dA1,dRefTime,dRefWeek,iLeapSeconds
                
    end

    methods
       
        function this = EphHeader(indata)
        % Navigation file header data, stored in the following parameters
        %
        % RINEX VERSION         : iRinexVersion
        % File Type             : sFileType
        % Name of Program       : sPGM
        % Name of Agency        : sRunBy
        % Date of File Creation : sDate
        % Comment               : sComment
        % Ionosphere param. A0  : dAlfaIon1
        % Ionosphere param. A1  : dAlfaIon2
        % Ionosphere param. A2  : dAlfaIon3
        % Ionosphere param. A3  : dAlfaIon4
        % Ionosphere param. B0  : dBetaIon1
        % Ionosphere param. B1  : dBetaIon2
        % Ionosphere param. B2  : dBetaIon3
        % Ionosphere param. B3  : dBetaIon4
        % Delta-UTC A0          : dA0
        % Delta-UTC A1          : dA1
        % Reference Time        : dRefTime
        % Reference Week        : dRefWeek
        % Leap second           : iLeapSeconds
        %
        %
        % The class contains the following functions
        %
        %   display(A)    - All data in A are displayed in fprintf style
        %   set(A,par,in) - Set one of the parameters in A
        %   get(A,par)    - Get one parameter
        %
        
        switch nargin
            
            case 0
            % no arguments
                this.iRinexVersion=0;this.sFileType='';this.sPGM='';this.sRunBy='';
                this.sDate='';this.sComment='';
                this.dAlfaIon1=0;this.dAlfaIon2=0;this.dAlfaIon3=0;this.dAlfaIon4=0;
                this.dBetaIon1=0;this.dBetaIon2=0;this.dBetaIon3=0;this.dBetaIon4=0
                this.dA0=0;this.dA1=0;this.dRefTime=0;this.dRefWeek=0;
                this.iLeapSeconds=0;
                
            case 1 
        
                 this.iRinexVersion=indata.iRinexVersion;
                 this.sFileType=indata.sFileType;this.sPGM=indata.sPGM;this.sRunBy=indata.sRunBy;
                this.sDate=indata.sDate;this.sComment=indata.sComment;
                this.dAlfaIon1=indata.dAlfaIon1;this.dAlfaIon2=indata.dAlfaIon2;
                this.dAlfaIon3=indata.dAlfaIon3;this.dAlfaIon4=indata.dAlfaIon4;
                this.dBetaIon1=indata.dBetaIon1;this.dBetaIon2=indata.dBetaIon2;
                this.dBetaIon3=indata.dBetaIon3;this.dBetaIon4=indata.dBetaIon4;
                this.dA0=indata.dA0;this.dA1=indata.dA1;this.dRefTime=indata.dRefTime;
                this.dRefWeek=indata.dRefWeek;
                this.iLeapSeconds=indata.iLeapSeconds;
              
        end
        end
        function a = get(b,par)

            a=b.(par);
        end
        function IonoPar = getIonoPar(this)
        % This function returns the broadcast Ionospheric parameters from a 
        % navigation file header
        
        IonoPar = [this.dAlfaIon1,this.dAlfaIon2,this.dAlfaIon3,this.dAlfaIon4;
                  this.dBetaIon1,this.dBetaIon2,this.dBetaIon3,this.dBetaIon4];
        end
        function b = set(b,par,in)
            b.(par)=in;
        end
        function display(b)
        
        fprintf('RINEX VERSION         : %1.0i \n',b.iRinexVersion);
        fprintf('File Type             : %s\n'   ,b.sFileType);
        fprintf('Name of Program       : %s\n'   ,b.sPGM);
        fprintf('Name of Agency        : %s\n'   ,b.sRunBy);
        fprintf('Date of File Creation : %s\n'   ,b.sDate);
        fprintf('Comment               : %s\n'   ,b.sComment);
        fprintf('Ionosphere param. A0  : %2.0i\n',b.dAlfaIon1);
        fprintf('Ionosphere param. A1  : %2.0i\n',b.dAlfaIon2);
        fprintf('Ionosphere param. A2  : %2.0i\n',b.dAlfaIon3);
        fprintf('Ionosphere param. A3  : %2.0i\n',b.dAlfaIon4);
        fprintf('Ionosphere param. B0  : %2.0i\n',b.dBetaIon1);
        fprintf('Ionosphere param. B1  : %2.0i\n',b.dBetaIon2);
        fprintf('Ionosphere param. B2  : %2.0i\n',b.dBetaIon3);
        fprintf('Ionosphere param. B3  : %2.0i\n',b.dBetaIon4);
        fprintf('Delta-UTC A0          : %2.0i\n',b.dA0);
        fprintf('Delta-UTC A1          : %2.0i\n',b.dA1);
        fprintf('Reference Time        : %2.0i\n',b.dRefTime);
        fprintf('Reference Week        : %2.0i\n',b.dRefWeek);
        fprintf('Leap second           : %2.0i\n\n',b.iLeapSeconds);
        end


        
                
    end
end


   
    
    
