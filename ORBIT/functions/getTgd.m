function [TGD] = getTgd(sFileName)
% This function collect all the time group delays from a broadcast (TGD)
% navigation file
%
% Input:
% sFileName = file name broadcast navigation file
%
% TGD

sz = zeros(32,10); %allocate space for 32 satellites
TGD = struct('vDateTime',DateTime,'cTGD',sz);
fid=fopen(sFileName);


%% ----------------------------------------------------------
rad=0; % Start value Counter
num=0; % Start value Counter
tline= fgetl(fid);
z=[]; % z is a vector with all the satellites in the file (see below)

while num<1 % Loop until END OF HEADER
    tline= fgetl(fid);
    matches = findstr(tline,'END OF HEADER');
    num = length(matches);
end

%% ----------------------------------------------------------
k=0;
while feof(fid) == 0; % Loop continues until the end of the file
    % PRN / EPOCH / SV CLK
        tline= fgetl(fid);
        rad=rad+1; % Take a new line from the file
        
        iPRN =str2num(tline(1:2));
         year=str2double(tline(3:5));
            % Convert the year to 4 digets
            if year<94;
                year=year+2000;
            else
                year=year+1900;
            end       

         ye=year;
         mo=str2num(tline(6:8));    
         da=str2num(tline(9:11));
         h=str2num(tline(12:14));
         m=str2num(tline(15:18));        
         s=str2num(tline(19:22));       
        
        if  k == 0 %abs(TGD.vDateTime(1) - DateTime(ye,mo,da,0,0,0)) >= 86400 % New day 
            k = 1;
            TGD.vDateTime(k) = DateTime(ye,mo,da,0,0,0);
        end
                 
        for i = 1:6
            tline= fgetl(fid);
            rad=rad+1;
        end          
    
        TGD.cTGD(iPRN,k)=str2num(tline(42:60));
        
        tline= fgetl(fid);

end %while

fid1 = fopen('TGD.txt','w');

fprintf(fid1,'TGD data from file :%s\n\nThe data covers %2d days\n\n',sFileName,k);

for j=1:k
    fprintf(fid1,'Year   month    day\n');
    fprintf(fid1,'%4d   %2d       %2d \n\n',ye,mo,da);
    fprintf(fid1,'iPRN    TGD \n\n');
    for i = 1:length(TGD.cTGD)
        if TGD.cTGD(i) ~= 0
            fprintf(fid1,'%2d     %2.12e\n',i,TGD.cTGD(i,j) );
        end
    end
end
%fclose('all');


% Return the the header and all data