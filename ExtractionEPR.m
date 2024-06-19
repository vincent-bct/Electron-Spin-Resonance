function [File_list,B,spc,mw_frequency,Comments,g_factor,Day,Time,nb_scans,mw_power,receiver_gain] = ExtractionEPR(varargin)

addpath(genpath('C:\Users\vince\Documents\MATLAB\easyspin-6.0.0-dev.53'));
[File_list,FilePath] = uigetfile({'*.spc;*.xlsx; *.xlsm',' Compatible Files (*.spc,*.xlsx,*.xlsm)'}, ...
                      'Select one or More Files', 'MultiSelect', 'on');
% Let's handle cases where:
% only one file is selected: Uigetfile then returns a char chain -->  conversion into a cell.
if ischar(File_list) == 1
    File_list = cellstr(File_list);
end
%Let's handle the case where no file is selected. Then uigetfile returns 0.
%As the case where multiple or only one file have been selected now both returns a cell array. 
% We can use this to distinguish the case where = 0
if ~iscell(File_list)
    return
end
nb_samples = length(File_list); % Nb of data files to be treated
%Initialization of the storages matrices
% B and spc matrices need to be initialized in the loop as their length will
% depend on the nb of points returned in the measure (ANZ in Bruker files)
B = [];
spc =[];
g_factor = [];
nb_scans = zeros(1,nb_samples);
mw_frequency = zeros(1,nb_samples);
mw_power = zeros(1,nb_samples);
receiver_gain = zeros(1,nb_samples);
Comments = cell(1,nb_samples);
Day = cell(1,nb_samples); 
Time = cell(1,nb_samples);
%Constants
Bohr_magneton = 9.274009994*1E-24; %Joules.Tesla−1
h_Planck = 6.62607015E-34; % Joules.sec
n_sup = 0;%Necessary to account for excel files with multiple curves
% q_el  = 1.602176634e-19; %Coulomb
% k_B  =1.380649E-23 ; % Joules.Kelvin-1

% Extraction loop: the data from each file is extracted
for i =1:1:nb_samples 
   FileName = strcat(char(FilePath), char(File_list(i)));%We get the name
   % We need to distinguish 2 cases: re-import from an .xlsx file or import with .spc files
   if endsWith(FileName,".spc")
       % The easyspin function eprload allows direct extraction from Bruker files
        % 'n' denotes a scaling by the nb of scans performed. (/n)
        % T denotes a scaling by the temperature in Kelvin, considering a Curie law (*T)
        % P denotes divides by sqrt(p_µwave) with p_µwave in mW
        clear B_i spc_i
        if nargin == 1
            [B_i(:,1),spc_i(:,1),param_i] = eprload(FileName,varargin{1});    
        elseif nargin ==0
            [B_i(:,1),spc_i(:,1),param_i] = eprload(FileName);
        end
%        nb_Xpts = param_i.ANZ; % We collect the resolution (nb of x points)
       % Let's take into account the case where files with multiple
       % resolutions are loaded at the same time:the B storage or B_i need to be extended
       if length(B_i) > length(B)
           diff =  length(B_i) - length(B);
           B = [B;nan(diff,size(B,2))];
           spc = [spc;nan(diff,size(B,2))];
       end
       if length(B_i) <length(B)
           diff = length(B) - length(B_i);
           B_i = [B_i;nan(diff,size(B_i,2))];
           spc_i = [spc_i;nan(diff,size(B_i,2))];
       end
       B = [B, B_i];
       spc = [spc, spc_i];
       nb_max_res = length(B);
       nb_scans(i)= param_i.JSD;
       mw_frequency(i)= param_i.MF*1e9;%En Hz
       mw_power(i) = param_i.MP;
       receiver_gain(i) = param_i.RRG;
       Comments(i) = {param_i.JCO};
       Day(i) = {param_i.JDA};
       Time(i) = {param_i.JTM};
       
       for j=1:1:nb_max_res % Calculating g, as B is expressed in Gauss, conversion to SI unit tesla (x1e-4)
           g_factor(j,i) = (h_Planck*mw_frequency(i))/(1e-4*B(j,i)*Bohr_magneton);
       end  
   elseif endsWith(FileName,".xlsx")||endsWith(FileName,".xlsm")
       % We need to take into account the case where multiple curves are
       % present in the excel files, so the storages have to be extended accordingly
       clear B_i spc_i Curves_data_i Parameters_i % Clear previous variables
       % Loading from excel file...
       Curves_data_i = readtable( FileName,'Sheet','Curves');
       Curves_data_i = table2array(Curves_data_i);%We get a matrix alternating B|spc|B_2|spc_2|...
       Parameters_i = readtable(FileName,'Sheet','Parameters');
       Parameters_i = table2cell(Parameters_i); % We get a cell array
       Nb_Curves = size(Curves_data_i,2)/4; % Normally, the input should be [B, B_corr, g, spc]
       if Nb_Curves > 1 %We need to extend storage for parameters too
           nb_scans = [nb_scans, zeros(1,Nb_Curves-1)];
           mw_frequency = [mw_frequency, zeros(1,Nb_Curves-1)];
           mw_power = [mw_power, zeros(1,Nb_Curves-1)];
           receiver_gain = [receiver_gain, zeros(1,Nb_Curves-1)];
           Comments = [Comments, cell(1,Nb_Curves-1)];
           Day = [Day, cell(1,Nb_Curves-1)];
           Time = [Time, cell(1,Nb_Curves-1)];
           File_list = [File_list(1:i-1),cell(1,Nb_Curves), File_list(i+1:end)];
       end
       for j = 1:Nb_Curves
            %Collecting the curves
            B =  [B, Curves_data_i(:,4*j-2)];% Already converted mili-Tesla to Gauss in save_fit (in case excel is from fit)
            spc = [spc, Curves_data_i(:,4*j)];
            g_factor = [g_factor,Curves_data_i(:,4*j-1)];
            nb_max_res = length(B);
            name = Parameters_i{j};
            mw_frequency(i+j-1+n_sup) = Parameters_i{j,3} ; %GHz
            nb_scans(i+j-1+n_sup)= Parameters_i{j,2};
            mw_power(i+j-1+n_sup) = Parameters_i{j,4};
            receiver_gain(i+j-1+n_sup) = Parameters_i{j,5};
            Day{i+j-1+n_sup} = Parameters_i{j,6};
            Time{i+j-1+n_sup} = Parameters_i{j,7};
            Comments{i+j-1+n_sup} = Parameters_i{j,8}; 
            File_list{1,i+j-1} = name;
       end
       n_sup = n_sup + Nb_Curves-1;
           
   end
end

end



