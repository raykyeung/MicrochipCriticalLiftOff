%% B2KCriticalLiftOff.m - [Script] Critical Lift-Off Analysis
%
% Author: Raymond Yeung
% Release date: 2025
% E-mail: ryeun003@ucr.edu
% B2K Group, Dept. of Bioengineering, Univ. of California, Riverside
% Victor G. J. Rodgers Dept. of Bioengineering, Univ. of California, Riverside
% William H. Grover, Dept. of Bioengineering, Univ. of California, Riverside
% Philip L. Brisk Dept. of Computer Science, Univ. of California, Riverside

%% Close all figures, clear all workspace variables, clear command window
close all;
clear all;
clc;

%% Elapsed time, tic (start)   
tic 

%% Flags
flagCombine = 1; %Combine or separate data points

%% Base path for exporting figures and files
%---[Get the script file name]---
[~, scriptName, ~] = fileparts(mfilename('fullpath'));

%---[Define the base export path location]---
baseExportPath = fullfile(getenv('USERPROFILE'),'Downloads');

%---[Define subfolder]---
exportSubFolder = scriptName;

%---[Construct full base path]---
fullBaseExportPath = fullfile(baseExportPath,exportSubFolder);

%---[Check if the subfolder exists]---
if ~exist(fullBaseExportPath, 'dir')
    %---[Create the subfolder if it doesn't exist]---
    mkdir(fullBaseExportPath);
else
    %---[If the subfolder exists, delete all files in it]---
    files = dir(fullfile(fullBaseExportPath, '*'));
    for k = 1:length(files)
        %---[Skip '.' current directory and '..' parent directory fields]---
        if ~files(k).isdir || (~strcmp(files(k).name, '.') && ~strcmp(files(k).name, '..'))
            delete(fullfile(fullBaseExportPath, files(k).name));
        end
    end
end

%% Project root path
proj = matlab.project.currentProject;
root = proj.RootFolder;

%% Setup for log file to list dependencies

% Current script name ('Batch1pTracking01')
currScriptName = mfilename; %current script name

% Date ('230613')
currDate = char(datetime('now','Format','yyMMdd'));

% Combine current script name and date ('Batch1pTracking01_230613')
currScriptNameAddDate = append(currScriptName,'_',currDate); %current script name with current date

%% nParameter definition using a structure (paramStruct)
%---[List of parameter names as a cell array of strings]---
paramNames = { ...
    'xParameter', ...               %[1]
    'yParameter', ...               %[2]
    'zParameterHW_D_maxSq', ...     %[3]
    'zParameterHW_D_VSq', ...       %[4]
    'aParameter', ...               %[5]
    'mParameter', ...               %[6]
    'rParameter', ...               %[7]
    'hD_maxParameter', ...          %[8]
    'wD_maxParameter', ...          %[9]
    'hD_VParameter',...             %[10]
    'wD_VParameter'};               %[11]

%---[Define the default parameter structure with all fields set to empty arrays]---
defaultParam = struct( ...
    'Str',               [], ...
    'Unit',              [], ...
    'LegendNumStr',      [], ...
    'LegendDenStr',      [], ...
    'LegendStr',         [], ...
    'Value',             [], ...
    'Uncertainty',       [], ...
    'ValueSorted',       [], ...
    'UncertaintySorted', [], ...
    'NewStr',            []);

%---[Create the paramStruct and assign defaultParam to each of the parameter names]---
paramStruct = struct();
for i = 1:length(paramNames)
    paramStruct.(paramNames{i}) = defaultParam;
end

%% Chosen xParameter %[1]
paramStruct.xParameter.Str = ["ArchimedesNumberUsingD_max"];
paramStruct.xParameter.Unit = [1];

%% Chosen yParameter %[2]
paramStruct.yParameter.Str = ["ShearReynoldsNumberUsingD_maxCRIT"];
paramStruct.yParameter.Unit = [1];

%% Chosen zParameterHW_D_maxSq %[3]
paramStruct.zParameterHW_D_maxSq.Str = ["RatioChannelWidthD_max" "RatioChannelHeightD_max"];
paramStruct.zParameterHW_D_maxSq.Unit = [1 1];

%% Chosen zParameterHW_D_VSq %[4]
paramStruct.zParameterHW_D_VSq.Str = ["RatioChannelWidthD_V" "RatioChannelHeightD_V"];
paramStruct.zParameterHW_D_VSq.Unit = [1 1];

%% Chosen aParameter %[5]
paramStruct.aParameter.Str = ["ChannelAspectRatio"];
paramStruct.aParameter.Unit = [1];

%% Chosen mParameter %[6]
paramStruct.mParameter.Str = ["ModifiedArchimedesNumberUsingD_V"];
paramStruct.mParameter.Unit = [1];

%% Chosen rParameter %[7]
paramStruct.rParameter.Str = ["ModifiedParticleReynoldsNumberUsingU_avgD_VCRIT"];
paramStruct.rParameter.Unit = [1];

%% Chosen hD_maxParameter %[8]
paramStruct.hD_maxParameter.Str = ["RatioChannelHeightD_max"];
paramStruct.hD_maxParameter.Unit = [1];

%% Chosen wD_maxParameter %[9]
paramStruct.wD_maxParameter.Str = ["RatioChannelWidthD_max"];
paramStruct.wD_maxParameter.Unit = [1];

%% Chosen hD_VParameter %[10]
paramStruct.hD_VParameter.Str = ["RatioChannelHeightD_V"];
paramStruct.hD_VParameter.Unit = [1];

%% Chosen wD_VParameter %[11]
paramStruct.wD_VParameter.Str = ["RatioChannelWidthD_V"];
paramStruct.wD_VParameter.Unit = [1];

%% Input all experimental data
%---[dataCase options: 'All', 'H>=W', 'H>W', 'H<=W', 'H<W', 'H=W']---
dataCase = 'All';
% dataCase = 'H>W';
% dataCase = 'H=W';
% dataCase = 'H<W';
switch dataCase
    case 'All'
        %--Arguments for AllParameters-- %1,2,3,4,5,6,7,8,9,10,11 [11,3]
        ChannelHeight = repmat([1.000;1.500;1.500;0.889;1.333;1.131;2.263;1.777;2.000;1.500;1.697],3,1);
        ChannelWidth = repmat([1.500;1.000;1.500;1.333;1.333;1.697;1.697;1.333;1.500;1.500;1.697],3,1);
        Solvent = [repmat({"Water"},11,1);repmat({"Isopropanol"},11,1);repmat({"Methanol"},11,1)];
        CritFlowRate = ...
            [9;17;21;8;14;13;36;19;25;17;25;...
            5;16;9;7;9;10;22;14;17;12;17;...
            9;14;17;7;12;12;35;20;31;16;25];
        sizeCase = [11,3];
    case 'H>=W'
        %--Arguments for AllParameters, H>=W-- %2,3,5,7,8,9,10,11 [8,3]
        %H>=W: 2,7,8,9
        % equals: 3,5,10,11
        ChannelHeight = repmat([1.500;1.500;1.333;2.263;1.777;2.000;1.500;1.697],3,1);
        ChannelWidth = repmat([1.000;1.500;1.333;1.697;1.333;1.500;1.500;1.697],3,1);
        Solvent = [repmat({"Water"},8,1);repmat({"Isopropanol"},8,1);repmat({"Methanol"},8,1)];
        CritFlowRate = ...
            [17;21;14;36;19;25;17;25;...
            16;9;9;22;14;17;12;17;...
            14;17;12;35;20;31;16;25];
        sizeCase = [8,3];
    case 'H>W'
        %--Arguments for AllParameters, H>=W-- %2,7,8,9 [4,3]
        %H>=W: 2,7,8,9
        ChannelHeight = repmat([1.500;2.263;1.777;2.000],3,1);
        ChannelWidth = repmat([1.000;1.697;1.333;1.500],3,1);
        Solvent = [repmat({"Water"},4,1);repmat({"Isopropanol"},4,1);repmat({"Methanol"},4,1)];
        CritFlowRate = ...
            [17;36;19;25;...
            16;22;14;17;...
            14;35;20;31];
        sizeCase = [4,3];
    case 'H<=W'
        %--Arguments for AllParameters, H<=W-- %1,3,4,5,6,10,11 [7,3]
        %H<=W: 1,4,6
        % equals: 3,5,10,11
        ChannelHeight = repmat([1.000;1.500;0.889;1.333;1.131;1.500;1.697],3,1);
        ChannelWidth = repmat([1.500;1.500;1.333;1.333;1.697;1.500;1.697],3,1);
        Solvent = [repmat({"Water"},7,1);repmat({"Isopropanol"},7,1);repmat({"Methanol"},7,1)];
        CritFlowRate = ...
            [9;21;8;14;13;17;25;...
            5;9;7;9;10;12;17;...
            9;17;7;12;12;16;25];
        sizeCase = [7,3];
    case 'H<W'
        %--Arguments for AllParameters, H<=W-- %1,4,6 [3,3]
        %H<=W: 1,4,6
        ChannelHeight = repmat([1.000;0.889;1.131],3,1);
        ChannelWidth = repmat([1.500;1.333;1.697],3,1);
        Solvent = [repmat({"Water"},3,1);repmat({"Isopropanol"},3,1);repmat({"Methanol"},3,1)];
        CritFlowRate = ...
            [9;8;13;...
            5;7;10;...
            9;7;12];
        sizeCase = [3,3];
    otherwise %'H=W'
        %--Arguments for AllParameters-- %3,5,10,11 [4,3]
        % equals: 3,5,10,11
        ChannelHeight = repmat([1.500;1.333;1.500;1.697],3,1);
        ChannelWidth = repmat([1.500;1.333;1.500;1.697],3,1);
        Solvent = [repmat({"Water"},4,1);repmat({"Isopropanol"},4,1);repmat({"Methanol"},4,1)];
        CritFlowRate = ...
            [21;14;17;25;...
            9;9;12;17;...
            17;12;16;25];
        sizeCase = [4,3];
end

%% Input all Patankar 2001 Data [x- Archimedes number \mathrm{Ar}_{D_\mathrm{max}}; y- Critical shear Reynolds number \mathrm{\check{Re}}_{\mathrm{s},D_\mathrm{max}} ]
%
pH_d48x = ...
[7.797568557
57.77446247
73.42404176
280.3895538
568.9501723
697.8831536
881.3575516
944.7237183
6362.997496];

pH_d48y = ...
[2.362916158
9.797845174
11.69956431
29.43761168
49.09698615
58.44756139
71.45364724
78.76068066
295.405979];

%
pH_d12x = ...
[1.784632472
5.211159696
10.52918523
24.05079338
189.3000071];

pH_d12y = ...
[0.995567877
1.97767971
2.975476548
4.935326824
23.64574459];

%
pH_d6x = ...
[40.75122784
188.6403511
369.4303919
647.0014619];

pH_d6y = ...
[9.81579469
29.40041915
49.02940407
78.43440554];

%
pH_d4x = ...
[26.77547675
459.0702087
263.6789934
133.0575067];

pH_d4y = ...
[9.802653163
78.579579
49.12107897
29.45428031];

%% Read all Excel file names [Data 01-11]
DFiles_01_11_sub1 = 'ExperimentalScripts';
DFiles_01_11_sub2 = 'COMSOLDataProcessing';
DFiles_01_11_sub3 = 'COMSOLDataExportExcel';
DFiles_01_11_sub4 = 'Data_01_11';

DFiles_01_11 = dir(fullfile(root,DFiles_01_11_sub1,DFiles_01_11_sub2,DFiles_01_11_sub3,DFiles_01_11_sub4));
DFiles_01_11 = DFiles_01_11(~ismember({DFiles_01_11.name},{'.','..'}));

%% Read all Excel files and store data in cell [Data 01-11]
C_data_01_11 = cell(1,length(DFiles_01_11));

for i = 1:length(DFiles_01_11)
    C_data_01_11{1,i} = readmatrix(DFiles_01_11(i).name);
end

%% C_Data_01_11 organization and order of exported data
%channel[1](10)*
%-surface[1](6)* *Stored in (60) separate Excel files
%--totalDrag(6)**
%--pressureDrag**
%--viscousDrag**
%--totalLift**
%--pressureLift**
%--viscousLift**
%---solvent(3)
%----xDistInletObj(5)
%-----botClear(3)

%6*3*5*3 = 270

%% Store C_Data_01_11 data in nested structure
%channel(10)
%-surface(6)
%--botClear(3)
%---xDistInletObj(5)
%----solvent(3)
%-----totalDrag
%-----pressureDrag
%-----viscousDrag
%-----totalLift
%-----pressureLift
%-----viscousLift

%% Parsing data_01-11 
channelCounter_01_11 = 1; %10 [done]
surfaceCounter_01_11 = 0; %6 [done]
botClearCounter_01_11 = 1; %3 [done] % 0.01, 0.05, 0.10
xDistInletObjCounter_01_11 = 1; %5 [done]
solventCounter_01_11 = 1; %3 [done]
evaluateCounter_01_11 = 0; %6 [done]

for i_Cell_data_01_11 = 1:length(C_data_01_11)

    surfaceCounter_01_11 = surfaceCounter_01_11 + 1;

    for i_Table_01_11 = 1:length(C_data_01_11{1,1})
        if rem(i_Table_01_11,6) == 0
            evaluateCounter_01_11 = 6;
        else
            evaluateCounter_01_11 = rem(i_Table_01_11,6);
        end

        switch evaluateCounter_01_11
            case 1
                channel_01_11(channelCounter_01_11).surface_01_11(surfaceCounter_01_11).botClear_01_11(botClearCounter_01_11).xDistInletObj_01_11(xDistInletObjCounter_01_11).solvent_01_11(solventCounter_01_11).totalDrag = C_data_01_11{1,i_Cell_data_01_11}(i_Table_01_11)*2;
            case 2
                channel_01_11(channelCounter_01_11).surface_01_11(surfaceCounter_01_11).botClear_01_11(botClearCounter_01_11).xDistInletObj_01_11(xDistInletObjCounter_01_11).solvent_01_11(solventCounter_01_11).pressureDrag = C_data_01_11{1,i_Cell_data_01_11}(i_Table_01_11)*2;
            case 3
                channel_01_11(channelCounter_01_11).surface_01_11(surfaceCounter_01_11).botClear_01_11(botClearCounter_01_11).xDistInletObj_01_11(xDistInletObjCounter_01_11).solvent_01_11(solventCounter_01_11).viscousDrag = C_data_01_11{1,i_Cell_data_01_11}(i_Table_01_11)*2;
            case 4
                channel_01_11(channelCounter_01_11).surface_01_11(surfaceCounter_01_11).botClear_01_11(botClearCounter_01_11).xDistInletObj_01_11(xDistInletObjCounter_01_11).solvent_01_11(solventCounter_01_11).totalLift = C_data_01_11{1,i_Cell_data_01_11}(i_Table_01_11)*2;
            case 5
                channel_01_11(channelCounter_01_11).surface_01_11(surfaceCounter_01_11).botClear_01_11(botClearCounter_01_11).xDistInletObj_01_11(xDistInletObjCounter_01_11).solvent_01_11(solventCounter_01_11).pressureLift = C_data_01_11{1,i_Cell_data_01_11}(i_Table_01_11)*2;
            case 6
                channel_01_11(channelCounter_01_11).surface_01_11(surfaceCounter_01_11).botClear_01_11(botClearCounter_01_11).xDistInletObj_01_11(xDistInletObjCounter_01_11).solvent_01_11(solventCounter_01_11).viscousLift = C_data_01_11{1,i_Cell_data_01_11}(i_Table_01_11)*2;
        end

        if rem(i_Table_01_11,6) == 0 %6
            solventCounter_01_11 = solventCounter_01_11 + 1;
        end
        if rem(i_Table_01_11,18) == 0 %6*3
            solventCounter_01_11 = 1;
        end

        if rem(i_Table_01_11,18) == 0 %6*3
            xDistInletObjCounter_01_11 = xDistInletObjCounter_01_11 + 1;
        end
        if rem(i_Table_01_11,90) == 0 %6*3*5
            xDistInletObjCounter_01_11 = 1;
        end

        if rem(i_Table_01_11,90) == 0 %6*3*5
            botClearCounter_01_11 = botClearCounter_01_11 + 1;
        end
    end

    botClearCounter_01_11 = 1; %3
    xDistInletObjCounter_01_11 = 1; %5
    solventCounter_01_11 = 1; %3
    evaluateCounter_01_11 = 0; %6

    if rem(i_Cell_data_01_11,6) == 0
        channelCounter_01_11 = channelCounter_01_11 + 1;
    end

    if rem(i_Cell_data_01_11,6) == 0
        surfaceCounter_01_11 = 0;
    end
end

%% channel_01_11 values
%channel(10)
%-surface(6)
%--botClear(3)
%---xDistInletObj(5)
%----solvent(3)

channelHeight_Val_01_11 = [1.000;1.500;1.500;0.889;1.333;1.131;2.263;1.777;2.000;1.697];
channelWidth_Val_01_11 = [1.500;1.000;1.500;1.333;1.333;1.697;1.697;1.333;1.500;1.697];
channel_Count_01_11 = [1 2 3 4 5 6 7 8 9 10];
surface_Str_01_11 = ["all" "top" "bottom" "front" "left" "back"];
botClear_Val_01_11 = [1E-5 5E-5 1E-4];
xDistInletObj_Val_01_11 = [1E-4 9E-4 0.0024 0.00301 0.023];
solvent_Str_01_11 = ["Isopropanol" "Water" "Methanol"];

%% Read all Excel file names [Data 00]
DFiles_00_sub1 = 'ExperimentalScripts';
DFiles_00_sub2 = 'COMSOLDataProcessing';
DFiles_00_sub3 = 'COMSOLDataExportExcel';
DFiles_00_sub4 = 'Data_00';

DFiles_00 = dir(fullfile(root,DFiles_00_sub1,DFiles_00_sub2,DFiles_00_sub3,DFiles_00_sub4));

DFiles_00 = DFiles_00(~ismember({DFiles_00.name},{'.','..'}));

%% Read all Excel files and store data in cell [Data 00]
C_data_00 = cell(1,length(DFiles_00));

for i = 1:length(DFiles_00)
    C_data_00{1,i} = readmatrix(DFiles_00(i).name);
end

%% C_Data_00 organization and order of exported data
%surface[1](6) *Stored in (6) separate Excel files
%-totalDrag(6)**
%-pressureDrag**
%-viscousDrag**
%-totalLift**
%-pressureLift**
%-viscousLift**
%----mesh(11)
%---botClear(5)
%--solvent*channel(3*10)

%30*5*11*6 = 9900

%% Store C_Data_00 data in nested structure

%channel(10)
%-surface(6)
%--botClear(5)
%---mesh(11)
%----solvent(3)
%-----totalDrag
%-----pressureDrag
%-----viscousDrag
%-----totalLift
%-----pressureLift
%-----viscousLift

%% Parsing data_00

surfaceCounter_00 = 0; %6 [done]
channelCounter_00 = 1; %10*{} {Channel 01-11; Isopropanol, Water, Methanol} [done]
solventCounter_00 = 1; %{}*3 {Channel 01-11; Isopropanol, Water, Methanol} [done]
botClearCounter_00 = 1; %5 [done]
meshCounter_00 = 1; %11 [done]
evaluateCounter_00 = 0; %6 [done]

for i_Cell_data_00 = 1:length(C_data_00) %6

    surfaceCounter_00 = surfaceCounter_00 + 1;

    for i_Table_00 = 1:length(C_data_00{1,1}) %9900
        if rem(i_Table_00,6) == 0
            evaluateCounter_00 = 6;
        else
            evaluateCounter_00 = rem(i_Table_00,6);
        end

        switch evaluateCounter_00
            case 1
                channel_00(channelCounter_00).surface_00(surfaceCounter_00).botClear_00(botClearCounter_00).mesh_00(meshCounter_00).solvent_00(solventCounter_00).totalDrag = C_data_00{1,i_Cell_data_00}(i_Table_00)*2;
            case 2
                channel_00(channelCounter_00).surface_00(surfaceCounter_00).botClear_00(botClearCounter_00).mesh_00(meshCounter_00).solvent_00(solventCounter_00).pressureDrag = C_data_00{1,i_Cell_data_00}(i_Table_00)*2;
            case 3
                channel_00(channelCounter_00).surface_00(surfaceCounter_00).botClear_00(botClearCounter_00).mesh_00(meshCounter_00).solvent_00(solventCounter_00).viscousDrag = C_data_00{1,i_Cell_data_00}(i_Table_00)*2;
            case 4
                channel_00(channelCounter_00).surface_00(surfaceCounter_00).botClear_00(botClearCounter_00).mesh_00(meshCounter_00).solvent_00(solventCounter_00).totalLift = C_data_00{1,i_Cell_data_00}(i_Table_00)*2;
            case 5
                channel_00(channelCounter_00).surface_00(surfaceCounter_00).botClear_00(botClearCounter_00).mesh_00(meshCounter_00).solvent_00(solventCounter_00).pressureLift = C_data_00{1,i_Cell_data_00}(i_Table_00)*2;
            case 6
                channel_00(channelCounter_00).surface_00(surfaceCounter_00).botClear_00(botClearCounter_00).mesh_00(meshCounter_00).solvent_00(solventCounter_00).viscousLift = C_data_00{1,i_Cell_data_00}(i_Table_00)*2;
        end

        if rem(i_Table_00,6) == 0 %6
            meshCounter_00 = meshCounter_00 + 1;
        end
        if rem(i_Table_00,66) == 0 %6*11
            meshCounter_00 = 1;
        end

        if rem(i_Table_00,66) == 0 %6*11
            botClearCounter_00 = botClearCounter_00 + 1;
        end
        if rem(i_Table_00,330) == 0 %6*11*5
            botClearCounter_00 = 1;
        end

        if rem(i_Table_00,330) == 0 %6*11*5
            solventCounter_00 = solventCounter_00 + 1;
        end
        if rem(i_Table_00,990) == 0 %6*11*5*3
            solventCounter_00 = 1;
        end

        if rem(i_Table_00,990) == 0 %6*11*5*3
            channelCounter_00 = channelCounter_00 + 1;
        end
    end

    channelCounter_00 = 1;
    solventCounter_00 = 1;
    botClearCounter_00 = 1;
    meshCounter_00 = 1;
    evaluateCounter = 0;
end

%% channel_00 values
%channel(10)
%-surface(6)
%--botClear(5)
%---mesh(11)
%----solvent(3)

channelHeight_Val_00 = [1.000;1.500;1.500;0.889;1.333;1.131;2.263;1.777;2.000;1.697];
channelWidth_Val_00 = [1.500;1.000;1.500;1.333;1.333;1.697;1.697;1.333;1.500;1.697];
channel_Count_00 = [1 2 3 4 5 6 7 8 9 10];
surface_Str_00 = ["all" "top" "bottom" "front" "left" "back"];
botClear_Val_00 = [1E-6 5E-6 1E-5 5E-5 1E-4];
mesh_Count_00 = [1 2 3 4 5 6 7 8 9 10 11];
solvent_Str_00 = ["Isopropanol" "Water" "Methanol"];

%% Parameter calculations - for loop
Parameters = cell(size(ChannelHeight,1),1);

for i = 1:size(ChannelHeight,1)
    Parameters{i,1} = B2KCalc.B2KAllParametersLiftOff(DimVar(UC(ChannelHeight(i),B2KCalc.B2KAllParametersLiftOff.ChannelHeightUncertainty), 'mm'),DimVar(UC(ChannelWidth(i),B2KCalc.B2KAllParametersLiftOff.ChannelWidthUncertainty), 'mm'),Solvent{i},DimVar(UC(CritFlowRate(i),B2KCalc.B2KAllParametersLiftOff.CritFlowRateUncertainty), 'mL/min'));
end

%% Parameter properties
parametersProperties = properties(Parameters{1,1}); % Determine property names of Parameters cell

%% Define parameters and strings
paramLists.propertyStr = ... %[40 count]
    ["ChannelHeight";
    "ChannelWidth";
    "Solvent";
    "FluidFlowRateCRIT";
    "ChannelLength";
    "RotationalDiameter";
    "EquivVolSpherDiameter"; %
    "HydraulicDiameter";
    "ChannelAspectRatio";
    "RatioD_maxChannelWidth";
    "RatioChannelWidthD_max";
    "RatioD_maxChannelHeight";
    "RatioChannelHeightD_max";
    "RatioD_maxSqandChannelHeightChannelWidth";
    "RatioChannelHeightChannelWidthandD_maxSq";
    "RatioD_VChannelWidth"; %
    "RatioChannelWidthD_V"; %
    "RatioD_VChannelHeight"; %
    "RatioChannelHeightD_V"; %
    "RatioD_VSqandChannelHeightChannelWidth"; %
    "RatioChannelHeightChannelWidthandD_VSq"; %
    "GravityForce";
    "BuoyantForce";
    "ArchimedesNumberUsingD_max";
    "ArchimedesNumberUsingD_V";
    "AverageFluidVelocityCRIT";
    "MaximumFluidVelocityCRIT";
    "ReynoldsNumberUsingU_avgCRIT";
    "ParticleReynoldsNumberUsingU_avgD_maxCRIT";
    "ParticleReynoldsNumberUsingU_avgD_VCRIT";
    "EntranceLengthTruskeyRectangularDuct";
    "EntranceLengthBergmanPipe";
    "TimeToSteadyState";
    "AvgWallShearRateCRIT";
    "ShearReynoldsNumberUsingD_maxCRIT";
    "ShearReynoldsNumberUsingD_VCRIT"
    "ModifiedParticleReynoldsNumberUsingU_avgD_maxCRIT";
    "ModifiedParticleReynoldsNumberUsingU_avgD_VCRIT";
    "ModifiedArchimedesNumberUsingD_V";
    "FroudeNumberUsingU_avgD_maxCRIT"];

paramLists.abbrev = ... %[40 count]
    ["$H$";                                                     % 'ChannelHeight'
    "$W$";                                                      % 'ChannelWidth'
    "$S$";                                                      % 'Solvent'
    "$\check{Q}$";                                              % 'FluidFlowRateCRIT'
    "$L$";                                                      % 'ChannelLength'
    "$d_\mathrm{max}$";                                         % 'RotationalDiameter'
    "$d_\mathrm{V}$";                                           % 'EquivVolSpherDiameter'
    "$D_{H}$";                                                  % 'HydraulicDiameter'
    "$AR$";                                                     % 'ChannelAspectRatio'
    "$d_\mathrm{max}/W$";                                       % 'RatioD_maxChannelWidth'
    "$W/d_\mathrm{max}$";                                       % 'RatioChannelWidthD_max'
    "$d_\mathrm{max}/H$";                                       % 'RatioD_maxChannelHeight'
    "$H/d_\mathrm{max}$";                                       % 'RatioChannelHeightD_max'
    "${d_\mathrm{max}}^2/{HW}$";                                % 'RatioD_maxSqandChannelHeightChannelWidth'
    "$HW/{d_\mathrm{max}}^2$";                                  % 'RatioChannelHeightChannelWidthandD_maxSq'
    "$d_\mathrm{V}/W$";                                         % 'RatioD_VChannelWidth'
    "$W/d_\mathrm{V}$";                                         % 'RatioChannelWidthD_V'
    "$d_\mathrm{V}/H$";                                         % 'RatioD_VChannelHeight'
    "$H/d_\mathrm{V}$";                                         % 'RatioChannelHeightD_V'
    "${d_\mathrm{V}}^2/{HW}$";                                  % 'RatioD_VSqandChannelHeightChannelWidth'
    "$HW/{d_\mathrm{V}}^2$";                                    % 'RatioChannelHeightChannelWidthandD_VSq'
    "$F_\mathrm{g}$";                                           % 'GravityForce'
    "$F_\mathrm{b}$";                                           % 'BuoyantForce'
    "$\mathrm{Ar}_{d_\mathrm{max}}$";                           % 'ArchimedesNumberUsingD_max'
    "$\mathrm{Ar}_{d_\mathrm{V}}$";                             % 'ArchimedesNumberUsingD_V'
    "$\langle \check{U} \rangle$";                              % 'AverageFluidVelocityCRIT'
    "$\check{U}_\mathrm{max}$";                                 % 'MaximumFluidVelocityCRIT'
    "$\mathrm{\check{Re}}$";                                    % 'ReynoldsNumberUsingU_avgCRIT'
    "$\mathrm{\check{Re}_{p,d_\mathrm{max}}}$";                 % 'ParticleReynoldsNumberUsingU_avgD_maxCRIT'
    "$\mathrm{\check{Re}_{p,d_\mathrm{V}}}$";                   % 'ParticleReynoldsNumberUsingU_avgD_VCRIT'
    "$L_\mathrm{e_{Truskey}}$";                                 % 'EntranceLengthTruskeyRectangularDuct'
    "$L_\mathrm{e_{Bergman}}$";                                 % 'EntranceLengthBergmanPipe'
    "$t_\mathrm{e}$";                                           % 'TimeToSteadyState'
    "$\langle \dot{\gamma}_\mathrm{w} \rangle$";                % 'AvgWallShearRateCRIT'
    "$\mathrm{\check{Re}}_{\mathrm{s},d_\mathrm{max}}$";        % 'ShearReynoldsNumberUsingD_maxCRIT'
    "$\mathrm{\check{Re}}_{\mathrm{s},d_\mathrm{V}}$";          % 'ShearReynoldsNumberUsingD_VCRIT'
    "$\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{max}}^*$";      % 'ModifiedParticleReynoldsNumberUsingU_avgD_maxCRIT'
    "$\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}^*$";        % 'ModifiedParticleReynoldsNumberUsingU_avgD_VCRIT'    
    "$\mathrm{Ar}_{d_\mathrm{V}}^*$";                           % 'ModifiedArchimedesNumberUsingD_V'
    "$\mathrm{Fr}_{d_\mathrm{max}}$"];                          % 'FroudeNumberUsingU_avgD_maxCRIT'

paramLists.str = ... %[40 count]
    ["Channel height";                                                                          % 'ChannelHeight'
    "Channel width";                                                                            % 'ChannelWidth'
    "Solvent";                                                                                  % 'Solvent'
    "Volumetric flow rate";                                                                     % 'FluidFlowRateCRIT'
    "Channel length";                                                                           % 'ChannelLength'
    "Particle rotational diameter";                                                             % 'RotationalDiameter'
    "Equivalent volume spherical particle diameter";                                            % 'EquivVolSpherDiameter'
    "Hydraulic diameter";                                                                       % 'HydraulicDiameter'
    "Channel aspect ratio";                                                                     % 'ChannelAspectRatio'
    "Ratio of rotational diameter and channel width";                                           % 'RatioD_maxChannelWidth'
    "Ratio of channel width and rotational diameter";                                           % 'RatioChannelWidthD_max'
    "Ratio of rotational diameter and channel height";                                          % 'RatioD_maxChannelHeight'
    "Ratio of channel height and rotational diameter";                                          % 'RatioChannelHeightD_max'
    "Ratio of rotational diameter squared to channel height and channel width";                 % 'RatioD_maxSqandChannelHeightChannelWidth'
    "Ratio of channel height and channel width to rotational diameter squared";                 % 'RatioChannelHeightChannelWidthandD_maxSq'
    "Ratio of equivalent volume spherical diameter and channel width";                          % 'RatioD_VChannelWidth'
    "Ratio of channel width and equivalent volume spherical diameter";                          % 'RatioChannelWidthD_V'
    "Ratio of equivalent volume spherical diameter and channel height";                         % 'RatioD_VChannelHeight'
    "Ratio of channel height and equivalent volume spherical diameter";                         % 'RatioChannelHeightD_V'
    "Ratio of equivalent volume spherical diameter squared to channel height and channel width";% 'RatioD_VSqandChannelHeightChannelWidth'
    "Ratio of channel height and channel width to equivalent volume spherical diameter squared";% 'RatioChannelHeightChannelWidthandD_VSq'
    "Gravity force";                                                                            % 'GravityForce'
    "Buoyant force";                                                                            % 'BuoyantForce'
    "Archimedes number";                                                                        % 'ArchimedesNumberUsingD_max'
    "Archimedes number based on equivalent volume spherical diameter";                          % 'ArchimedesNumberUsingD_V'
    "Average fluid velocity";                                                                   % 'AverageFluidVelocityCRIT'
    "Maximum fluid velocity"                                                                    % 'MaximumFluidVelocityCRIT'
    "Reynolds number";                                                                          % 'ReynoldsNumberUsingU_avgCRIT'
    "Particle Reynolds number based on rotational diameter";                                    % 'ParticleReynoldsNumberUsingU_avgD_maxCRIT'
    "Particle Reynolds bumber based on equivalent volume spherical diameter";                   % 'ParticleReynoldsNumberUsingU_avgD_VCRIT'
    "Entrance length, Truskey";                                                                 % 'EntranceLengthTruskeyRectangularDuct'
    "Entrance length, Bergman";                                                                 % 'EntranceLengthBergmanPipe'
    "Time to reach steady-state";                                                               % 'TimeToSteadyState'
    "Average wall shear rate";                                                                  % 'AvgWallShearRateCRIT'
    "Critical shear Reynolds number";                                                           % 'ShearReynoldsNumberUsingD_maxCRIT'
    "Critical shear Reynolds number based on equivalent volume spherical diameter";             % 'ShearReynoldsNumberUsingD_VCRIT'
    "Modified critical particle Reynolds number";                                               % 'ModifiedParticleReynoldsNumberUsingU_avgD_maxCRIT'
    "Modified critical particle Reynolds number";                                               % 'ModifiedParticleReynoldsNumberUsingU_avgD_VCRIT'    
    "Modified Archimedes number";                                                               % 'ModifiedArchimedesNumberUsingD_V'
    "Froude number based on rotational diameter"];                                              % 'FroudeNumberUsingU_avgD_maxCRIT'

paramLists.dim = ... %[40 count]
    ["(m)";                                                                 % 'ChannelHeight'
    "(m)";                                                                  % 'ChannelWidth'
    "";                                                                     % 'Solvent'
    "(m^3/s)";                                                              % 'FluidFlowRateCRIT'
    "(m)";                                                                  % 'ChannelLength'
    "(m)";                                                                  % 'RotationalDiameter'
    "(m)";                                                                  % 'EquivVolSpherDiameter'
    "(m)";                                                                  % 'HydraulicDiameter'
    "";                                                                     % 'ChannelAspectRatio'
    "";                                                                     % 'RatioD_maxChannelWidth'
    "";                                                                     % 'RatioChannelWidthD_max'
    "";                                                                     % 'RatioD_maxChannelHeight'
    "";                                                                     % 'RatioChannelHeightD_max'
    "";                                                                     % 'RatioD_maxSqandChannelHeightChannelWidth'
    "";                                                                     % 'RatioChannelHeightChannelWidthandD_maxSq'
    "";                                                                     % 'RatioD_VChannelWidth'
    "";                                                                     % 'RatioChannelWidthD_V'
    "";                                                                     % 'RatioD_VChannelHeight'
    "";                                                                     % 'RatioChannelHeightD_V'
    "";                                                                     % 'RatioD_VSqandChannelHeightChannelWidth'
    "";                                                                     % 'RatioChannelHeightChannelWidthandD_VSq'
    "(N)";                                                                  % 'GravityForce'
    "(N)";                                                                  % 'BuoyantForce'
    "";                                                                     % 'ArchimedesNumberUsingD_max'
    "";                                                                     % 'ArchimedesNumberUsingD_V'
    "(m/s)";                                                                % 'AverageFluidVelocityCRIT'
    "(m/s)";                                                                % 'MaximumFluidVelocityCRIT'
    "";                                                                     % 'ReynoldsNumberUsingU_avgCRIT'
    "";                                                                     % 'ParticleReynoldsNumberUsingU_avgD_maxCRIT'
    "";                                                                     % 'ParticleReynoldsNumberUsingU_avgD_VCRIT'
    "(m)";                                                                  % 'EntranceLengthTruskeyRectangularDuct'
    "(m)";                                                                  % 'EntranceLengthBergmanPipe'
    "(s)";                                                                  % 'TimeToSteadyState'
    "(s^-1)";                                                               % 'AvgWallShearRateCRIT'
    "";                                                                     % 'ShearReynoldsNumberUsingD_maxCRIT'
    "";                                                                     % 'ShearReynoldsNumberUsingD_VCRIT'
    "";                                                                     % 'ModifiedParticleReynoldsNumberUsingU_avgD_maxCRIT'
    "";                                                                     % 'ModifiedParticleReynoldsNumberUsingU_avgD_VCRIT'
    "";                                                                     % 'ModifiedArchimedesNumberUsingD_V'
    ""];                                                                    % 'FroudeNumberUsingU_avgD_maxCRIT'

% Rev unit array [$H$ $W$ $\check{Q}$ $L$ $D_{max}$ $D_{V}$ $D_{H}$ $F_{g}$ $F_{b}$ $\mathrm{Ar_{D_{max}}}$ $\mathrm{Ar_{D_{V}}}$ $\langle \check{U} \rangle$  $\check{U}_{max}$ $\mathrm{\check{Re}}$ $\mathrm{\check{Re}_{p,D_{max}}}$ $\mathrm{\check{Re}_{p,D_{V}}}$ $L_{e_{Truskey}}$ $L_{e_{Bergman}}$ $t_{e}$ $\langle \check {\dot{\gamma}}_{w}\rangle$ $\mathrm{\check{R}_{D_{max}}}$ $\mathrm{\check{R}_{D_{V}}}$ $\mathrm{\check{Re}}^*_{p,D_{max}}$ $\mathrm{\check{Re}}^*_{p,D_{V}}$ $\mathrm{Ar}^*_{D_{V}}$ $\mathrm{Fr_{D_{max}}}$]
paramLists.unitArray.('ChannelHeight') =                                      [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('ChannelWidth') =                                       [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('Solvent') =                                            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('FluidFlowRateCRIT') =                                  [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('ChannelLength') =                                      [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('RotationalDiameter') =                                 [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('EquivVolSpherDiameter') =                              [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('HydraulicDiameter') =                                  [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('ChannelAspectRatio') =                                 [1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('RatioD_maxChannelWidth') =                             [0 -1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('RatioChannelWidthD_max') =                             [0 1 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('RatioD_maxChannelHeight') =                            [-1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('RatioChannelHeightD_max') =                            [1 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('RatioD_maxSqandChannelHeightChannelWidth') =           [-1 -1 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('RatioChannelHeightChannelWidthandD_maxSq') =           [1 1 0 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('RatioD_VChannelWidth') =                               [0 -1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('RatioChannelWidthD_V') =                               [0 1 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('RatioD_VChannelHeight') =                              [-1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('RatioChannelHeightD_V') =                              [1 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('RatioD_VSqandChannelHeightChannelWidth') =             [-1 -1 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('RatioChannelHeightChannelWidthandD_VSq') =             [1 1 0 0 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('GravityForce') =                                       [0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('BuoyantForce') =                                       [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('ArchimedesNumberUsingD_max') =                         [0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('ArchimedesNumberUsingD_V') =                           [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('AverageFluidVelocityCRIT') =                           [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('MaximumFluidVelocityCRIT') =                           [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('ReynoldsNumberUsingU_avgCRIT') =                       [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('ParticleReynoldsNumberUsingU_avgD_maxCRIT') =          [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('ParticleReynoldsNumberUsingU_avgD_VCRIT') =            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('EntranceLengthTruskeyRectangularDuct') =               [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
paramLists.unitArray.('EntranceLengthBergmanPipe') =                          [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0];
paramLists.unitArray.('TimeToSteadyState') =                                  [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0];
paramLists.unitArray.('AvgWallShearRateCRIT') =                               [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
paramLists.unitArray.('ShearReynoldsNumberUsingD_maxCRIT') =                  [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0];
paramLists.unitArray.('ShearReynoldsNumberUsingD_VCRIT') =                    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0];
paramLists.unitArray.('ModifiedParticleReynoldsNumberUsingU_avgD_maxCRIT') =  [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]; 
paramLists.unitArray.('ModifiedParticleReynoldsNumberUsingU_avgD_VCRIT') =    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0];
paramLists.unitArray.('ModifiedArchimedesNumberUsingD_V') =                   [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0];
paramLists.unitArray.('FroudeNumberUsingU_avgD_maxCRIT') =                    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];

% Prelim Unit array [H W S Q L D_{H} AR D_{max}/W W/D_{max} F_{G} F_{B} Ar_{D_{max}} Ar_{EqvSphD}} <v> Re Re_{p_{D_{max}}} Re_{p_{EqvSphD}} L_{e_{Truskey}} L_{e_{Bergman}} t_{e} \dot{\gamma} R_{c_{D_{max}} R_{EqvSphD}]
% Unit array [H W Q L D_{max} D_{V} D_{H} F_{G} F_{B} Ar_{D_{max}} Ar_{EqvSphD}} <v> Re Re_{p_{D_{max}}} Re_{p_{EqvSphD}} L_{e_{Truskey}} L_{e_{Bergman}} t_{e} \dot{\gamma} R_{c_{D_{max}} R_{EqvSphD}]
% Rev unit array [$H$ $W$ $\check{Q}$ $L$ $D_{max}$ $D_{V} $D_{H}$ $F_{g}$ $F_{b}$ $\mathrm{Ar_{D_{max}}}$  $\mathrm{Ar_{D_{V}}}$ $\langle \check{U} \rangle$  $\check{U}_{max}$ $\mathrm{\check{Re}}$ $\mathrm{\check{Re}_{p,D_{max}}}$ $\mathrm{\check{Re}_{p,D_{V}}}$ $L_{e_{Truskey}}$ $L_{e_{Bergman}}$ $t_{e}$ $\langle \check {\dot{\gamma}}_{w}\rangle$ $\mathrm{\check{R}_{D_{max}}}$ $\mathrm{\check{R}_{D_{V}}}$ $\mathrm{\check{Re}_{p,D_{max}}^*}$ $\mathrm{\check{Re}_{p,D_{V}}^*}$ $\mathrm Ar_{D_{V}}^*$ $\mathrm{Fr_{D_{max}}}$]
paramLists.unitDef = ...
    ["$H$";                                                     % 'ChannelHeight'
    "$W$";                                                      % 'ChannelWidth'
    "$\check{Q}$";                                              % 'FluidFlowRateCRIT'
    "$L$";                                                      % 'ChannelLength'
    "$d_\mathrm{max}$";                                         % 'RotationalDiameter'
    "$d_\mathrm{V}$";                                           % 'EquivVolSpherDiameter'
    "$D_{H}$";                                                  % 'HydraulicDiameter'
    "$F_\mathrm{g}$";                                           % 'GravityForce'
    "$F_\mathrm{b}$";                                           % 'BuoyantForce'
    "$\mathrm{Ar}_{d_\mathrm{max}}$";                           % 'ArchimedesNumberUsingD_max'
    "$\mathrm{Ar}_{d_\mathrm{V}}$";                             % 'ArchimedesNumberUsingD_V'
    "$\langle \check{U} \rangle$";                              % 'AverageFluidVelocityCRIT'
    "$\check{U}_\mathrm{max}$";                                 % 'MaximumFluidVelocityCRIT'
    "$\mathrm{\check{Re}}$";                                    % 'ReynoldsNumberUsingU_avgCRIT'
    "$\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{max}}}$";       % 'ParticleReynoldsNumberUsingU_avgD_maxCRIT'
    "$\mathrm{\check{Re}_{\mathrm{p},d_\mathrm{V}}}$";          % 'ParticleReynoldsNumberUsingU_avgD_VCRIT'
    "$L_\mathrm{e_{Truskey}}$";                                 % 'EntranceLengthTruskeyRectangularDuct'
    "$L_\mathrm{e_{Bergman}}$";                                 % 'EntranceLengthBergmanPipe'
    "$t_\mathrm{e}$";                                           % 'TimeToSteadyState'
    "$\langle \dot{\gamma}_\mathrm{w} \rangle$";                % 'AvgWallShearRateCRIT'
    "$\mathrm{\check{Re}}_{\mathrm{s},d_\mathrm{max}}$";        % 'ShearReynoldsNumberUsingD_maxCRIT'
    "$\mathrm{\check{Re}}_{\mathrm{s},d_\mathrm{V}}$";          % 'ShearReynoldsNumberUsingD_VCRIT'
    "$\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{max}}^*$";      % 'ModifiedParticleReynoldsNumberUsingU_avgD_maxCRIT'
    "$\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}^*$";        % 'ModifiedParticleReynoldsNumberUsingU_avgD_VCRIT'    
    "$\mathrm{Ar}_{d_\mathrm{V}}^*$";                           % 'ModifiedArchimedesNumberUsingD_V'
    "$\mathrm{Fr}_{d_\mathrm{max}}$"];                          % 'FroudeNumberUsingU_avgD_maxCRIT'

paramLists.unitStr = ... % currently not used
    ["Channel height (m)";                                                                          % 'ChannelHeight'
    "Channel width (m)";                                                                            % 'ChannelWidth'
    "Volumetric flow rate (m^3/s)";                                                                 % 'FluidFlowRateCRIT'
    "Channel length (m)";                                                                           % 'ChannelLength'
    "Particle rotational diameter (m)";                                                             % 'RotationalDiameter'
    "Equivalent volume spherical particle diameter (m)";                                            % 'EquivVolSpherDiameter'
    "Hydraulic diameter (m)";                                                                       % 'HydraulicDiameter'
    "Gravity force (N)";                                                                            % 'GravityForce'
    "Buoyant force (N)";                                                                            % 'BuoyantForce'
    "Archimedes number (-)";                                                                        % 'ArchimedesNumberUsingD_max'
    "Archimedes number based on equivalent volume spherical diameter (-)";                          % 'ArchimedesNumberUsingD_V'
    "Average fluid velocity (m/s)";                                                                 % 'AverageFluidVelocityCRIT'
    "Maximum fluid velocity (m/s)";                                                                 % 'MaximumFluidVelocityCRIT'
    "Reynolds number (-)";                                                                          % 'ReynoldsNumberUsingU_avgCRIT'
    "Particle Reynolds number based on rotational diameter (-)";                                    % 'ParticleReynoldsNumberUsingU_avgD_maxCRIT'
    "Particle Reynolds number based on equivalent volume spherical diameter (-)";                   % 'ParticleReynoldsNumberUsingU_avgD_VCRIT'
    "Entrance length, Truskey (m)";                                                                 % 'EntranceLengthTruskeyRectangularDuct'
    "Entrance length, Bergman (m)";                                                                 % 'EntranceLengthBergmanPipe'
    "Time to reach steady-state (s)";                                                               % 'TimeToSteadyState'
    "Average wall shear rate (s^-1)";                                                               % 'AvgWallShearRateCRIT'
    "Critical shear Reynolds number (-)";                                                           % 'ShearReynoldsNumberUsingD_maxCRIT'
    "Critical shear Reynolds number based on equivalent volume spherical diameter (-)";             % 'ShearReynoldsNumberUsingD_VCRIT'
    "Modified critical particle Reynolds number (-)";                                               % 'ModifiedParticleReynoldsNumberUsingU_avgD_maxCRIT'
    "Modified critical particle Reynolds number (-)";                                               % 'ModifiedParticleReynoldsNumberUsingU_avgD_VCRIT'
    "Modified Archimedes number (-)";                                                               % 'ModifiedArchimedesNumberUsingD_V'
    "Froude number based on rotational diameter (-)"];                                              % 'FroudeNumberUsingU_avgD_maxCRIT'

paramLists.unitPropertyStr = ...
    ["ChannelHeight";
    "ChannelWidth";
    "FluidFlowRateCRIT";
    "ChannelLength";
    "RotationalDiameter";
    "EquivVolSpherDiameter";
    "HydraulicDiameter";
    "GravityForce";
    "BuoyantForce";
    "ArchimedesNumberUsingD_max";
    "ArchimedesNumberUsingD_V";
    "AverageFluidVelocityCRIT";
    "MaximumFluidVelocityCRIT";
    "ReynoldsNumberUsingU_avgCRIT";
    "ParticleReynoldsNumberUsingU_avgD_maxCRIT";
    "ParticleReynoldsNumberUsingU_avgD_VCRIT";
    "EntranceLengthTruskeyRectangularDuct";
    "EntranceLengthBergmanPipe";
    "TimeToSteadyState";
    "AvgWallShearRateCRIT";
    "ShearReynoldsNumberUsingD_maxCRIT";
    "ShearReynoldsNumberUsingD_VCRIT";
    "ModifiedParticleReynoldsNumberUsingU_avgD_maxCRIT";
    "ModifiedParticleReynoldsNumberUsingU_avgD_VCRIT";
    "ModifiedArchimedesNumberUsingD_V";
    "FroudeNumberUsingU_avgD_maxCRIT"];

%% Generate table of parameters
%---[Preallocate boolean array to flag which Parameters Property will have an associated uncertainty]---
parametersBooleanExtra = NaN(size(parametersProperties,1),size(parametersProperties,2));

%---[Evaluate boolean array]---
for q = 1:length(parametersProperties)
    parentPropertyClassBool = Parameters{1,1}.(parametersProperties{q,1});
    if isa(parentPropertyClassBool,'DimVar') || isa(parentPropertyClassBool,'UC')
        parametersBooleanExtra(q,1) = 1;
    else
        parametersBooleanExtra(q,1) = 0;
    end
end

%---[Preallocate string array]---
tableVarNames = strings(length(parametersProperties)+nnz(parametersBooleanExtra),1);

%---[Evaluate table variable names]---
countTableVarNames = 0;
for r = 1:length(parametersProperties)
    countTableVarNames = countTableVarNames + 1;
    tableVarNames(countTableVarNames) = parametersProperties{r,1};
    tableVarTypes(countTableVarNames) = {'double'};
    if parametersBooleanExtra(r,1) == 1
        countTableVarNames = countTableVarNames + 1;
        addTableVarNamesStr = strcat(parametersProperties{r,1},"Uncertainty");
        tableVarNames(countTableVarNames) = addTableVarNamesStr;
        tableVarTypes(countTableVarNames) = {'double'};
    else
        tableVarTypes(countTableVarNames) = {'string'};
    end
end

sz = [length(Parameters) length(tableVarNames)];

channelHeightStr = arrayfun(@(a) sprintf('H%0.3f',a),ChannelHeight, 'UniformOutput', false);
channelWidthStr = arrayfun(@(b) sprintf('W%0.3f',b),ChannelWidth, 'UniformOutput', false);
solventChar = cellfun(@(c) convertStringsToChars(c),Solvent, 'UniformOutput', false);
critFlowRateStr = arrayfun(@(d) sprintf('cfr%0.3f',d),CritFlowRate, 'UniformOutput', false);

tableRowNames = cellfun(@(e,f,g,h) strcat(e," ",f," ",g," ",h),channelHeightStr,channelWidthStr,Solvent,critFlowRateStr,'UniformOutput', false);

critLiftOffTable = table('Size',sz,'VariableTypes',tableVarTypes,'VariableNames',tableVarNames);

%---[Fill table with parameters data]---
countCritLiftOffTable = 0;

for t = 1:length(Parameters)
    for u = 1:length(parametersProperties)
        countCritLiftOffTable = countCritLiftOffTable + 1;
        parentPropertyClassTable = Parameters{1,1}.(parametersProperties{u,1});
        if isa(parentPropertyClassTable,'DimVar')
            critLiftOffTable(t,countCritLiftOffTable) = {Parameters{t,1}.(parametersProperties{u,1}).Value.Value};
        elseif isa(parentPropertyClassTable,'UC')
            critLiftOffTable(t,countCritLiftOffTable) = {Parameters{t,1}.(parametersProperties{u,1}).Value};
        else
            critLiftOffTable(t,countCritLiftOffTable) = {Parameters{t,1}.(parametersProperties{u,1})};
        end

        if parametersBooleanExtra(u,1) == 1
            countCritLiftOffTable = countCritLiftOffTable + 1;
            if isa(parentPropertyClassTable,'DimVar')
                critLiftOffTable(t,countCritLiftOffTable) = {Parameters{t,1}.(parametersProperties{u,1}).Value.Err};
            else % isa(parentPropertyClassTable,'UC')
                critLiftOffTable(t,countCritLiftOffTable) = {Parameters{t,1}.(parametersProperties{u,1}).Err};
            end
        end
    end
    countCritLiftOffTable = 0;
end

%% Find duplicate center point data
%---[Find rows with more than one occurrence, but cannot distinguish which rows correspond to which set of data]---
[~, ia, ic] = unique([critLiftOffTable.ChannelHeight,critLiftOffTable.ChannelWidth,critLiftOffTable.Solvent], 'rows');
duplrowidx = setdiff(1:size([critLiftOffTable.ChannelHeight,critLiftOffTable.ChannelWidth,critLiftOffTable.Solvent],1), ia( sum(bsxfun(@eq,ic,(1:max(ic))))<=1 ));

%---[Find which row pairs are duplicates]---
duplZeros = zeros(length(duplrowidx),2);
for a = 1:length(duplrowidx)
    duplZeros(a,:) = find(ic == ic(duplrowidx(a)));
end

%---[Remove doubles of row pairs containing duplicates]---
[U,I,J] = unique(duplZeros, 'rows', 'first');
hasDuplicates = size(U,1) < size(duplZeros,1);
ixDuplRows = setdiff(1:size(duplZeros,1), I);
duplRowValues = duplZeros(ixDuplRows,:); %each row is a pair of duplicate center points data

%% Combine duplicate center point data
%---[Create a new cell array to store averaged objects]---
ParametersAvg = cell(length(duplRowValues), 1);

for w = 1:length(duplRowValues)
    %---[Compute averaged properties]---
    avgChannelHeight = (Parameters{duplRowValues(w,1)}.ChannelHeight + ...
                        Parameters{duplRowValues(w,2)}.ChannelHeight) / 2;
                    
    avgChannelWidth  = (Parameters{duplRowValues(w,1)}.ChannelWidth + ...
                        Parameters{duplRowValues(w,2)}.ChannelWidth) / 2;
                    
    solvent = Parameters{duplRowValues(w,1)}.Solvent;

    avgCritFlowRate  = (Parameters{duplRowValues(w,1)}.CritFlowRate + ...
                        Parameters{duplRowValues(w,2)}.CritFlowRate) / 2;

    ParametersAvg{w,1} = B2KCalc.B2KAllParametersLiftOff(avgChannelHeight,avgChannelWidth,solvent,avgCritFlowRate);
end

%% Reshape and resize Parameters cell and critLiftOffTable
% Start with a copy of the original cell array
ParametersCombined = Parameters;

% Replace the entries corresponding to the first run of each duplicate pair with the newly averaged objects
for y = 1:length(duplRowValues)
    ParametersCombined(duplRowValues(y,1)) = ParametersAvg(y);
end

% Delete the entries corresponding to the second run
for z = length(duplRowValues):-1:1
    ParametersCombined(duplRowValues(z,2)) = [];
end

%% Generate combined table of parameters
% Evaluate table variable names
szCombined = [length(ParametersCombined) length(tableVarNames)];

% Make copies
ChannelHeightCombined = ChannelHeight;
ChannelWidthCombined = ChannelWidth;
SolventCombined = Solvent;
CritFlowRateCombined = CritFlowRate;

for rr = length(duplRowValues):-1:1
    % Replace rows corresponding to first run with combined data
    CritFlowRateCombined(duplRowValues(rr,1),1) = (CritFlowRate(duplRowValues(rr,1),1) + CritFlowRate(duplRowValues(rr,2),1)) / 2;

    % Delete rows corresponding to second run
    ChannelHeightCombined(duplRowValues(rr,2)) = [];
    ChannelWidthCombined(duplRowValues(rr,2)) = []; 
    SolventCombined(duplRowValues(rr,2)) = [];
    CritFlowRateCombined(duplRowValues(rr,2)) = [];
end

channelHeightStrCombined = arrayfun(@(aa) sprintf('H%0.3f',aa),ChannelHeightCombined, 'UniformOutput', false);
channelWidthStrCombined = arrayfun(@(bb) sprintf('W%0.3f',bb),ChannelWidthCombined, 'UniformOutput', false);
solventCharCombined = cellfun(@(cc) convertStringsToChars(cc),SolventCombined, 'UniformOutput', false);
critFlowRateStrCombined = arrayfun(@(dd) sprintf('cfr%0.3f',dd),CritFlowRateCombined, 'UniformOutput', false);

tableRowNamesCombined = cellfun(@(ee,ff,gg,hh) strcat(ee," ",ff," ",gg," ",hh),channelHeightStrCombined,channelWidthStrCombined,SolventCombined,critFlowRateStrCombined,'UniformOutput', false);

critLiftOffTableCombined = table('Size',szCombined,'VariableTypes',tableVarTypes,'VariableNames',tableVarNames);

% Fill table with parameters data
countCritLiftOffTableCombined = 0;

for tt = 1:length(ParametersCombined)
    for uu = 1:length(parametersProperties)
        countCritLiftOffTableCombined = countCritLiftOffTableCombined + 1;
        parentPropertyClassTableCombined = ParametersCombined{1,1}.(parametersProperties{uu,1});
        if isa(parentPropertyClassTableCombined,'DimVar')
            critLiftOffTableCombined(tt,countCritLiftOffTableCombined) = {ParametersCombined{tt,1}.(parametersProperties{uu,1}).Value.Value};
        elseif isa(parentPropertyClassTableCombined,'UC')
            critLiftOffTableCombined(tt,countCritLiftOffTableCombined) = {ParametersCombined{tt,1}.(parametersProperties{uu,1}).Value};
        else
            critLiftOffTableCombined(tt,countCritLiftOffTableCombined) = {ParametersCombined{tt,1}.(parametersProperties{uu,1})};
        end

        if parametersBooleanExtra(uu,1) == 1
            countCritLiftOffTableCombined = countCritLiftOffTableCombined + 1;
            if isa(parentPropertyClassTableCombined,'DimVar')
                critLiftOffTableCombined(tt,countCritLiftOffTableCombined) = {ParametersCombined{tt,1}.(parametersProperties{uu,1}).Value.Err};
            else % isa(parentPropertyClassTable,'UC')
                critLiftOffTableCombined(tt,countCritLiftOffTableCombined) = {ParametersCombined{tt,1}.(parametersProperties{uu,1}).Err};
            end
        end
    end
    countCritLiftOffTableCombined = 0;
end

%% Evaluate (n)Parameters in paramStruct
sortingParamName = 'zParameterHW_D_maxSq';
% sortingParamName = 'zParameterHW_D_VSq';
[paramStruct,refSize,refSizeEach] = evaluateParameters(flagCombine, Parameters, ParametersCombined, dataCase, paramStruct, paramLists, sortingParamName);

%% lineSpec for arrayMarker and arrayColor
lineSpec = [paramStruct.plotData.arrayMarker, paramStruct.plotData.arrayColor];

%% Figure 2) 2D plot of critical shear Re number vs Archimedes number vs. HW/(D_{max})^2
%---[Figure and Axis Handles]---
fig2 = figure(2);
ax1_fig2 = gca;

% Custom axis and legend font size
AxisFontSizeMultiplier_fig2 = 6/4.9; %1
AxisLabelFontSizeMultiplier_fig2 = 9/7; %1
LegendLabelFontSizeMultiplier_fig2 = 6/4.9; %1
LegendFontSizeMultiplier_fig2 = 6/4.9; %1

% Custom figure size scale factors (to account for limitations in exported fig size using exportgraphics MATLAB R2024b)
expGraphicsScaleFactor_fig2 = (14/11.8533).*(14/13.8289);

%---[Get Figure Properties]---
figProps = getDefaultFigProperties('1.5Column',AxisFontSizeMultiplier_fig2,AxisLabelFontSizeMultiplier_fig2,LegendLabelFontSizeMultiplier_fig2,LegendFontSizeMultiplier_fig2);

%---[Set Figure Size]---
set(fig2,'Units','centimeters','Position',[figProps.Figure.Position.Left figProps.Figure.Position.Bottom figProps.Figure.Position.Width.*expGraphicsScaleFactor_fig2 figProps.Figure.Position.Width*figProps.Figure.Position.HWRatio.*expGraphicsScaleFactor_fig2])
pos = get(fig2,'Position');
set(fig2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%
lineArchimedesNumValue = [10071 1695.3  25020];
xl = xline(lineArchimedesNumValue,'--',{'\textbf{10070}' '\textbf{1700}' '\textbf{25000}'},'LabelVerticalAlignment','top','Color',"#117ec2",'Interpreter','latex');
set(xl,'FontSize',figProps.Legend.Fonts.FontSize)
hold on

%% Set figure properties

%---[Axis properties]---
%----[Font]----
ax1_fig2.FontSize = figProps.Axis.Font.FontSize;
%----[Ticks]----
ax1_fig2.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;
%----[Rulers]----
xlim([0.8 55000])
ylim([0.8 5000])

ax1_fig2.XAxis.Exponent = 0; %Scientific notation
xtickformat('%.0f') %Scientific notation
ax1_fig2.YAxis.Exponent = 0; %Scientific notation
ytickformat('%.0f') %Scientific notation
set(gca,'XScale','log');
set(gca,'YScale','log');
%----[Grids]----
grid on
ax1_fig2.GridLineStyle = figProps.Axis.Grids.GridLineStyle;
ax1_fig2.GridLineWidth = figProps.Axis.Grids.GridLineWidth;
ax1_fig2.GridColor = figProps.Axis.Grids.GridColor;
ax1_fig2.MinorGridLineStyle = figProps.Axis.Grids.MinorGridLineStyle;
ax1_fig2.MinorGridLineWidth = figProps.Axis.Grids.MinorGridLineWidth;
ax1_fig2.MinorGridColor = figProps.Axis.Grids.MinorGridColor;
%----[Labels]----
xlabel(paramStruct.xParameter.NewStr,'Interpreter',figProps.Axis.Labels.LabelInterpreter)

ax1_fig2.XLabel.FontSize = figProps.Axis.Labels.LabelFontSize;
ax1_fig2.XLabel.FontWeight = figProps.Axis.Labels.LabelFontWeight;

ylabel(paramStruct.yParameter.NewStr,'Interpreter',figProps.Axis.Labels.LabelInterpreter)

ax1_fig2.YLabel.FontSize = figProps.Axis.Labels.LabelFontSize;
ax1_fig2.YLabel.FontWeight = figProps.Axis.Labels.LabelFontWeight;

%% Drawing original points as ghost markers
for m = 1:refSize
    pOriginal(m) = scatter(paramStruct.xParameter.ValueSorted(m),paramStruct.yParameter.ValueSorted(m),figProps.Scatter.Markers.SizeData);
    set(pOriginal(m),{'Marker','MarkerFaceColor','MarkerFaceAlpha'},{lineSpec(m,1),lineSpec(m,2),figProps.Scatter.Markers.MarkerFaceAlpha})
    hold on
end
set(pOriginal,'MarkerEdgeColor',figProps.Scatter.Markers.MarkerEdgeColor)

%% Plotting and comparing Patankar 2001 data

pnew1 = scatter(pH_d48x,pH_d48y);
set(pnew1,{'Marker','MarkerFaceColor'},{"^",figProps.Scatter.Markers.MarkerEdgeColor2})
set(pnew1,'MarkerEdgeColor',"#5e4fa2")
set(pnew1,'SizeData',figProps.Scatter.Markers.SizeData2)
set(pnew1,'LineWidth',figProps.Scatter.Markers.LineWidth2/2)
pfit1 = polyfit(log(pH_d48x),log(pH_d48y),1);
y1 = polyval(pfit1,log(pH_d48x));
pfitline1 = plot(pH_d48x,exp(y1),'-');
set(pfitline1,'Color',"#5e4fa2")

pnew2 = scatter(pH_d12x,pH_d12y);
set(pnew2,{'Marker','MarkerFaceColor'},{">",figProps.Scatter.Markers.MarkerEdgeColor2})
set(pnew2,'MarkerEdgeColor',"#66c2a5")
set(pnew2,'SizeData',figProps.Scatter.Markers.SizeData2)
set(pnew2,'LineWidth',figProps.Scatter.Markers.LineWidth2/2)
pfit2 = polyfit(log(pH_d12x),log(pH_d12y),1);
y2 = polyval(pfit2,log(pH_d12x));
pfitline2 = plot(pH_d12x,exp(y2),'-');
set(pfitline2,'Color',"#66c2a5")

pnew3 = scatter(pH_d6x,pH_d6y);
set(pnew3,{'Marker','MarkerFaceColor'},{"v",figProps.Scatter.Markers.MarkerEdgeColor2})
set(pnew3,'MarkerEdgeColor',"#fdae61")
set(pnew3,'SizeData',figProps.Scatter.Markers.SizeData2)
set(pnew3,'LineWidth',figProps.Scatter.Markers.LineWidth2/2)
pfit3 = polyfit(log(pH_d6x),log(pH_d6y),1);
y3 = polyval(pfit3,log(pH_d6x));
pfitline3 = plot(pH_d6x,exp(y3),'-');
set(pfitline3,'Color',"#fdae61")

pnew4 = scatter(pH_d4x,pH_d4y);
set(pnew4,{'Marker','MarkerFaceColor'},{"<",figProps.Scatter.Markers.MarkerEdgeColor2})
set(pnew4,'MarkerEdgeColor',"#9e0142")
set(pnew4,'SizeData',figProps.Scatter.Markers.SizeData2)
set(pnew4,'LineWidth',figProps.Scatter.Markers.LineWidth2/2)
pfit4 = polyfit(log(pH_d4x),log(pH_d4y),1);
y4 = polyval(pfit4,log(pH_d4x));
pfitline4 = plot(pH_d4x,exp(y4),'-');
set(pfitline4,'Color',"#9e0142")

%%
uistack(pOriginal,'up',30)

%% Legend parameters
% Specify location of legend (northwest; start from bottom up) % north, south, east, west, northeast, northwest, southeast, southwest
loc_leg1_fig2 = 'northwest';

% Define offet and spacing for legends
leg1_axis_offset = 0.2;
leg1_xSpacing_centimeters = 0.2;

%% Zero legend [Title: Flat-plate, experimental] (unpositioned)
%---Get axis position ('centimeters')---
set(ax1_fig2,'Units','centimeters')
drawnow
get_ax1_fig2_pos = get(ax1_fig2,'Position'); %[figtom width height]
drawnow
top_pos_centimeter_ax1_fig2 = get_ax1_fig2_pos(2) + get_ax1_fig2_pos(4);
drawnow
left_pos_centimeter_ax1_fig2 = get_ax1_fig2_pos(1);
drawnow
right_pos_centimeter_ax1_fig2 = get_ax1_fig2_pos(1) + get_ax1_fig2_pos(3);
drawnow
bottom_pos_centimeter_ax1_fig2 = get_ax1_fig2_pos(2);

% Add initial features(unpositioned; 'Units','centimeters')
an_fig2 = annotation('textbox','String','Flat-plate, experimental','FontSize',figProps.Legend.Labels.Title.FontSize,'VerticalAlignment','bottom','HorizontalAlignment','center','BackgroundColor',[1 1 1],'EdgeColor',[1 1 1],'FaceAlpha',1,'Interpreter','latex','Position',[0.5 0.5 0.5 0.5],'FitBoxToText','on');
drawnow
set(an_fig2,'Units','centimeters');
drawnow
get_pos_an_fig2 = get(an_fig2,'Position'); %[x_begin y_begin length height]

anTop_fig2 = annotation('line',[0.5, 0.75], [0.95, 0.95], 'Color', 'black'); %dummy line
drawnow
set(anTop_fig2,'Units','centimeters');
drawnow
get_pos_anTop_fig2 = get(anTop_fig2,'Position'); %[x_begin y_begin dx dy]

anBottom_fig2 = annotation('line',[0.5, 0.75], [0.90, 0.90], 'Color', 'black'); %dummy line
drawnow
set(anBottom_fig2,'Units','centimeters');
drawnow
get_pos_anBottom_fig2 = get(anBottom_fig2,'Position'); %[x_begin y_begin dx dy]

%% First legend [HW/(D_{max})^2, set progressive color order]

%---Legends strings - Dynamic---ds
zLabelStrRounded = B2KFormatValuesToStrOneSigFigMaxUnc(paramStruct.zParameterHW_D_maxSq.ValueSorted,paramStruct.zParameterHW_D_maxSq.UncertaintySorted);

% Create dummy plot to generate first legend since default legend for scatter cannot modify legend marker size
for m = 1:refSize
    pDummy(m) = plot(NaN,NaN,'MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none');
    set(pDummy(m),{'Marker','MarkerFaceColor'},{lineSpec(m,1),[lineSpec(m,2)]})
    hold on
end
set(pDummy,'MarkerEdgeColor',figProps.Line.Markers.MarkerEdgeColor)

% Plot first legend
lh1_fig2 = legend(gca,pDummy(1+refSizeEach:refSizeEach*2),zLabelStrRounded,'FontSize',figProps.Legend.Fonts.FontSize,'NumColumns',1,'Interpreter',figProps.Legend.Interpreter); %circle marker
lh1_fig2.Units = 'centimeters';
lh1_fig2.Location = 'northwest';
drawnow
get_pos_lh1_fig2 = get(lh1_fig2,'Position');
lh1title_fig2 = title(lh1_fig2,paramStruct.zParameterHW_D_maxSq.LegendStr,'FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
drawnow
get_pos_lh1andtitle_fig2 = get(lh1_fig2,'Position'); %[left bottom width height]
get_pos_lh1title_width_fig2 = get_pos_lh1andtitle_fig2(3);
get_pos_lh1title_height_fig2 = get_pos_lh1andtitle_fig2(4) - get_pos_lh1_fig2(4); 

lh1_fig2.EdgeColor = figProps.Legend.ColorAndStyling.EdgeColor;

% Adjust horizontal spacing for legend lh1
lh1_fig2.ItemTokenSize(1) = figProps.Legend.Fonts.ItemTokenSize;

% Set position of first legend
lh1titleTop_fig2 = annotation('line',[0.5, 0.75], [0.85, 0.85], 'Color', 'black');
drawnow
set(lh1titleTop_fig2,'Units','centimeters');
drawnow
get_pos_lh1titleTop_fig2 = get(lh1titleTop_fig2,'Position'); %[x_begin y_begin dx dy]

lh1titleBottom_fig2 = annotation('line',[0.5, 0.75], [0.80, 0.80], 'Color', 'black');
drawnow
set(lh1titleBottom_fig2,'Units','centimeters');
drawnow
get_pos_lh1titleBottom_fig2 = get(lh1titleBottom_fig2,'Position'); %[x_begin y_begin dx dy]

lh1titleBottom2_fig2 = annotation('line',[0.5, 0.75], [0.75, 0.75], 'Color', 'black');
drawnow
set(lh1titleBottom2_fig2,'Units','centimeters');
drawnow
get_pos_lh1titleBottom2_fig2 = get(lh1titleBottom2_fig2,'Position'); %[x_begin y_begin dx dy]

%% Second legend [H/D, dependent color order]
% Copy the axes and plot the second legend
ax2_fig2 = copyobj(ax1_fig2,gcf); %Copies ah1 object (but duplicates data points)
delete(get(ax2_fig2,'Children')) %Deletes duplicate data points

% Remove second axis visibility 
set(ax2_fig2, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

%---Legends strings
zLabelStrRounded_lh2 = B2KFormatValuesToStrOneSigFigMaxUnc(paramStruct.hD_maxParameter.ValueSorted,paramStruct.hD_maxParameter.UncertaintySorted);

% Plot second legend
hold on
lh2dump1_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig2);
lh2dump2_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig2);
lh2dump3_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig2);
lh2dump4_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig2);
lh2dump5_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig2);
lh2dump6_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig2);
lh2dump7_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig2);
lh2dump8_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig2);
lh2dump9_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig2);
lh2dump10_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig2);
hold off

lh2_fig2 = legend(ax2_fig2,[lh2dump1_fig2 lh2dump2_fig2 lh2dump3_fig2 lh2dump4_fig2 lh2dump5_fig2 lh2dump6_fig2 lh2dump7_fig2 lh2dump8_fig2 lh2dump9_fig2 lh2dump10_fig2],zLabelStrRounded_lh2,'FontSize',figProps.Legend.Fonts.FontSize,'NumColumns',1,'Interpreter',figProps.Legend.Interpreter);
lh2_fig2.Units = 'centimeters';
lh2_fig2.Location = 'north';
drawnow
get_pos_lh2_fig2 = get(lh2_fig2,'Position');
lh2title_fig2 = title(lh2_fig2,paramStruct.hD_maxParameter.LegendStr,'FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
drawnow
get_pos_lh2andtitle_fig2 = get(lh2_fig2,'Position'); %[left bottom width height]
get_pos_lh2title_width_fig2 = get_pos_lh2andtitle_fig2(3);
get_pos_lh2title_height_fig2 = get_pos_lh2andtitle_fig2(4) - get_pos_lh2_fig2(4); 

lh2_fig2.EdgeColor = figProps.Legend.ColorAndStyling.EdgeColor;

% Set background color of second legend to white (instead of grey by default)
set(lh2_fig2,'Color',get(ax1_fig2,'Color'));

% Adjust horizontal spacing for legend lh2
lh2_fig2.ItemTokenSize(1) = 0; %figProps.Legend.Fonts.ItemTokenSize

lh2titleTop_fig2 = annotation('line',[0.5, 0.75], [0.70, 0.70], 'Color', 'black');
drawnow
set(lh2titleTop_fig2,'Units','centimeters');
drawnow
get_pos_lh2titleTop_fig2 = get(lh2titleTop_fig2,'Position');

lh2titleBottom_fig2 = annotation('line',[0.5, 0.75], [0.65, 0.65], 'Color', 'black');
drawnow
set(lh2titleBottom_fig2,'Units','centimeters');
drawnow
get_pos_lh2titleBottom_fig2 = get(lh2titleBottom_fig2,'Position');

lh2titleBottom2_fig2 = annotation('line',[0.5, 0.75], [0.60, 0.60], 'Color', 'black');
drawnow
set(lh2titleBottom2_fig2,'Units','centimeters');
drawnow
get_pos_lh2titleBottom2_fig2 = get(lh2titleBottom2_fig2,'Position');

%% Third legend [AR] [TESTING]

% Copy the axes and plot the third legend
ax3_fig2 = copyobj(ax2_fig2,gcf); %Copies ah1 object (but duplicates data points)
delete(get(ax3_fig2,'Children')) %Deletes duplicate data points

% Remove third axis visibility 
set(ax3_fig2, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

%---Legends strings
zLabelStrRounded_lh3 = B2KFormatValuesToStrOneSigFigMaxUnc(paramStruct.aParameter.ValueSorted,paramStruct.aParameter.UncertaintySorted);

% Plot third legend
hold on
lh3dump1_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig2);
lh3dump2_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig2);
lh3dump3_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig2);
lh3dump4_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig2);
lh3dump5_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig2);
lh3dump6_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig2);
lh3dump7_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig2);
lh3dump8_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig2);
lh3dump9_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig2);
lh3dump10_fig2 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig2);
hold off

lh3_fig2 = legend(ax3_fig2,[lh3dump1_fig2 lh3dump2_fig2 lh3dump3_fig2 lh3dump4_fig2 lh3dump5_fig2 lh3dump6_fig2 lh3dump7_fig2 lh3dump8_fig2 lh3dump9_fig2 lh3dump10_fig2],zLabelStrRounded_lh3,'FontSize',figProps.Legend.Fonts.FontSize,'NumColumns',1,'Interpreter','latex');
lh3_fig2.Units = 'centimeters';
lh3_fig2.Location = 'northeast';
drawnow
get_pos_lh3_fig2 = get(lh3_fig2,'Position');
lh3title_fig2 = title(lh3_fig2,paramStruct.aParameter.LegendStr,'FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
drawnow
get_pos_lh3andtitle_fig2 = get(lh3_fig2,'Position'); %[left bottom width height]
get_pos_lh3title_width_fig2 = get_pos_lh3andtitle_fig2(3);
get_pos_lh3title_height_fig2 = get_pos_lh3andtitle_fig2(4) - get_pos_lh3_fig2(4); 

lh3_fig2.EdgeColor = figProps.Legend.ColorAndStyling.EdgeColor;

% Set background color of third legend to white (instead of grey by default)
set(lh3_fig2,'Color',get(ax1_fig2,'Color'));

% Adjust horizontal spacing for legend lh3
lh3_fig2.ItemTokenSize(1) = 0; %figProps.Legend.Fonts.ItemTokenSize

lh3titleTop_fig2 = annotation('line',[0.5, 0.75], [0.55, 0.55], 'Color', 'black');
drawnow
set(lh3titleTop_fig2,'Units','centimeters');
drawnow
get_pos_lh3titleTop_fig2 = get(lh3titleTop_fig2,'Position');

lh3titleBottom_fig2 = annotation('line',[0.5, 0.75], [0.50, 0.50], 'Color', 'black');
drawnow
set(lh3titleBottom_fig2,'Units','centimeters');
drawnow
get_pos_lh3titleBottom_fig2 = get(lh3titleBottom_fig2,'Position');


lh3titleBottom2_fig2 = annotation('line',[0.5, 0.75], [0.45, 0.45], 'Color', 'black');
drawnow
set(lh3titleBottom2_fig2,'Units','centimeters');
drawnow
get_pos_lh3titleBottom2_fig2 = get(lh3titleBottom2_fig2,'Position');

%% Fourth legend [Solvent]

% Create array for markers and labels for fourth legend
lh4_fig2_arrayMarkerdump = ["o" "square" "diamond"];
lh4_fig2_typeSolvent = cellstr(["Isopropanol","Water","Methanol"]);

% Copy the axes and plot the fourth legend
ax4_fig2 = copyobj(ax3_fig2,gcf); %Copies ah1 object (but duplicates data points)
delete(get(ax4_fig2,'Children')) %Deletes duplicate data points

% Remove fourth axis visibility 
set(ax4_fig2, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

% Plot "dummy" data for fourth legend; tried using single variable with loop but does not work
hold on
lh4dump1_fig2 = plot(NaN,'Marker',lh4_fig2_arrayMarkerdump(1),'DisplayName',lh4_fig2_typeSolvent{1},'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none','Parent',ax4_fig2);
lh4dump2_fig2 = plot(NaN,'Marker',lh4_fig2_arrayMarkerdump(2),'DisplayName',lh4_fig2_typeSolvent{2},'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none','Parent',ax4_fig2);
lh4dump3_fig2 = plot(NaN,'Marker',lh4_fig2_arrayMarkerdump(3),'DisplayName',lh4_fig2_typeSolvent{3},'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none','Parent',ax4_fig2);
hold off

lh4_fig2 = legend(ax4_fig2,[lh4dump1_fig2 lh4dump2_fig2 lh4dump3_fig2],lh4_fig2_typeSolvent,'FontSize',figProps.Legend.Fonts.FontSize,'Interpreter',figProps.Legend.Interpreter);
lh4_fig2.Units = 'centimeters';
lh4_fig2.Location = 'southwest';
drawnow
get_pos_lh4_fig2 = get(lh4_fig2,'Position');
lh4title_fig2 = title(lh4_fig2,'Solvent','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
drawnow
get_pos_lh4andtitle_fig2 = get(lh4_fig2,'Position'); %[left bottom width height]
get_pos_lh4title_width_fig2 = get_pos_lh4andtitle_fig2(3);
get_pos_lh4title_height_fig2 = get_pos_lh4andtitle_fig2(4) - get_pos_lh4_fig2(4); 

lh4_fig2.EdgeColor = figProps.Legend.ColorAndStyling.EdgeColor;

% Set background color of fourth legend to white (instead of grey by default)
set(lh4_fig2,'Color',get(ax1_fig2,'Color'));

% Adjust horizontal spacing for legend lh4
lh4_fig2.ItemTokenSize(1) = figProps.Legend.Fonts.ItemTokenSize;

lh4titleTop_fig2 = annotation('line',[0.5, 0.75], [0.40, 0.40], 'Color', 'black');
drawnow
set(lh4titleTop_fig2,'Units','centimeters');
drawnow
get_pos_lh4titleTop_fig2 = get(lh4titleTop_fig2,'Position');

lh4titleBottom_fig2 = annotation('line',[0.5, 0.75], [0.35, 0.35], 'Color', 'black');
drawnow
set(lh4titleBottom_fig2,'Units','centimeters');
drawnow
get_pos_lh4titleBottom_fig2 = get(lh4titleBottom_fig2,'Position');

lh4titleBottom2_fig2 = annotation('line',[0.5, 0.75], [0.30, 0.30], 'Color', 'black');
drawnow
set(lh4titleBottom2_fig2,'Units','centimeters');
drawnow
get_pos_lh4titleBottom2_fig2 = get(lh4titleBottom2_fig2,'Position');

%% Fifth legend [Patankar 2001, Spherical Data]

% Create array for markers and labels for fifth legend
lh5_fig2_arrayMarker = ["none" "<" "v" ">" "^"];
lh5_fig2_label = cellstr(["$\ H/d_\mathrm{max}$","$\ \ \ \ \ 4$","$\ \ \ \ \ 6$","$\ \ \ \ 12$","$\ \ \ \ 48$"]);

% Copy the axes and plot the fifth legend
ax5_fig2 = copyobj(ax4_fig2,gcf); %Copies ah1 object (but duplicates data points)
delete(get(ax5_fig2,'Children')) %Deletes duplicate data points

% Remove fifth axis visibility 
set(ax5_fig2, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

% Plot "dummy" data for fifth legend; tried using single variable with loop but does not work
hold on
lh5dump1_fig2 = plot(NaN,'Marker',lh5_fig2_arrayMarker(1),'DisplayName',lh5_fig2_label{1},'LineStyle',"none",'Marker',"none",'LineStyle','none','Parent',ax5_fig2);
lh5dump2_fig2 = plot(NaN,'Marker',lh5_fig2_arrayMarker(2),'DisplayName',lh5_fig2_label{2},'LineWidth',figProps.Scatter.Markers.LineWidth2/2,'MarkerEdgeColor',"#9e0142",'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',sqrt(figProps.Scatter.Markers.SizeData2),'LineStyle','none','Parent',ax5_fig2);
lh5dump3_fig2 = plot(NaN,'Marker',lh5_fig2_arrayMarker(3),'DisplayName',lh5_fig2_label{3},'LineWidth',figProps.Scatter.Markers.LineWidth2/2,'MarkerEdgeColor',"#fdae61",'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',sqrt(figProps.Scatter.Markers.SizeData2),'LineStyle','none','Parent',ax5_fig2);
lh5dump4_fig2 = plot(NaN,'Marker',lh5_fig2_arrayMarker(4),'DisplayName',lh5_fig2_label{4},'LineWidth',figProps.Scatter.Markers.LineWidth2/2,'MarkerEdgeColor',"#66c2a5",'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',sqrt(figProps.Scatter.Markers.SizeData2),'LineStyle','none','Parent',ax5_fig2);
lh5dump5_fig2 = plot(NaN,'Marker',lh5_fig2_arrayMarker(5),'DisplayName',lh5_fig2_label{5},'LineWidth',figProps.Scatter.Markers.LineWidth2/2,'MarkerEdgeColor',"#5e4fa2",'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',sqrt(figProps.Scatter.Markers.SizeData2),'LineStyle','none','Parent',ax5_fig2);
hold off

lh5_fig2 = legend(ax5_fig2,[lh5dump1_fig2 lh5dump2_fig2 lh5dump3_fig2 lh5dump4_fig2 lh5dump5_fig2],lh5_fig2_label,'FontSize',figProps.Legend.Fonts.FontSize,'Interpreter',figProps.Legend.Interpreter);
lh5_fig2.Units = 'centimeters';
lh5_fig2.Location = 'southeast';
drawnow
get_pos_lh5_fig2 = get(lh5_fig2,'Position');
lh5title_fig2 = title(lh5_fig2,{'Patankar et al. (2001),','sphere, 2-D numerical'},'FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
drawnow
get_pos_lh5andtitle_fig2 = get(lh5_fig2,'Position'); %[left bottom width height]
get_pos_lh5title_width_fig2 = get_pos_lh5andtitle_fig2(3);
get_pos_lh5title_height_fig2 = get_pos_lh5andtitle_fig2(4) - get_pos_lh5_fig2(4); 

lh5_fig2.EdgeColor = figProps.Legend.ColorAndStyling.EdgeColor;

% Set background color of fourth legend to white (instead of grey by default)
set(lh5_fig2,'Color',get(ax1_fig2,'Color'));

% Adjust horizontal spacing for legend lh4
lh5_fig2.ItemTokenSize(1) = figProps.Legend.Fonts.ItemTokenSize;

lh5titleTop_fig2 = annotation('line',[0.5, 0.75], [0.25, 0.25], 'Color', 'black');
set(lh5titleTop_fig2,'Units','centimeters');
drawnow
get_pos_lh5titleTop_fig2 = get(lh5titleTop_fig2,'Position');

lh5titleBottom_fig2 = annotation('line',[0.5, 0.75], [0.20, 0.20], 'Color', 'black');
set(lh5titleBottom_fig2,'Units','centimeters');
drawnow
get_pos_lh5titleBottom_fig2 = get(lh5titleBottom_fig2,'Position');

lh5titleBottom2_fig2 = annotation('line',[0.5, 0.75], [0.15, 0.15], 'Color', 'black');
set(lh5titleBottom2_fig2,'Units','centimeters');
drawnow
get_pos_lh5titleBottom2_fig2 = get(lh5titleBottom2_fig2,'Position');

lh5titleBottom3_fig2 = annotation('line',[0.5, 0.75], [0.10, 0.10], 'Color', 'black');
set(lh5titleBottom3_fig2,'Units','centimeters');
drawnow
get_pos_lh5titleBottom3_fig2 = get(lh5titleBottom3_fig2,'Position');

%% Evaluate proper size and positioning
drawnow
get_pos_lh1_fig2 = get(lh1_fig2,'Position');
get_pos_lh2_fig2 = get(lh2_fig2,'Position');
get_pos_lh3_fig2 = get(lh3_fig2,'Position');
get_pos_lh4_fig2 = get(lh4_fig2,'Position');

total_linewidth_lh0_fig2 = get_pos_lh1_fig2(3) + get_pos_lh2_fig2(3) + get_pos_lh3_fig2(3) + leg1_xSpacing_centimeters + get_pos_lh4_fig2(3);

%% Reposition Zero legend [Title: Flat-plate, experimental]
drawnow
set(an_fig2,'Position',[left_pos_centimeter_ax1_fig2 + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - leg1_axis_offset, total_linewidth_lh0_fig2, get_pos_an_fig2(4)]);

drawnow
set(anTop_fig2,'Position',[left_pos_centimeter_ax1_fig2 + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - leg1_axis_offset, total_linewidth_lh0_fig2, 0]);

drawnow
set(anBottom_fig2,'Position',[left_pos_centimeter_ax1_fig2 + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - leg1_axis_offset, total_linewidth_lh0_fig2, 0]);

%% Reposition First legend [HW/(D_{max})^2, set progressive color order]
drawnow
set(lh1_fig2,'Position',[left_pos_centimeter_ax1_fig2 + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - get_pos_lh1_fig2(4) - leg1_axis_offset, get_pos_lh1_fig2(3), get_pos_lh1_fig2(4)]);
set(lh1titleTop_fig2,'Position',[left_pos_centimeter_ax1_fig2 + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - leg1_axis_offset, get_pos_lh1title_width_fig2, 0]);
set(lh1titleBottom_fig2,'Position',[left_pos_centimeter_ax1_fig2 + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - get_pos_lh1_fig2(4) - leg1_axis_offset, get_pos_lh1title_width_fig2, 0]);
set(lh1titleBottom2_fig2,'Position',[left_pos_centimeter_ax1_fig2 + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - get_pos_lh1title_height_fig2 - leg1_axis_offset, get_pos_lh1title_width_fig2, 0]);

%% Reposition Second legend [H/D, dependent color order]
drawnow
set(lh2_fig2,'Position',[left_pos_centimeter_ax1_fig2 + get_pos_lh1_fig2(3) + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - get_pos_lh2_fig2(4) - leg1_axis_offset, get_pos_lh2_fig2(3), get_pos_lh2_fig2(4)]);
set(lh2titleTop_fig2,'Position',[left_pos_centimeter_ax1_fig2 + get_pos_lh1title_width_fig2 + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - leg1_axis_offset, get_pos_lh2_fig2(3), 0]);
set(lh2titleBottom_fig2,'Position',[left_pos_centimeter_ax1_fig2 + get_pos_lh1title_width_fig2 + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - get_pos_lh2_fig2(4) - leg1_axis_offset, get_pos_lh2_fig2(3), 0]);
set(lh2titleBottom2_fig2,'Position',[left_pos_centimeter_ax1_fig2 + get_pos_lh1title_width_fig2 + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - get_pos_lh2title_height_fig2 - leg1_axis_offset, get_pos_lh2_fig2(3), 0]);

%% Reposition Third legend [AR]
drawnow
set(lh3_fig2,'Position',[left_pos_centimeter_ax1_fig2 + get_pos_lh1_fig2(3) + get_pos_lh2_fig2(3) + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - get_pos_lh3_fig2(4) - leg1_axis_offset, get_pos_lh3_fig2(3), get_pos_lh3_fig2(4)]);
set(lh3titleTop_fig2,'Position',[left_pos_centimeter_ax1_fig2 + get_pos_lh1title_width_fig2 + get_pos_lh2_fig2(3) + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - leg1_axis_offset, get_pos_lh3_fig2(3), 0]);
set(lh3titleBottom_fig2,'Position',[left_pos_centimeter_ax1_fig2 + get_pos_lh1title_width_fig2 + get_pos_lh2_fig2(3) + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - get_pos_lh3_fig2(4) - leg1_axis_offset, get_pos_lh3_fig2(3), 0]);
set(lh3titleBottom2_fig2,'Position',[left_pos_centimeter_ax1_fig2 + get_pos_lh1title_width_fig2 + get_pos_lh2_fig2(3) + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - get_pos_lh3title_height_fig2 - leg1_axis_offset, get_pos_lh3_fig2(3), 0]);

%% Reposition Fourth legend [Solvent]
set(lh4_fig2,'Position',[left_pos_centimeter_ax1_fig2 + get_pos_lh1_fig2(3) + get_pos_lh2_fig2(3) + get_pos_lh3_fig2(3) + leg1_xSpacing_centimeters + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - get_pos_lh4_fig2(4) - leg1_axis_offset, get_pos_lh4_fig2(3), get_pos_lh4_fig2(4)]);
set(lh4titleTop_fig2,'Position',[left_pos_centimeter_ax1_fig2 + get_pos_lh1title_width_fig2 + get_pos_lh2_fig2(3) + get_pos_lh3_fig2(3) + leg1_xSpacing_centimeters + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - leg1_axis_offset, get_pos_lh4_fig2(3), 0]);
set(lh4titleBottom_fig2,'Position',[left_pos_centimeter_ax1_fig2 + get_pos_lh1title_width_fig2 + get_pos_lh2_fig2(3) + get_pos_lh3_fig2(3) + leg1_xSpacing_centimeters + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - get_pos_lh4_fig2(4) - leg1_axis_offset, get_pos_lh4_fig2(3), 0]);
set(lh4titleBottom2_fig2,'Position',[left_pos_centimeter_ax1_fig2 + get_pos_lh1title_width_fig2 + get_pos_lh2_fig2(3) + get_pos_lh3_fig2(3) + leg1_xSpacing_centimeters + leg1_axis_offset, top_pos_centimeter_ax1_fig2 - get_pos_an_fig2(4) - get_pos_lh4title_height_fig2 - leg1_axis_offset, get_pos_lh4_fig2(3), 0]);

%% Reposition Fifth legend [Patankar 2001, Spherical Data]
set(lh5_fig2,'Position',[right_pos_centimeter_ax1_fig2 - get_pos_lh5title_width_fig2 - leg1_axis_offset, bottom_pos_centimeter_ax1_fig2 + leg1_axis_offset, get_pos_lh5title_width_fig2, get_pos_lh5andtitle_fig2(4)]);
set(lh5titleTop_fig2,'Position',[right_pos_centimeter_ax1_fig2 - get_pos_lh5title_width_fig2 - leg1_axis_offset, bottom_pos_centimeter_ax1_fig2 + get_pos_lh5andtitle_fig2(4) + leg1_axis_offset, get_pos_lh5title_width_fig2, 0]);
set(lh5titleBottom_fig2,'Position',[right_pos_centimeter_ax1_fig2 - get_pos_lh5title_width_fig2 - leg1_axis_offset, bottom_pos_centimeter_ax1_fig2 + leg1_axis_offset, get_pos_lh5title_width_fig2, 0]);
set(lh5titleBottom2_fig2,'Position',[right_pos_centimeter_ax1_fig2 - get_pos_lh5title_width_fig2 - leg1_axis_offset, bottom_pos_centimeter_ax1_fig2 + get_pos_lh5andtitle_fig2(4) - get_pos_lh5title_height_fig2 + leg1_axis_offset, get_pos_lh5title_width_fig2, 0]);
set(lh5titleBottom3_fig2,'Position',[right_pos_centimeter_ax1_fig2 - (4/5)*get_pos_lh5title_width_fig2 - leg1_axis_offset, bottom_pos_centimeter_ax1_fig2 + get_pos_lh5andtitle_fig2(4) - get_pos_lh5title_height_fig2 - get_pos_lh5_fig2(4)/5 + leg1_axis_offset, (3/5)*get_pos_lh5title_width_fig2, 0]);

%% Figure 4) 3-D Bar graph of critical lift-off data

%---[Figure and Axis Handles]---
fig4 = figure(4);
ax1_fig4 = gca;

% Custom axis and legend font size
AxisFontSizeMultiplier_fig4 = (8.5716/7); %1
AxisLabelFontSizeMultiplier_fig4 = (6/4.8999)*(9/8.5716); %1
LegendLabelFontSizeMultiplier_fig4 = 1; %1
LegendFontSizeMultiplier_fig4 = (8.5716/7); %1

% Custom figure size scale factors (to account for limitations in exported fig size using exportgraphics MATLAB R2024b)
expGraphicsScaleFactor_fig4 = (9/6.7028).*(9/9.1017);

%---[Set Figure Properties]---
figProps = getDefaultFigProperties('1Column',AxisFontSizeMultiplier_fig4,AxisLabelFontSizeMultiplier_fig4,LegendLabelFontSizeMultiplier_fig4,LegendFontSizeMultiplier_fig4);

%---[Set Figure Size4
set(fig4,'Units','centimeters','Position',[3 3 figProps.Figure.Position.Width.*expGraphicsScaleFactor_fig4 figProps.Figure.Position.Width*figProps.Figure.Position.HWRatio.*expGraphicsScaleFactor_fig4])
pos = get(fig4,'Position');
set(fig4,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%---[Prepare data]---
b1x = paramStruct.xParameter.ValueSorted;
[b1xSorted,b1xSortedIdx] = sortrows(b1x);
b1xSortedRoundedCellStr1 = B2KFormatValuesToStrOneSigFigMaxUnc(b1xSorted(1:10),paramStruct.xParameter.UncertaintySorted(b1xSortedIdx(1:10)));
b1xSortedRoundedCellStr2 = B2KFormatValuesToStrOneSigFigMaxUnc(b1xSorted(11:20),paramStruct.xParameter.UncertaintySorted(b1xSortedIdx(11:20)));
b1xSortedRoundedCellStr3 = B2KFormatValuesToStrOneSigFigMaxUnc(b1xSorted(21:30),paramStruct.xParameter.UncertaintySorted(b1xSortedIdx(21:30)));
b1xSortedRoundedCellStr = vertcat(b1xSortedRoundedCellStr1,b1xSortedRoundedCellStr2,b1xSortedRoundedCellStr3);

B1x = reshape(b1xSortedRoundedCellStr,paramStruct.plotData.reshapeSize);

b1y = paramStruct.zParameterHW_D_maxSq.ValueSorted;
b1ySorted = NaN(size(b1y,1),size(b1y,2));
for b1yidx = 1:length(b1y)
    b1ySorted(b1yidx) = b1y(b1xSortedIdx(b1yidx));
end
b1ySortedRounded = round(b1ySorted,1,"decimals");
b1ySortedRoundedCellStr = B2KFormatValuesToStrOneSigFigMaxUnc(paramStruct.zParameterHW_D_maxSq.ValueSorted,paramStruct.zParameterHW_D_maxSq.UncertaintySorted);
B1y = reshape(b1ySortedRoundedCellStr,paramStruct.plotData.reshapeSize);

b1z = paramStruct.yParameter.ValueSorted;
b1zSorted = NaN(size(b1z,1),size(b1z,2));
for b1zidx = 1:length(b1z)
    b1zSorted(b1zidx) = b1z(b1xSortedIdx(b1zidx));
end
B1z = reshape(b1zSorted,paramStruct.plotData.reshapeSize);

bar3h = bar3(B1z);
bar3h(1).FaceColor = '#4CC1A1';
bar3h(2).FaceColor = '#E69F00';
bar3h(3).FaceColor = '#0072B2';

% Set the z-axis limits and major ticks
ax1_fig4.ZLim = [0 1200];
ax1_fig4.ZTick = 0:200:1200;

% Enable minor ticks and specify the minor tick values
ax1_fig4.ZAxis.MinorTick = 'on';
ax1_fig4.ZAxis.MinorTickValues = 100:200:1100;

% Enable minor grid lines for the z-axis
ax1_fig4.ZMinorGrid = 'on';  % Turn on minor grid lines

ax1_fig4.GridLineStyle = '-';
ax1_fig4.MinorGridLineStyle = ':';

ax1_fig4.GridLineWidth = 0.75;
ax1_fig4.MinorGridLineWidth = 0.5;

view(-70,25)

%% Plot error bars 

b1zError = paramStruct.yParameter.UncertaintySorted;
b1zErrorSorted = NaN(size(b1zError,1),size(b1zError,2));
for b1zErroridx = 1:length(b1zError)
    b1zErrorSorted(b1zErroridx) = b1zError(b1xSortedIdx(b1zErroridx));
end
B1zError = reshape(b1zErrorSorted,paramStruct.plotData.reshapeSize);

errorCapWidth_fig4 = figProps.Line.ErrorLine.CapWidth;

[bar3ErrorRows,bar3ErrorCols] = size(B1zError);

hold on;
% Loop through each bar3 and add corresponding error bar with cap
for bar3ErrorRowsIdx = 1:bar3ErrorRows
    for bar3ErrorColsIdx = 1:bar3ErrorCols
        bar3Errorx = [bar3ErrorColsIdx, bar3ErrorColsIdx];
        bar3Errory = [bar3ErrorRowsIdx, bar3ErrorRowsIdx];
        bar3Errorz = [B1z(bar3ErrorRowsIdx,bar3ErrorColsIdx),B1z(bar3ErrorRowsIdx,bar3ErrorColsIdx) + B1zError(bar3ErrorRowsIdx,bar3ErrorColsIdx)];

        % Plot the vertical error bar
        line(bar3Errorx,bar3Errory,bar3Errorz,'color',figProps.Line.ErrorLine.Color,'Linewidth',figProps.Line.ErrorLine.LineWidth);

        % Add horizontal caps at the top of the error bar
        % X-direction cap
        line([bar3ErrorColsIdx - errorCapWidth_fig4, bar3ErrorColsIdx + errorCapWidth_fig4],[bar3ErrorRowsIdx, bar3ErrorRowsIdx],[bar3Errorz(2), bar3Errorz(2)],'color',figProps.Line.ErrorLine.Color,'LineWidth',figProps.Line.ErrorLine.LineWidth);

        % Y-direction cap
        line([bar3ErrorColsIdx, bar3ErrorColsIdx],[bar3ErrorRowsIdx - errorCapWidth_fig4, bar3ErrorRowsIdx + errorCapWidth_fig4],[bar3Errorz(2), bar3Errorz(2)],'color',figProps.Line.ErrorLine.Color,'LineWidth',figProps.Line.ErrorLine.LineWidth);
    end
end
hold off;

%%

%---[Axis properties]---
%----[Font]----
ax1_fig4.FontSize = figProps.Axis.Font.FontSize;
%----[Ticks]----
set(ax1_fig4,'xticklabel',{B1x{1,:}},'yticklabel',{B1y{:,1}})
ax1_fig4.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;
%----[Rulers]----
%----[Grids]----
%----[Labels]----
xlabel({'Archimedes';'number,';'$\mathrm{Ar}_{D_\mathrm{max}}$'},'Interpreter',figProps.Axis.Labels.LabelInterpreter)
ylabel({'Area domain size,';'$HW/D_\mathrm{max}^2$'},'Interpreter',figProps.Axis.Labels.LabelInterpreter)

zlabel({'Critical shear';'Reynolds number, $\mathrm{\check{Re}}_{\mathrm{s},D_\mathrm{max}}$'},'Interpreter',figProps.Axis.Labels.LabelInterpreter)
ax1_fig4.XLabel.FontSize = figProps.Axis.Labels.LabelFontSize;
ax1_fig4.YLabel.FontSize = figProps.Axis.Labels.LabelFontSize;
ax1_fig4.ZLabel.FontSize = figProps.Axis.Labels.LabelFontSize;
ax1_fig4.XLabel.FontWeight = figProps.Axis.Labels.LabelFontWeight;
ax1_fig4.YLabel.FontWeight = figProps.Axis.Labels.LabelFontWeight;
ax1_fig4.ZLabel.FontWeight = figProps.Axis.Labels.LabelFontWeight;

%---[Unofficial properties]---
ytickangle(45)
xtickangle(-45)

%% Figure 5) Stacked bar graph of critical lift-off data
fig5 = figure(5);
ax5 = gca;

% Custom axis and legend font size
AxisFontSizeMultiplier_fig5 = (8.5716/7);
AxisLabelFontSizeMultiplier_fig5 = (6/4.8999)*(9/8.5716); %1
LegendLabelFontSizeMultiplier_fig5 = (6/4.8999);
LegendFontSizeMultiplier_fig5 = (8.5716/7);

% Custom figure size scale factors (to account for limitations in exported fig size using exportgraphics MATLAB R2024b)
expGraphicsScaleFactor_fig5 = (9/7.9375).*(9/9.0311);

%---[Set Figure Properties]---
figProps = getDefaultFigProperties('1Column',AxisFontSizeMultiplier_fig5,AxisLabelFontSizeMultiplier_fig5,LegendLabelFontSizeMultiplier_fig5,LegendFontSizeMultiplier_fig5);

%---[Set Figure Size5
set(fig5,'Units','centimeters','Position',[3 3 figProps.Figure.Position.Width.*expGraphicsScaleFactor_fig5 figProps.Figure.Position.Width.*figProps.Figure.Position.HWRatio.*expGraphicsScaleFactor_fig5])
pos = get(fig5,'Position');
set(fig5,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%---[Figure 3) Settings]---
ylimVal_fig5 = [0 1200];

w3_25020 = 0.5;
w2_10071 = 0.35;
w1_1695 = 0.2;

barColor3_25020 = '#0072B2';
barColor2_10071 = '#E69F00';
barColor1_1695 = '#4CC1A1';
channelARTextBoxMargin = 1;
%---
numTiles_fig5 = 5;
mainTile_fig5 = tiledlayout(numTiles_fig5,1); %Outer layout
tiledlayout_fig5 = tiledlayout(mainTile_fig5,1,10); %Inner layout
xlabel(tiledlayout_fig5,'Area domain size, $HW/D_\mathrm{max}^2$','Interpreter',figProps.Axis.Labels.LabelInterpreter,'FontSize',figProps.Axis.Labels.LabelFontSize)
ylabel(tiledlayout_fig5,{'Critical shear';'Reynolds number, $\mathrm{\check{Re}}_{\mathrm{s},D_\mathrm{max}}$'},'Interpreter',figProps.Axis.Labels.LabelInterpreter,'FontSize',figProps.Axis.Labels.LabelFontSize)
tiledlayout_fig5.Layout.Tile = 1;
tiledlayout_fig5.Layout.TileSpan = [numTiles_fig5-1 1];
tiledlayout_fig5.TileSpacing = "none";
tiledlayout_fig5.Padding = "tight";

%% Organize data for bar chart 

% Number of sets
numSets = length(paramStruct.aParameter.ValueSorted) / 10;

% Initialize cell arrays for groups
aGreaterThan1 = cell(numSets, 1);
aEqualTo1 = cell(numSets, 1);
aLessThan1 = cell(numSets, 1);

xGreaterThan1 = cell(numSets, 1);
xEqualTo1 = cell(numSets, 1);
xLessThan1 = cell(numSets, 1);

yGreaterThan1 = cell(numSets, 1);
yEqualTo1 = cell(numSets, 1);
yLessThan1 = cell(numSets, 1);

yGreaterThan1Error = cell(numSets, 1);
yEqualTo1Error = cell(numSets, 1);
yLessThan1Error = cell(numSets, 1);

zGreaterThan1 = cell(numSets, 1);
zEqualTo1 = cell(numSets, 1);
zLessThan1 = cell(numSets, 1);

zGreaterThan1Rounded = cell(numSets, 1);
zEqualTo1Rounded = cell(numSets, 1);
zLessThan1Rounded = cell(numSets, 1);

% Iterate through each set of 10 values
for setIdx = 1:numSets
    % Determine the range of indices for the current set
    startIdx = (setIdx - 1) * 10 + 1;
    endIdx = setIdx * 10;

    % Get the current set of values
    aValues = paramStruct.aParameter.ValueSorted(startIdx:endIdx);
    xValues = paramStruct.xParameter.ValueSorted(startIdx:endIdx);
    yValues = paramStruct.yParameter.ValueSorted(startIdx:endIdx);
    zValues = paramStruct.zParameterHW_D_maxSq.ValueSorted(startIdx:endIdx);

    yValuesError = paramStruct.yParameter.UncertaintySorted(startIdx:endIdx);

    % Initialize temporary arrays for the current set
    temp_aGreaterThan1 = [];
    temp_aEqualTo1 = [];
    temp_aLessThan1 = [];

    temp_xGreaterThan1 = [];
    temp_xEqualTo1 = [];
    temp_xLessThan1 = [];

    temp_yGreaterThan1 = [];
    temp_yEqualTo1 = [];
    temp_yLessThan1 = [];

    temp_zGreaterThan1 = [];
    temp_zEqualTo1 = [];
    temp_zLessThan1 = [];

    % Initialize temporary arrays for the current set for error-uncertainty
    temp_yGreaterThan1_Error = [];
    temp_yEqualTo1_Error = [];
    temp_yLessThan1_Error = [];

    % Classify the values based on aParameter.ValueSorted
    for i = 1:10
        if aValues(i) > 1
            temp_aGreaterThan1 = [temp_aGreaterThan1; aValues(i)];
            temp_xGreaterThan1 = [temp_xGreaterThan1; xValues(i)];
            temp_yGreaterThan1 = [temp_yGreaterThan1; yValues(i)];
            temp_zGreaterThan1 = [temp_zGreaterThan1; zValues(i)];

            temp_yGreaterThan1_Error = [temp_yGreaterThan1_Error; yValuesError(i)];
        elseif aValues(i) == 1
            temp_aEqualTo1 = [temp_aEqualTo1; aValues(i)];
            temp_xEqualTo1 = [temp_xEqualTo1; xValues(i)];
            temp_yEqualTo1 = [temp_yEqualTo1; yValues(i)];
            temp_zEqualTo1 = [temp_zEqualTo1; zValues(i)];

            temp_yEqualTo1_Error = [temp_yEqualTo1_Error; yValuesError(i)];
        else
            temp_aLessThan1 = [temp_aLessThan1; aValues(i)];
            temp_xLessThan1 = [temp_xLessThan1; xValues(i)];
            temp_yLessThan1 = [temp_yLessThan1; yValues(i)];
            temp_zLessThan1 = [temp_zLessThan1; zValues(i)];

            temp_yLessThan1_Error = [temp_yLessThan1_Error; yValuesError(i)];
        end
    end

    % Sort the groups based on zParameter
    [~, sortIdx] = sort(temp_zGreaterThan1);
    temp_aGreaterThan1 = temp_aGreaterThan1(sortIdx);
    temp_xGreaterThan1 = temp_xGreaterThan1(sortIdx);
    temp_yGreaterThan1 = temp_yGreaterThan1(sortIdx);
    temp_zGreaterThan1 = temp_zGreaterThan1(sortIdx);

    temp_yGreaterThan1_Error = temp_yGreaterThan1_Error(sortIdx);

    [~, sortIdx] = sort(temp_zEqualTo1);
    temp_aEqualTo1 = temp_aEqualTo1(sortIdx);
    temp_xEqualTo1 = temp_xEqualTo1(sortIdx);
    temp_yEqualTo1 = temp_yEqualTo1(sortIdx);
    temp_zEqualTo1 = temp_zEqualTo1(sortIdx);

    temp_yEqualTo1_Error = temp_yEqualTo1_Error(sortIdx);

    [~, sortIdx] = sort(temp_zLessThan1);
    temp_aLessThan1 = temp_aLessThan1(sortIdx);
    temp_xLessThan1 = temp_xLessThan1(sortIdx);
    temp_yLessThan1 = temp_yLessThan1(sortIdx);
    temp_zLessThan1 = temp_zLessThan1(sortIdx);

    temp_yLessThan1_Error = temp_yLessThan1_Error(sortIdx);

    % Store the sorted arrays in the corresponding cells
    aGreaterThan1{setIdx} = temp_aGreaterThan1;
    aEqualTo1{setIdx} = temp_aEqualTo1;
    aLessThan1{setIdx} = temp_aLessThan1;

    xGreaterThan1{setIdx} = temp_xGreaterThan1;
    xEqualTo1{setIdx} = temp_xEqualTo1;
    xLessThan1{setIdx} = temp_xLessThan1;

    yGreaterThan1{setIdx} = temp_yGreaterThan1;
    yEqualTo1{setIdx} = temp_yEqualTo1;
    yLessThan1{setIdx} = temp_yLessThan1;

    yGreaterThan1Error{setIdx} = temp_yGreaterThan1_Error;
    yEqualTo1Error{setIdx} = temp_yEqualTo1_Error;
    yLessThan1Error{setIdx} = temp_yLessThan1_Error;

    zGreaterThan1{setIdx} = temp_zGreaterThan1;
    zEqualTo1{setIdx} = temp_zEqualTo1;
    zLessThan1{setIdx} = temp_zLessThan1;

    zGreaterThan1Rounded{setIdx} = round(zGreaterThan1{setIdx},2,"significant");
    zEqualTo1Rounded{setIdx} = round(zEqualTo1{setIdx},2,"significant");
    zLessThan1Rounded{setIdx} = round(zLessThan1{setIdx},2,"significant");
end

% Create axis tick labels for zGreaterThan1
zGreaterThan1_labels = {'2.1';'3.3';'4.2';'5.3'};

% Create axis tick labels for zEqualTo1
zEqualTo1_labels = {'2.5';'3.1';'4.0'};

% Create axis tick labels for zLessThan1
zLessThan1_labels = {'1.6';'2.1';'2.7'};

%% H>W

%
tile1_HgreaterW_Fig5 = nexttile(tiledlayout_fig5,[1 4]);
ax1_Fig5 = gca;
ax1_Fig5.XTick = [1 2 3 4];
x_tile1_loc = [1 2 3 4];
ylim(ylimVal_fig5)
ax1_Fig5.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;
ax1_Fig5.XTickLabels = zGreaterThan1_labels;
ax1_Fig5.XTickLabelRotation = 45;

set(ax1_Fig5,'YMinorTick','on');
set(ax1_Fig5,'YMinorGrid','on');
ax1_fig5_yAx = get(ax1_Fig5,'YAxis');

ax1_fig5_yAx.MinorTickValues = 100:200:1100;

ax1_fig5.GridLineStyle = '-';
ax1_fig5.MinorGridLineStyle = ':';

ax1_fig5.GridLineWidth = 0.75;
ax1_fig5.MinorGridLineWidth = 0.5;

text(2.5,1100,'$H > W$','VerticalAlignment','middle', 'HorizontalAlignment','center','Interpreter','latex','FontSize',figProps.Legend.Fonts.FontSize,'BackgroundColor','w','EdgeColor','k','Margin',channelARTextBoxMargin)
grid on
ax1_Fig5.FontSize = figProps.Axis.Font.FontSize;
box off

%
x_tile1_25020 = zGreaterThan1{3};
y_tile1_25020 = yGreaterThan1{3};
y_tile1_25020_Error = yGreaterThan1Error{3};

hold on
b3_tile1 = bar(x_tile1_loc,y_tile1_25020,w3_25020,'FaceColor',barColor3_25020);
b3_tile1_error = errorbar(x_tile1_loc,y_tile1_25020,y_tile1_25020_Error,y_tile1_25020_Error);
b3_tile1_error.Color = figProps.Line.ErrorLine.Color;
b3_tile1_error.LineStyle = 'none';
b3_tile1_error.LineWidth = figProps.Line.ErrorLine.LineWidth;
b3_tile1_error.CapSize = figProps.Line.ErrorBar.CapSize;
hold off

%
x_tile1_10071 = zGreaterThan1{1};
y_tile1_10071 = yGreaterThan1{1};
y_tile1_10071_Error = yGreaterThan1Error{1};

hold on
b2_tile1 = bar(x_tile1_loc,y_tile1_10071,w2_10071,'FaceColor',barColor2_10071);
b2_tile1_error = errorbar(x_tile1_loc,y_tile1_10071,y_tile1_10071_Error,y_tile1_10071_Error);
b2_tile1_error.Color = figProps.Line.ErrorLine.Color;
b2_tile1_error.LineStyle = 'none';
b2_tile1_error.LineWidth = figProps.Line.ErrorLine.LineWidth;
b2_tile1_error.CapSize = figProps.Line.ErrorBar.CapSize .* (w2_10071/w3_25020);
hold off

%
x_tile1_1695 = zGreaterThan1{2};
y_tile1_1695 = yGreaterThan1{2};
y_tile1_1695_Error = yGreaterThan1Error{2};

hold on
b1_tile1 = bar(x_tile1_loc,y_tile1_1695,w1_1695,'FaceColor',barColor1_1695);
b1_tile1_error = errorbar(x_tile1_loc,y_tile1_1695,y_tile1_1695_Error,y_tile1_1695_Error);
b1_tile1_error.Color = figProps.Line.ErrorLine.Color;
b1_tile1_error.LineStyle = 'none';
b1_tile1_error.LineWidth = figProps.Line.ErrorLine.LineWidth;
b1_tile1_error.CapSize = figProps.Line.ErrorBar.CapSize .* (w1_1695/w2_10071);
hold off

%% H=W

%
tile2_HgreaterW_Fig5 = nexttile(tiledlayout_fig5,[1 3]);
ax2_Fig5 = gca;
ax2_Fig5.XTick = [1 2 3];
x_tile2_loc = [1 2 3];
ylim(ylimVal_fig5)
ax2_Fig5.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;
ax2_Fig5.YTickLabel = [];
ax2_Fig5.YAxis.TickLength = [0,0];
ax2_Fig5.XTickLabels = zEqualTo1_labels;
ax2_Fig5.XTickLabelRotation = 45;

set(ax2_Fig5,'YMinorTick','on');
set(ax2_Fig5,'YMinorGrid','on');
ax2_fig5_yAx = get(ax2_Fig5,'YAxis');

ax2_fig5_yAx.MinorTickValues = 100:200:1100;

ax2_fig5.GridLineStyle = '-';
ax2_fig5.MinorGridLineStyle = ':';

ax2_fig5.GridLineWidth = 0.75;
ax2_fig5.MinorGridLineWidth = 0.5;

text(2,1100,'$H = W$','VerticalAlignment','middle', 'HorizontalAlignment','center','Interpreter','latex','FontSize',figProps.Legend.Fonts.FontSize,'BackgroundColor','w','EdgeColor','k','Margin',channelARTextBoxMargin)
grid on
ax2_Fig5.FontSize = figProps.Axis.Font.FontSize;
box off

%
x_tile2_25020 = zEqualTo1{3};
y_tile2_25020 = yEqualTo1{3};
y_tile2_25020_Error = yEqualTo1Error{3};

hold on
b3_tile2 = bar(x_tile2_loc,y_tile2_25020,w3_25020,'FaceColor',barColor3_25020);
b3_tile2_error = errorbar(x_tile2_loc,y_tile2_25020,y_tile2_25020_Error,y_tile2_25020_Error);
b3_tile2_error.Color = figProps.Line.ErrorLine.Color;
b3_tile2_error.LineStyle = 'none';
b3_tile2_error.LineWidth = figProps.Line.ErrorLine.LineWidth;
b3_tile2_error.CapSize = figProps.Line.ErrorBar.CapSize;
hold off

%
x_tile2_10071 = zEqualTo1{1};
y_tile2_10071 = yEqualTo1{1};
y_tile2_10071_Error = yEqualTo1Error{1};

hold on
b2_tile2 = bar(x_tile2_loc,y_tile2_10071,w2_10071,'FaceColor',barColor2_10071);
b2_tile2_error = errorbar(x_tile2_loc,y_tile2_10071,y_tile2_10071_Error,y_tile2_10071_Error);
b2_tile2_error.Color = figProps.Line.ErrorLine.Color;
b2_tile2_error.LineStyle = 'none';
b2_tile2_error.LineWidth = figProps.Line.ErrorLine.LineWidth;
b2_tile2_error.CapSize = figProps.Line.ErrorBar.CapSize .* (w2_10071/w3_25020);
hold off

%
x_tile2_1695 = zEqualTo1{2};
y_tile2_1695 = yEqualTo1{2};
y_tile2_1695_Error = yEqualTo1Error{2};

hold on
b1_tile2 = bar(x_tile2_loc,y_tile2_1695,w1_1695,'FaceColor',barColor1_1695);
b1_tile2_error = errorbar(x_tile2_loc,y_tile2_1695,y_tile2_1695_Error,y_tile2_1695_Error);
b1_tile2_error.Color = figProps.Line.ErrorLine.Color;
b1_tile2_error.LineStyle = 'none';
b1_tile2_error.LineWidth = figProps.Line.ErrorLine.LineWidth;
b1_tile2_error.CapSize = figProps.Line.ErrorBar.CapSize .* (w1_1695/w2_10071);
hold off

%% H<W

%
tile3_HgreaterW_Fig5 = nexttile(tiledlayout_fig5,[1 3]);
ax3_Fig5 = gca;
ax3_Fig5.XTick = [1 2 3];
x_tile3_loc = [1 2 3];
ylim(ylimVal_fig5)
ax3_Fig5.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;
ax3_Fig5.YTickLabel = [];
ax3_Fig5.YAxis.TickLength = [0,0];
ax3_Fig5.XTickLabels = zLessThan1_labels;
ax3_Fig5.XTickLabelRotation = 45;

set(ax3_Fig5,'YMinorTick','on');
set(ax3_Fig5,'YMinorGrid','on');
ax3_fig5_yAx = get(ax3_Fig5,'YAxis');

ax3_fig5_yAx.MinorTickValues = 100:200:1100;

ax3_fig5.GridLineStyle = '-';
ax3_fig5.MinorGridLineStyle = ':';

ax3_fig5.GridLineWidth = 0.75;
ax3_fig5.MinorGridLineWidth = 0.5;

text(2,1100,'$H < W$','VerticalAlignment','middle', 'HorizontalAlignment','center','Interpreter','latex','FontSize',figProps.Legend.Fonts.FontSize,'BackgroundColor','w','EdgeColor','k','Margin',channelARTextBoxMargin)
grid on
ax3_Fig5.FontSize = figProps.Axis.Font.FontSize;
box off

%
x_tile3_25020 = zLessThan1{3};
y_tile3_25020 = yLessThan1{3};
y_tile3_25020_Error = yLessThan1Error{3};

hold on
b3_tile3 = bar(x_tile3_loc,y_tile3_25020,w3_25020,'FaceColor',barColor3_25020);
b3_tile3_error = errorbar(x_tile3_loc,y_tile3_25020,y_tile3_25020_Error,y_tile3_25020_Error);
b3_tile3_error.Color = figProps.Line.ErrorLine.Color;
b3_tile3_error.LineStyle = 'none';
b3_tile3_error.LineWidth = figProps.Line.ErrorLine.LineWidth;
b3_tile3_error.CapSize = figProps.Line.ErrorBar.CapSize;
hold off

%
x_tile3_10071 = zLessThan1{1};
y_tile3_10071 = yLessThan1{1};
y_tile3_10071_Error = yLessThan1Error{1};

hold on
b2_tile3 = bar(x_tile3_loc,y_tile3_10071,w2_10071,'FaceColor',barColor2_10071);
b2_tile3_error = errorbar(x_tile3_loc,y_tile3_10071,y_tile3_10071_Error,y_tile3_10071_Error);
b2_tile3_error.Color = figProps.Line.ErrorLine.Color;
b2_tile3_error.LineStyle = 'none';
b2_tile3_error.LineWidth = figProps.Line.ErrorLine.LineWidth;
b2_tile3_error.CapSize = figProps.Line.ErrorBar.CapSize .* (w2_10071/w3_25020);
hold off

%
x_tile3_1695 = zLessThan1{2};
y_tile3_1695 = yLessThan1{2};
y_tile3_1695_Error = yLessThan1Error{2};

hold on
b1_tile3 = bar(x_tile3_loc,y_tile3_1695,w1_1695,'FaceColor',barColor1_1695);
b1_tile3_error = errorbar(x_tile3_loc,y_tile3_1695,y_tile3_1695_Error,y_tile3_1695_Error);
b1_tile3_error.Color = figProps.Line.ErrorLine.Color;
b1_tile3_error.LineStyle = 'none';
b1_tile3_error.LineWidth = figProps.Line.ErrorLine.LineWidth;
b1_tile3_error.CapSize = figProps.Line.ErrorBar.CapSize .* (w1_1695/w2_10071);
hold off

%%
tileLeg = nexttile(mainTile_fig5,numTiles_fig5);
hold on
bar(nan,nan,w1_1695,'FaceColor',barColor1_1695);
bar(nan,nan,w2_10071,'FaceColor',barColor2_10071);
bar(nan,nan,w3_25020,'FaceColor',barColor3_25020);

axis off
box off
axFig3Tile5 = gca;
axFig3Tile5.Visible = 'off';
set(axFig3Tile5,'xtick',[])
set(axFig3Tile5,'XColor', 'none','YColor','none')

leg_Fig5 = legend({'1700{   }','10070{   }','25000{   }'},'Orientation','horizontal','Interpreter','latex','Location','north','FontSize',figProps.Legend.Fonts.FontSize);
title(leg_Fig5,paramStruct.xParameter.NewStr,'Interpreter','latex','FontSize',figProps.Legend.Labels.Title.FontSize)
set(leg_Fig5,'Position',[0.2800   0.0025    0.5000    0.1200])

%% Plot settings

colorOrderRegular = [1 2 3 4 5 6 7 8 9 10];
colorOrderRatioHWandDmaxSq = [4 1 2 5 6 3 8 10 9 7];
colorOrderRatioHandDmax = [4 1 6 5 2 3 10 8 9 7];

arrayColorCOMSOLSorted = [4 1 6 5 2 3 10 8 9 7];

%% [ALL FORCE EVAL SETUP] Force evaluation plots [ALL LOOPS]
% % channel_01_11
%channel(10)
%-surface(6)
%--botClear(3)
%---xDistInletObj(5)
%----solvent(3)
%-----totalDrag
%-----pressureDrag
%-----viscousDrag
%-----totalLift
%-----pressureLift
%-----viscousLift

% % channel_00
%channel(10)
%-surface(6)
%--botClear(5)
%---mesh(11)
%----solvent(3)
%-----totalDrag
%-----pressureDrag
%-----viscousDrag
%-----totalLift
%-----pressureLift
%-----viscousLift

%ALL PARAMETERS
%-eval_data_01_11 = 'xDistInletObj'; %{'totalDrag','pressureDrag','viscousDrag','totalLift','pressureLift','viscousLift'}
%-xParameterUsed = 'botClear'; %{'xDistInletObj','botClear'}
%- [index of constant value for other xParameter]

indices_HgreaterW = [2 8 9 7];
indices_HequalW = [5 3 10]; %[5 3 11];
indices_HlesserW = [4 1 6];

xTileIdxs_Cell_ChannelAR_ALL_01_11 = {indices_HgreaterW,indices_HequalW,indices_HlesserW};
yTileIdxs_Cell_Solvent_ALL_01_11 = {1 2 3};

xTileIdxs_Cell_ChannelAR_ALL_00 = {indices_HgreaterW,indices_HequalW,indices_HlesserW};
yTileIdxs_Cell_Solvent_ALL_00 = {1 2 3};

xDistInletObj_array = [1E-4 9E-4 0.0024 0.00301 0.023]; %[m]

idxXParameterHeldConstant = 1;

%% Figure 6) Drag Force vs. X-Distance between particle and inlet

%---[Figure and Axis Handles]---
fig6 = figure(6);
ax6 = gca;

%---[Set Figure Properties]---
figProps = getDefaultFigProperties('2Column',1,1,1,1);
% figProps = getDefaultFigProperties('1Column',1,1,1,1); %doesn't work; too small

figProps.Figure.Position.Width = 19; %19, 19, 20, 24, 19
figProps.Figure.Position.HWRatio = 0.70; %0.65, 0.70, 0.70, 0.70, 0.70

%---[Set Figure Size6 - 
set(fig6,'Units','centimeters','Position',[1 1 figProps.Figure.Position.Width figProps.Figure.Position.Width*figProps.Figure.Position.HWRatio])
pos = get(fig6,'Position');
set(fig6,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%---[Parameters]---
eval_data_01_11 = 'totalDrag';
% % xParameterUsed = {xDistInletObj','botClear'};
xParameterUsed = 'xDistInletObj';

genForcePlot(xParameterUsed, channel_01_11, eval_data_01_11, xTileIdxs_Cell_ChannelAR_ALL_01_11, colorOrderRatioHWandDmaxSq, paramStruct.plotData.arrayColor, colorOrderRegular, figProps);

%% Figure 7) Lift Force vs. X-Distance between particle and inlet

%---[Figure and Axis Handles]---
fig7 = figure(7);
ax7 = gca;

%---[Set Figure Properties]---
figProps = getDefaultFigProperties('2Column',1,1,1,1);

figProps.Figure.Position.Width = 19; %19, 19, 20, 24, 19
figProps.Figure.Position.HWRatio = 0.70; %0.65, 0.70, 0.70, 0.70, 0.70

%---[Set Figure Size7 - 
set(fig7,'Units','centimeters','Position',[1 1 figProps.Figure.Position.Width figProps.Figure.Position.Width*figProps.Figure.Position.HWRatio])
pos = get(fig7,'Position');
set(fig7,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%---[Parameters]---
eval_data_01_11 = 'totalLift';
xParameterUsed = 'xDistInletObj';

genForcePlot(xParameterUsed, channel_01_11, eval_data_01_11, xTileIdxs_Cell_ChannelAR_ALL_01_11, colorOrderRatioHWandDmaxSq, paramStruct.plotData.arrayColor, colorOrderRegular, figProps);

%% Figure 8) Drag Force vs. Bottom clearance between particle and inlet

%---[Figure and Axis Handles]---
fig8 = figure(8);
ax8 = gca;

%---[Set Figure Properties]---
figProps = getDefaultFigProperties('2Column',1,1,1,1);

%---[Set Figure Size8 - 
set(fig8,'Units','centimeters','Position',[3 3 figProps.Figure.Position.Width figProps.Figure.Position.Width*figProps.Figure.Position.HWRatio])
pos = get(fig8,'Position');
set(fig8,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%---[Parameters]---
eval_data_01_11 = 'totalDrag';
xParameterUsed = 'botClear';

genForcePlot(xParameterUsed, channel_01_11, eval_data_01_11, xTileIdxs_Cell_ChannelAR_ALL_01_11, colorOrderRatioHWandDmaxSq, paramStruct.plotData.arrayColor, colorOrderRegular, figProps);

%% Figure 9) Lift Force vs. Bottom clearance between particle and inlet

%---[Figure and Axis Handles]---
fig9 = figure(9);
ax9 = gca;

%---[Set Figure Properties]---
figProps = getDefaultFigProperties('2Column',1,1,1,1);

%---[Set Figure Size9 - 
set(fig9,'Units','centimeters','Position',[3 3 figProps.Figure.Position.Width figProps.Figure.Position.Width*figProps.Figure.Position.HWRatio])
pos = get(fig9,'Position');
set(fig9,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%---[Parameters]---
eval_data_01_11 = 'totalLift';
xParameterUsed = 'botClear';

genForcePlot(xParameterUsed, channel_01_11, eval_data_01_11, xTileIdxs_Cell_ChannelAR_ALL_01_11, colorOrderRatioHWandDmaxSq, paramStruct.plotData.arrayColor, colorOrderRegular, figProps);

%% Figure 71) Mesh independence plot [Lift force vs. xDistInletObj (5) vs. mesh config (9)]

meshIndepData1_04_Isopropanol_Lift_xDistObjInlet = readmatrix("Isopropanol_LiftForce_xDistObjInlet.xlsx");
meshIndepData2_04_Water_Lift_xDistObjInlet = readmatrix("Water_LiftForce_xDistObjInlet.xlsx");
meshIndepData3_04_Methanol_Lift_xDistObjInlet = readmatrix("Methanol_LiftForce_xDistObjInlet.xlsx");
meshIndepData4_04_Isopropanol_Drag_xDistObjInlet = readmatrix("Isopropanol_DragForce_xDistObjInlet.xlsx");
meshIndepData5_04_Water_Drag_xDistObjInlet = readmatrix("Water_DragForce_xDistObjInlet.xlsx");
meshIndepData6_04_Methanol_Drag_xDistObjInlet = readmatrix("Methanol_DragForce_xDistObjInlet.xlsx");

meshLegStr = {'Mesh \#1 (Extremely coarse)','Mesh \#2 (Extra coarse)','Mesh \#3 (Coarser)','Mesh \#4 (Coarse)','Mesh \#5 (Normal)','Mesh \#6 (Fine)','Mesh \#7 (Finer)','Mesh \#8 (Extra fine)','Mesh \#9 (Extremely fine)'};

%---[Figure and Axis Handles]---
fig71 = figure(71);
ax71 = gca;

%---[Set Figure Properties]---
figProps = getDefaultFigProperties('2Column',1,1,1,1);

%---[Set Figure Size71]--- 
set(fig71,'Units','centimeters','Position',[3 3 figProps.Figure.Position.Width figProps.Figure.Position.Width*figProps.Figure.Position.HWRatio])
pos = get(fig71,'Position');
set(fig71,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

% Your six datasets & titles
dataAll = { ...
    meshIndepData1_04_Isopropanol_Lift_xDistObjInlet, ...
    meshIndepData2_04_Water_Lift_xDistObjInlet, ...
    meshIndepData3_04_Methanol_Lift_xDistObjInlet, ...
    meshIndepData4_04_Isopropanol_Drag_xDistObjInlet, ...
    meshIndepData5_04_Water_Drag_xDistObjInlet, ...
    meshIndepData6_04_Methanol_Drag_xDistObjInlet  ...
};
titles = { ...
    'Isopropanol Lift','Water Lift','Methanol Lift', ...
    'Isopropanol Drag','Water Drag','Methanol Drag'  ...
};

mainTileSizeY_fig71 = 124;
mainTileSizeX_fig71 = 190;
Tile1Size_fig71    = [54 125];
GapSize_fig71      = [16 125];
Tile2Size_fig71    = [54 125];
Tile3Size_fig71    = [124 65];

xticks_fig71 = 0:10:30;
xminorticks_fig71 = 5:10:25;
yticks_fig71 = 0:0.5E-6:3E-6;

% Font & marker settings
FontSizeMultiplier_fig71 = 1.75;     % adjust all font sizes
LegFontSizeMultiplier_fig71 = 1.3;
markerSize_fig71         = 5;       % size of markers
lineWidth_fig71          = 1.1;     % width of plot lines

%% 1) Create the “big” outer grid
outer = tiledlayout(fig71, ...
    mainTileSizeY_fig71, mainTileSizeX_fig71, ...
    'TileSpacing','none','Padding','none');

%% 2) Top row plots: 1×3 layout in tile #1
t1 = tiledlayout(outer, 1, 3, ...
    'TileSpacing','compact','Padding','compact');
t1.Layout.Tile     = 1;               
t1.Layout.TileSpan = Tile1Size_fig71;  

for i = 1:3
    ax = nexttile(t1, i);
    x  = dataAll{i}(1:5,1)*1000;
    hold(ax,'on');
    for m = 2:10
        y = dataAll{i}(1:5,m);
        c = lineSpec{m-1,2};  % extract color for this series
        plot(ax, x, y, '-s', ...
             'Color',           c, ...
             'MarkerFaceColor', c, ...
             'MarkerSize',      markerSize_fig71, ...
             'MarkerEdgeColor', figProps.Scatter.Markers.MarkerEdgeColor, ...
             'LineWidth',       lineWidth_fig71, ...
             'DisplayName',     meshLegStr{m-1});
    end
    hold(ax,'off');
    grid(ax,'on');
    ylim(ax, [0 2.5e-6]);
    xlim(ax, [0 25]);
    ax.XTick = xticks_fig71;
    ax.XAxis.MinorTick = 'on';
    ax.XAxis.MinorTickValues = xminorticks_fig71;
    ax.YTick = yticks_fig71;
    ax.TickLength = [0.02 0.050];
    
    % apply font settings
    ax.FontSize             = figProps.Axis.Font.FontSize * FontSizeMultiplier_fig71;
    ax.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;

    ax.XLabel.Interpreter   = figProps.Axis.Labels.LabelInterpreter;
    ax.XLabel.FontSize      = figProps.Axis.Labels.LabelFontSize * FontSizeMultiplier_fig71;
    ax.XLabel.FontWeight    = figProps.Axis.Labels.LabelFontWeight;
    
    ax.YLabel.Interpreter   = figProps.Axis.Labels.LabelInterpreter;
    ax.YLabel.FontSize      = figProps.Axis.Labels.LabelFontSize * FontSizeMultiplier_fig71;
    ax.YLabel.FontWeight    = figProps.Axis.Labels.LabelFontWeight;
end

% shared Y-label for top row
ylabel(t1, 'Lift force (N)', ...
    'Interpreter', figProps.Axis.Labels.LabelInterpreter, ...
    'FontSize',    figProps.Axis.Labels.LabelFontSize * FontSizeMultiplier_fig71, ...
    'FontWeight',  figProps.Axis.Labels.LabelFontWeight);

%% 3) Blank gap between rows
gapLay = tiledlayout(outer, 1, 1, ...
    'TileSpacing','none','Padding','none');
gapLay.Layout.Tile     = Tile1Size_fig71(1)*mainTileSizeX_fig71 + 1;
gapLay.Layout.TileSpan = GapSize_fig71;
nexttile(gapLay);  axis off

%% 4) Bottom row plots: 1×3 layout after the gap
t2 = tiledlayout(outer, 1, 3, ...
    'TileSpacing','compact','Padding','compact');
t2.Layout.Tile     = (Tile1Size_fig71(1) + GapSize_fig71(1)) * mainTileSizeX_fig71 + 1;
t2.Layout.TileSpan = Tile2Size_fig71;

for i = 1:3
    ax = nexttile(t2, i);
    x  = dataAll{i+3}(1:5,1)*1000;
    hold(ax,'on');
    for m = 2:10
        y = dataAll{i+3}(1:5,m);
        c = lineSpec{m-1,2};
        plot(ax, x, y, '-s', ...
             'Color',           c, ...
             'MarkerFaceColor', c, ...
             'MarkerSize',      markerSize_fig71, ...
             'MarkerEdgeColor', figProps.Scatter.Markers.MarkerEdgeColor, ...
             'LineWidth',       lineWidth_fig71, ...
             'DisplayName',     meshLegStr{m-1});
    end
    hold(ax,'off');
    grid(ax,'on');
    ylim(ax, [0 2.5e-6]);
    xlim(ax, [0 25]);
    ax.XTick = xticks_fig71;
    ax.XAxis.MinorTick = 'on';
    ax.XAxis.MinorTickValues = xminorticks_fig71;
    ax.YTick = yticks_fig71;
    ax.TickLength = [0.02 0.050];
    
    % apply font settings
    ax.FontSize             = figProps.Axis.Font.FontSize * FontSizeMultiplier_fig71;
    ax.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;

    ax.XLabel.Interpreter   = figProps.Axis.Labels.LabelInterpreter;
    ax.XLabel.FontSize      = figProps.Axis.Labels.LabelFontSize * FontSizeMultiplier_fig71;
    ax.XLabel.FontWeight    = figProps.Axis.Labels.LabelFontWeight;
    
    ax.YLabel.Interpreter   = figProps.Axis.Labels.LabelInterpreter;
    ax.YLabel.FontSize      = figProps.Axis.Labels.LabelFontSize * FontSizeMultiplier_fig71;
    ax.YLabel.FontWeight    = figProps.Axis.Labels.LabelFontWeight;
end

% shared Y-label for bottom row
ylabel(t2, 'Drag force (N)', ...
    'Interpreter', figProps.Axis.Labels.LabelInterpreter, ...
    'FontSize',    figProps.Axis.Labels.LabelFontSize * FontSizeMultiplier_fig71, ...
    'FontWeight',  figProps.Axis.Labels.LabelFontWeight);

%% 5) Shared X-label at bottom
xlabel(outer, 'X-distance between particle and inlet (mm)', ...
    'Interpreter', figProps.Axis.Labels.LabelInterpreter, ...
    'FontSize',    figProps.Axis.Labels.LabelFontSize * FontSizeMultiplier_fig71, ...
    'FontWeight',  figProps.Axis.Labels.LabelFontWeight);

%% 6) Legend block: 1×1 layout at right, spanning full height
legLay = tiledlayout(outer, 1, 1, ...
    'TileSpacing','none','Padding','none');
legLay.Layout.Tile     = Tile1Size_fig71(2) + 1;
legLay.Layout.TileSpan = Tile3Size_fig71;

legAx = nexttile(legLay);
axis(legAx, 'off');
hLines = flipud(findall(t1,'Type','Line'));
hLeg = legend(legAx, hLines, meshLegStr, ...
    'Location','east', ...
    'Orientation','vertical','Interpreter',figProps.Legend.Interpreter);
hLeg.FontSize        = figProps.Legend.Fonts.FontSize * LegFontSizeMultiplier_fig71;

%% Figure 74

fig74 = figure(74);
clf(fig74);

%---[Get default figure properties and multipliers]---
figProps                    = getDefaultFigProperties('2Column',1,1,1,1);
FontSizeMultiplier_fig74    = 2.0*(12/8.3225)*(12/11.5987);
LegFontSizeMultiplier_fig74 = 1.5*(9/6.2419)*(9/8.699);

%---[Set Figure Size in cm]---
figWidth  = figProps.Figure.Position.Width;
figHeight = figWidth * figProps.Figure.Position.HWRatio;
set(fig74, ...
    'Units','centimeters', ...
    'Position',[3, 3, figWidth, figHeight], ...
    'PaperPositionMode','auto', ...
    'PaperUnits','centimeters', ...
    'PaperSize',[figWidth, figHeight]);

%---[prepare bar data]---
bar1_fig74    = ((dataAll{1,1}(1:5,7)-dataAll{1,1}(1:5,8))./dataAll{1,1}(1:5,8)).*100;
bar2_fig74    = ((dataAll{1,2}(1:5,7)-dataAll{1,2}(1:5,8))./dataAll{1,2}(1:5,8)).*100;
bar3_fig74    = ((dataAll{1,3}(1:5,7)-dataAll{1,3}(1:5,8))./dataAll{1,3}(1:5,8)).*100;
barAll1_fig74 = [bar1_fig74,bar2_fig74,bar3_fig74]';
barAll2_fig74 = [bar1_fig74;bar2_fig74;bar3_fig74];

groupIdx   = [1,2,3];
condNames1 = {'0.10','0.90','2.40','3.01','23.00'};
groupNames = {'Isopropanol','Water','Methanol'};

%---[Define bar colors]---
barColors01_fig74 = [0.30, 0.76, 0.63;...
                     0.90, 0.62, 0;...
                     0, 0.45, 0.70];
barColors02_fig74 = [0.83, 0.37, 0];

%---[normalized margins & 3:1 width ratio]---
lm    = 0.08; rm = 0.05;
bm    = 0.12; tm = 0.08;
ratio = [3,1];
tw    = 1-lm-rm;
th    = 1-bm-tm;
w1    = tw*ratio(1)/sum(ratio);
w2    = tw*ratio(2)/sum(ratio);
pos1  = [lm,      bm, w1, th];
pos2  = [lm + w1, bm, w2, th];

%---define common Y-tick values intervals---
tickVals = -15.0:5.0:15.0;
ylimVals = [-15.0,15.0];
minortickVals = -12.5:5.0:12.5;

%% Left axes (Fill Performance)
ax1 = axes(...
  'Units','normalized', ...
  'Position',pos1, ...
  'Box','off', ...
  'Layer','top' ...
);
h1 = dabarplot(barAll1_fig74, ...
    'groups',groupIdx, ...
    'xtlabels',condNames1, ...
    'errorbars',0, ...
    'colors',barColors01_fig74);
ylim(ax1,ylimVals);
yline(ax1,0,'k--','HandleVisibility','off');
ax1.XGrid = 'off';
ax1.YGrid = 'on';
ax1.YTick = tickVals;
ax1.YAxis.MinorTick = 'on';
ax1.YAxis.MinorTickValues = minortickVals;
ax1.YMinorGrid = 'on';
ylabel(ax1,'Error in lift force (\%)');

% Apply figProps to ax1
ax1.FontSize             = figProps.Axis.Font.FontSize * FontSizeMultiplier_fig74;
ax1.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;
ax1.YLabel.Interpreter   = figProps.Axis.Labels.LabelInterpreter;
ax1.YLabel.FontSize      = figProps.Axis.Labels.LabelFontSize * FontSizeMultiplier_fig74;
ax1.YLabel.FontWeight    = figProps.Axis.Labels.LabelFontWeight;

hLeg1 = legend(ax1,h1.br(1,1:3),groupNames, ...
    'Location','northwest', ...
    'Interpreter',figProps.Legend.Interpreter);
hLeg1.FontSize = figProps.Legend.Fonts.FontSize * LegFontSizeMultiplier_fig74;

% Separator line (no legend entry)
hold(ax1,'on');
xsep = ax1.XLim(2);
plot(ax1,[xsep xsep],ax1.YLim,'k-','LineWidth',1,'HandleVisibility','off');
hold(ax1,'off');

%% Right axes (Average Performance)
ax2 = axes(...
  'Units','normalized', ...
  'Position',pos2, ...
  'Box','off', ...
  'Layer','top' ...
);
h2 = dabarplot(barAll2_fig74, ...
    'errorbars',0, ...
    'xtlabels',{'Mean'}, ...
    'colors',barColors02_fig74);
ylim(ax2,ylimVals);
yline(ax2,0,'k--','HandleVisibility','off');
ax2.XGrid = 'off';
ax2.YGrid = 'on';
ax2.YTick = tickVals;
ax2.YAxis.MinorTick = 'on';
ax2.YAxis.MinorTickValues = minortickVals;
ax2.YMinorGrid = 'on';
ax2.YAxis.Color = 'none';
ax2.YAxis.TickLength = [0 0];

% Apply figProps to ax2
ax2.FontSize             = figProps.Axis.Font.FontSize * FontSizeMultiplier_fig74;
ax2.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;

% Link y-axes and set ticks
linkaxes([ax1,ax2],'y');
ax2.YTick = tickVals;

% Align tick‐label gap & tick length, hide ax2 XTickLabels
ax2.XRuler.TickLabelGapOffset = ax1.XRuler.TickLabelGapOffset;
ax2.XAxis.TickLength          = ax1.XAxis.TickLength;
ax2.XTickLabel = {};

%% Custom errorbar for ax2
mu_fig74 = mean(barAll2_fig74,'omitnan');
se_fig74 = std(barAll2_fig74,'omitnan')/sqrt(numel(barAll2_fig74));
cx_fig74 = mean(h2.br.XData);
hold(ax2,'on');
errorbar(ax2,cx_fig74,mu_fig74,se_fig74,'k','LineWidth',1,'CapSize',6,'Marker','none','HandleVisibility','off');
hold(ax2,'off');

%% Match bar widths
p_ref   = h1.br(1,1);
xr      = unique(p_ref.XData);
targetW = xr(2)-xr(1);
for p = h2.br(:)
  oldX = p.XData; xc = mean(oldX); oldW = max(oldX)-min(oldX);
  p.XData = xc + (oldX-xc)*(targetW/oldW);
end

hLeg2 = legend(ax2,h2.br,{'Mean'}, ...
    'Location','north', ...
    'Interpreter',figProps.Legend.Interpreter);
hLeg2.FontSize = figProps.Legend.Fonts.FontSize * LegFontSizeMultiplier_fig74;

%% Add xlabel before tight inset calculation
xlabel(ax1,'X-distance between particle and inlet (mm)', ...
    'Interpreter',figProps.Axis.Labels.LabelInterpreter, ...
    'FontSize',figProps.Axis.Labels.LabelFontSize*FontSizeMultiplier_fig74, ...
    'FontWeight',figProps.Axis.Labels.LabelFontWeight);

%% Auto-tight margins & reposition (accounts for xlabel)
drawnow;
ti1 = ax1.TightInset;
newPos1 = [ ...
  pos1(1)+ti1(1), ...
  pos1(2)+ti1(2), ...
  pos1(3)-ti1(1)-ti1(3), ...
  pos1(4)-ti1(2)-ti1(4)  ...
];
ax1.Position = newPos1;

newW2   = newPos1(3)*(ratio(2)/ratio(1));
newPos2 = [newPos1(1)+newPos1(3), newPos1(2), newW2, newPos1(4)];
ax2.Position = newPos2;

%% Figure 75

fig75 = figure(75);
clf(fig75);

%---[Get default figure properties and multipliers]---
figProps                    = getDefaultFigProperties('2Column',1,1,1,1);
FontSizeMultiplier_fig75    = 2.0*(12/8.3225)*(12/11.5987);
LegFontSizeMultiplier_fig75 = 1.5*(9/6.2419)*(9/8.699);

%---[Set Figure Size in cm]---
figWidth  = figProps.Figure.Position.Width;
figHeight = figWidth * figProps.Figure.Position.HWRatio;
set(fig75, ...
    'Units','centimeters', ...
    'Position',[3, 3, figWidth, figHeight], ...
    'PaperPositionMode','auto', ...
    'PaperUnits','centimeters', ...
    'PaperSize',[figWidth, figHeight]);

%---[prepare your bar data]---
bar4_fig75    = ((dataAll{1,4}(1:5,7)-dataAll{1,4}(1:5,8))./dataAll{1,4}(1:5,8)).*100;
bar5_fig75    = ((dataAll{1,5}(1:5,7)-dataAll{1,5}(1:5,8))./dataAll{1,5}(1:5,8)).*100;
bar6_fig75    = ((dataAll{1,6}(1:5,7)-dataAll{1,6}(1:5,8))./dataAll{1,6}(1:5,8)).*100;
barAll1_fig75 = [bar4_fig75,bar5_fig75,bar6_fig75]';
barAll2_fig75 = [bar4_fig75;bar5_fig75;bar6_fig75];

groupIdx   = [1,2,3];
condNames1 = {'0.10','0.90','2.40','3.01','23.00'};
groupNames = {'Isopropanol','Water','Methanol'};

%---[Define bar colors]---
barColors01_fig75 = [0.30, 0.76, 0.63;...
                     0.90, 0.62, 0;...
                     0, 0.45, 0.70];
barColors02_fig75 = [0.83, 0.37, 0];

%---[normalized margins & 3:1 width ratio]---
lm    = 0.08; rm = 0.05;
bm    = 0.12; tm = 0.08;
ratio = [3,1];
tw    = 1-lm-rm;
th    = 1-bm-tm;
w1    = tw*ratio(1)/sum(ratio);
w2    = tw*ratio(2)/sum(ratio);
pos1  = [lm,      bm, w1, th];
pos2  = [lm + w1, bm, w2, th];

%---define common Y-tick values intervals---
tickVals = -15.0:5.0:15.0;
ylimVals = [-15.0,15.0];
minortickVals = -12.5:5.0:12.5;

%% Left axes (Fill Performance)
ax1 = axes(...
  'Units','normalized', ...
  'Position',pos1, ...
  'Box','off', ...
  'Layer','top' ...
);
h1 = dabarplot(barAll1_fig75, ...
    'groups',groupIdx, ...
    'xtlabels',condNames1, ...
    'errorbars',0, ...
    'colors',barColors01_fig75);
ylim(ax1,ylimVals);
yline(ax1,0,'k--','HandleVisibility','off');
ax1.XGrid = 'off';
ax1.YGrid = 'on';
ax1.YTick = tickVals;
ax1.YAxis.MinorTick = 'on';
ax1.YAxis.MinorTickValues = minortickVals;
ax1.YMinorGrid = 'on';
ylabel(ax1,'Error in drag force (\%)');

% Apply figProps to ax1
ax1.FontSize             = figProps.Axis.Font.FontSize * FontSizeMultiplier_fig75;
ax1.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;
ax1.YLabel.Interpreter   = figProps.Axis.Labels.LabelInterpreter;
ax1.YLabel.FontSize      = figProps.Axis.Labels.LabelFontSize * FontSizeMultiplier_fig75;
ax1.YLabel.FontWeight    = figProps.Axis.Labels.LabelFontWeight;

hLeg1 = legend(ax1,h1.br(1,1:3),groupNames, ...
    'Location','northwest', ...
    'Interpreter',figProps.Legend.Interpreter);
hLeg1.FontSize = figProps.Legend.Fonts.FontSize * LegFontSizeMultiplier_fig75;

% Separator line (no legend entry)
hold(ax1,'on');
xsep = ax1.XLim(2);
plot(ax1,[xsep xsep],ax1.YLim,'k-','LineWidth',1,'HandleVisibility','off');
hold(ax1,'off');

%% Right axes (Average Performance)
ax2 = axes(...
  'Units','normalized', ...
  'Position',pos2, ...
  'Box','off', ...
  'Layer','top' ...
);
h2 = dabarplot(barAll2_fig75, ...
    'errorbars',0, ...
    'xtlabels',{'Mean'}, ...
    'colors',barColors02_fig75);
ylim(ax2,ylimVals);
yline(ax2,0,'k--','HandleVisibility','off');
ax2.XGrid = 'off';
ax2.YGrid = 'on';
ax2.YTick = tickVals;
ax2.YAxis.MinorTick = 'on';
ax2.YAxis.MinorTickValues = minortickVals;
ax2.YMinorGrid = 'on';
ax2.YAxis.Color = 'none';
ax2.YAxis.TickLength = [0 0];

% Apply figProps to ax2
ax2.FontSize             = figProps.Axis.Font.FontSize * FontSizeMultiplier_fig75;
ax2.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;

% Link y-axes and set ticks
linkaxes([ax1,ax2],'y');
ax2.YTick = tickVals;

% Align tick‐label gap & tick length, hide ax2 XTickLabels
ax2.XRuler.TickLabelGapOffset = ax1.XRuler.TickLabelGapOffset;
ax2.XAxis.TickLength          = ax1.XAxis.TickLength;
ax2.XTickLabel = {};

%% Custom errorbar for ax2
mu_fig75 = mean(barAll2_fig75,'omitnan');
se_fig75 = std(barAll2_fig75,'omitnan')/sqrt(numel(barAll2_fig75));
cx_fig75 = mean(h2.br.XData);
hold(ax2,'on');
errorbar(ax2,cx_fig74,mu_fig75,se_fig75,'k','LineWidth',1,'CapSize',6,'Marker','none','HandleVisibility','off');
hold(ax2,'off');

%% Match bar widths
p_ref   = h1.br(1,1);
xr      = unique(p_ref.XData);
targetW = xr(2)-xr(1);
for p = h2.br(:)
  oldX = p.XData; xc = mean(oldX); oldW = max(oldX)-min(oldX);
  p.XData = xc + (oldX-xc)*(targetW/oldW);
end

hLeg2 = legend(ax2,h2.br,{'Mean'}, ...
    'Location','north', ...
    'Interpreter',figProps.Legend.Interpreter);
hLeg2.FontSize = figProps.Legend.Fonts.FontSize * LegFontSizeMultiplier_fig75;

%% Add xlabel before tight inset calculation
xlabel(ax1,'X-distance between particle and inlet (mm)', ...
    'Interpreter',figProps.Axis.Labels.LabelInterpreter, ...
    'FontSize',figProps.Axis.Labels.LabelFontSize*FontSizeMultiplier_fig75, ...
    'FontWeight',figProps.Axis.Labels.LabelFontWeight);

%% Auto-tight margins & reposition (accounts for xlabel)
drawnow;
ti1 = ax1.TightInset;
newPos1 = [ ...
  pos1(1)+ti1(1), ...
  pos1(2)+ti1(2), ...
  pos1(3)-ti1(1)-ti1(3), ...
  pos1(4)-ti1(2)-ti1(4)  ...
];
ax1.Position = newPos1;

newW2   = newPos1(3)*(ratio(2)/ratio(1));
newPos2 = [newPos1(1)+newPos1(3), newPos1(2), newW2, newPos1(4)];
ax2.Position = newPos2;

%% Figure 16 - Generalized master curve 
fig16 = figure(16);
ax1_fig16 = gca;

% Custom axis and legend font size
AxisFontSizeMultiplier_fig16 = (9/6.9999); %1
AxisLabelFontSizeMultiplier_fig16 = (9/6.9999); %1
LegendLabelFontSizeMultiplier_fig16 = (6/4.8999); %1
LegendFontSizeMultiplier_fig16 = (6/4.8999); %1

% Custom figure size scale factors (to account for limitations in exported fig size using exportgraphics MATLAB R2024b)
expGraphicsScaleFactor_fig16 = 1;

% Custom tile sizing
mainTileSizeY_fig16 = 100;
maintTileSizeX_fig16 = 100;
Tile1Size_fig16 = [100 55];
Tile2Size_fig16 = [99 45];
Tile3Size_fig16 = [1 45];

%---[Get Figure Properties]---
figProps = getDefaultFigProperties('1Column',AxisFontSizeMultiplier_fig16,AxisLabelFontSizeMultiplier_fig16,LegendLabelFontSizeMultiplier_fig16,LegendFontSizeMultiplier_fig16);

%---[Set Figure Size]--- 
set(fig16,'Units','centimeters','Position',[3 3 9 7.3]) %[3 3 9 5.85]
pos = get(fig16,'Position');
set(fig16,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%% Tiledlayout setup
mainTiledLayout_fig16 = tiledlayout(mainTileSizeY_fig16,maintTileSizeX_fig16);
set(mainTiledLayout_fig16,'Units','centimeters')
mainTiledLayout_fig16.TileSpacing = "none";
mainTiledLayout_fig16.Padding = "tight";

%%
tileax1_fig16 = nexttile(mainTiledLayout_fig16,Tile1Size_fig16);
set(tileax1_fig16,'Units','centimeters')

% Plot generalized master curve for pickup from a layer of particles in liquid
% Define the equation parameters
a_Kalman = 2.7; % coefficient
fricCoeff = 0.70;

% Define the bounds for x
x_min_Kalman = 1e-3;
x_max_Kalman = 1e8;
x_Kalman = [x_min_Kalman x_max_Kalman];

% Compute the corresponding y values
y_Kalman = a_Kalman .* (x_Kalman .* fricCoeff).^(3/7);
KalmanPlotLineObj = plot(x_Kalman, y_Kalman, 'k-', 'LineWidth', 0.5);

hold on

% Plot +30% and -30% (as specified by Kalman 2005; 90% of measured points bounded within limits of +_ 30%)
yplus_Kalman = 1.3 .* a_Kalman .* (x_Kalman .* fricCoeff).^(3/7);
plot(x_Kalman, yplus_Kalman, 'k--', 'LineWidth', 0.25);

yminus_Kalman = 0.7 .* a_Kalman .* (x_Kalman .* fricCoeff).^(3/7);
plot(x_Kalman, yminus_Kalman, 'k--', 'LineWidth', 0.25);

%% Drawing original points as scatter markers
for m = 1:refSize
    pOriginal(m) = scatter(paramStruct.mParameter.ValueSorted(m),paramStruct.rParameter.ValueSorted(m),figProps.Scatter.Markers.SizeData);
    set(pOriginal(m),{'Marker','MarkerFaceColor','MarkerFaceAlpha'},{lineSpec(m,1),lineSpec(m,2),figProps.Scatter.Markers.MarkerFaceAlpha})
    hold on
end
set(pOriginal,'MarkerEdgeColor',figProps.Scatter.Markers.MarkerEdgeColor)

%%
uistack(pOriginal,'up',30)

%% Set figure properties

%----[Markers]----
set(p,'MarkerEdgeColor',figProps.Line.Markers.MarkerEdgeColor)

%---[Axis properties]---
%----[Font]----
tileax1_fig16.FontSize = figProps.Axis.Font.FontSize;
%----[Ticks]----
tileax1_fig16.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;
%----[Rulers]----
xlim([10 10000])
xticks(logspace(1,4,4));
ylim([1 200])
yticks(logspace(0,3,4));

tileax1_fig16.XAxis.Exponent = 0; %Scientific notation
xtickformat('%.0f') %Scientific notation
tileax1_fig16.YAxis.Exponent = 0; %Scientific notation
ytickformat('%.0f') %Scientific notation
set(gca,'XScale','log');
set(gca,'YScale','log');
%----[Grids]----
grid on
tileax1_fig16.GridLineStyle = figProps.Axis.Grids.GridLineStyle;
tileax1_fig16.GridLineWidth = figProps.Axis.Grids.GridLineWidth;
tileax1_fig16.GridColor = figProps.Axis.Grids.GridColor;
tileax1_fig16.MinorGridLineStyle = figProps.Axis.Grids.MinorGridLineStyle;
tileax1_fig16.MinorGridLineWidth = figProps.Axis.Grids.MinorGridLineWidth;
tileax1_fig16.MinorGridColor = figProps.Axis.Grids.MinorGridColor;
%----[Labels]----
xlabel(paramStruct.mParameter.NewStr,'Interpreter',figProps.Axis.Labels.LabelInterpreter)

tileax1_fig16.XLabel.FontSize = figProps.Axis.Labels.LabelFontSize;
tileax1_fig16.XLabel.FontWeight = figProps.Axis.Labels.LabelFontWeight;

ylabel({'Modified critical';'particle Reynolds number, $\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}^*$'},'Interpreter',figProps.Axis.Labels.LabelInterpreter)

tileax1_fig16.YLabel.FontSize = figProps.Axis.Labels.LabelFontSize;
tileax1_fig16.YLabel.FontWeight = figProps.Axis.Labels.LabelFontWeight;

%----[Box Styling]----
set(tileax1_fig16, 'Box', 'Off')

%% Add text annotation parallel to the KalmanPlotLineObj line

%---[Choose a point (x, y) on the line to place the annotation]---
x_annot_main = 700;
y_annot_main = 60;
%---[Add the text annotation]---
str_main = '$\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}^* = 2.7{{\mathrm{Ar}}_{d_\mathrm{V}}^*}^{3/7}$';
eqnTextAnnot_fig16 = B2KTextAlignedToLineLogLogPlot(tileax1_fig16,x_annot_main,y_annot_main,str_main,3/7,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter','latex');

%---[Choose a point (x, y) on the line to place the annotation]---
x_annot_plus = 30;
y_annot_plus = 13.75;
%---[Add the text annotation]---
str_plus = '$\mathrm{+}$30\%';
plusTextAnnot = B2KTextAlignedToLineLogLogPlot(tileax1_fig16,x_annot_plus,y_annot_plus,str_plus,3/7,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter','latex');

%---[Choose a point (x, y) on the line to place the annotation]---
x_annot_minus = 30;
y_annot_minus = 4;
%---[Add the text annotation]---
str_minus = '$\mathrm{-}$30\%';
minusTextAnnot = B2KTextAlignedToLineLogLogPlot(tileax1_fig16,x_annot_minus, y_annot_minus, str_minus,3/7,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter','latex');

%% Text annotation for assumed friction coefficient
%---[Choose a point (x, y) on the line to place the annotation]---
x_annot_fric = 1300;
y_annot_fric = 5;
%---[Add the text annotation]---
str_fric = sprintf('Assumed\n$\\mu_{s} = %.2f$', fricCoeff);
fricCoeffTextAnnot = text(x_annot_fric, y_annot_fric, str_fric, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', figProps.Legend.Labels.Title.FontSize, 'Interpreter', 'latex');

%% Legend parameters
%---[Specify location of legend (northwest; start from bottom up) % north, south, east, west, northeast, northwest, southeast, southwest]---
loc_leg1_fig16 = 'northwest';

%---[Define offet and spacing for legends]---
leg1_axis_offset_x_fig16 = 0.1;
leg1_axis_offset_y_fig16 = 0.0;
leg1_xSpacing_centimeters = 0.2;

%% nexttile - legend (1)
% Second legend tile (spans 7 x 4)
tileax2a_fig16 = nexttile(mainTiledLayout_fig16,Tile2Size_fig16);
set(tileax2a_fig16,'Units','centimeters')
hold on

axis off
box off
tileax2a_fig16.Visible = 'off';
set(tileax2a_fig16,'xtick',[]);
set(tileax2a_fig16,'XColor','none','YColor','none');

yShiftAnnotation_fig16b = -0.0225;

%% Get initial position
%---[Get axis position ('centimeters')]---
set(tileax2a_fig16,'Units','centimeters')
drawnow
get_ax2a_fig16_pos = get(tileax2a_fig16,'Position'); %[figtom width height]
drawnow
top_pos_centimeter_ax2a_fig16 = get_ax2a_fig16_pos(2) + get_ax2a_fig16_pos(4);
drawnow
left_pos_centimeter_ax2a_fig16 = get_ax2a_fig16_pos(1);
drawnow
right_pos_centimeter_ax2a_fig16 = get_ax2a_fig16_pos(1) + get_ax2a_fig16_pos(3);
drawnow
bottom_pos_centimeter_ax2a_fig16 = get_ax2a_fig16_pos(2);

%---[Add initial features(unpositioned; 'Units','centimeters')]---
get_pos_an_fig16 = [0 0 0 0]; %[x_begin y_begin length height]

get_pos_anTop_fig16 = [0 0 0 0]; %[x_begin y_begin dx dy]

get_pos_anBottom_fig16 = [0 0 0 0]; %[x_begin y_begin dx dy]

%% First legend [HW/{D_V}^2, set progressive color order]
%---[Legends strings - Dynamic]---
zLabelStrRounded_lh1_fig16 = B2KFormatValuesToStrOneSigFigMaxUnc(paramStruct.zParameterHW_D_VSq.ValueSorted,paramStruct.zParameterHW_D_VSq.UncertaintySorted);

%---[Create dummy plot to generate first legend since default legend for scatter cannot modify legend marker size]---
for m = 1:refSize
    pDummy(m) = plot(tileax2a_fig16,NaN,NaN,'MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none');
    set(pDummy(m),{'Marker','MarkerFaceColor'},{lineSpec(m,1),[lineSpec(m,2)]})
    hold on
end
set(pDummy,'MarkerEdgeColor',figProps.Line.Markers.MarkerEdgeColor)

%---[Plot first legend]---
hold on
lh1_fig16 = legend(tileax2a_fig16,pDummy(1+refSizeEach:refSizeEach*2),zLabelStrRounded_lh1_fig16,'FontSize',figProps.Legend.Labels.Title.FontSize,'NumColumns',1,'Interpreter',figProps.Legend.Interpreter); %circle marker
lh1_fig16.Units = 'centimeters';
lh1_fig16.Location = 'northwest';
drawnow
get_pos_lh1_fig16 = get(lh1_fig16,'Position');
lh1_fig16_title = title(lh1_fig16,paramStruct.zParameterHW_D_VSq.LegendStr,'FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
drawnow
get_pos_lh1andtitle_fig16 = get(lh1_fig16,'Position'); %[left bottom width height]
get_pos_lh1title_width_fig16 = get_pos_lh1andtitle_fig16(3);
get_pos_lh1title_height_fig16 = get_pos_lh1andtitle_fig16(4) - get_pos_lh1_fig16(4);

lh1_fig16.EdgeColor = figProps.Legend.ColorAndStyling.EdgeColor;

%---[Adjust horizontal spacing for legend lh1]---
lh1_fig16.ItemTokenSize(1) = figProps.Legend.Fonts.ItemTokenSize;

%---[Set position of first legend]---
lh1titleTop_fig16 = annotation('line',[0.5, 0.75], [0.85, 0.85], 'Color', 'black');
drawnow
set(lh1titleTop_fig16,'Units','centimeters');
drawnow
get_pos_lh1titleTop_fig16 = get(lh1titleTop_fig16,'Position'); %[x_begin y_begin dx dy]

lh1titleBottom_fig16 = annotation('line',[0.5, 0.75], [0.80, 0.80], 'Color', 'black');
drawnow
set(lh1titleBottom_fig16,'Units','centimeters');
drawnow
get_pos_lh1titleBottom_fig16 = get(lh1titleBottom_fig16,'Position'); %[x_begin y_begin dx dy]

lh1titleBottom2_fig16 = annotation('line',[0.5, 0.75], [0.75, 0.75], 'Color', 'black');
drawnow
set(lh1titleBottom2_fig16,'Units','centimeters');
drawnow
get_pos_lh1titleBottom2_fig16 = get(lh1titleBottom2_fig16,'Position'); %[x_begin y_begin dx dy]

%% Second legend [H/D, dependent color order] [TESTING]
%---[Copy the axes and plot the second legend]---
tileax2b_fig16 = copyobj(tileax2a_fig16,gcf); %Copies ah1 object (but duplicates data points)
delete(get(tileax2b_fig16,'Children')) %Deletes duplicate data points

%---[Remove second axis visibility]---
set(tileax2b_fig16, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

%---[Legends strings - dynamic]---
zLabelStrRounded_lh2_fig16 = B2KFormatValuesToStrOneSigFigMaxUnc(paramStruct.hD_VParameter.ValueSorted,paramStruct.hD_VParameter.UncertaintySorted);

%---[Plot second legend]---
hold on
lh2dump1_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2b_fig16);
lh2dump2_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2b_fig16);
lh2dump3_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2b_fig16);
lh2dump4_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2b_fig16);
lh2dump5_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2b_fig16);
lh2dump6_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2b_fig16);
lh2dump7_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2b_fig16);
lh2dump8_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2b_fig16);
lh2dump9_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2b_fig16);
lh2dump10_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2b_fig16);
hold off

lh2_fig16 = legend(tileax2b_fig16,[lh2dump1_fig16 lh2dump2_fig16 lh2dump3_fig16 lh2dump4_fig16 lh2dump5_fig16 lh2dump6_fig16 lh2dump7_fig16 lh2dump8_fig16 lh2dump9_fig16 lh2dump10_fig16],zLabelStrRounded_lh2_fig16,'FontSize',figProps.Legend.Labels.Title.FontSize,'NumColumns',1,'Interpreter',figProps.Legend.Interpreter);
lh2_fig16.Units = 'centimeters';
lh2_fig16.Location = 'north';
drawnow
get_pos_lh2_fig16 = get(lh2_fig16,'Position');
% lh2title_fig16 = title(lh2_fig16,paramStruct.hD_maxParameter.LegendStr,'FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
lh2title_fig16 = title(lh2_fig16,"$H/d_\mathrm{V}$",'FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
drawnow
get_pos_lh2andtitle_fig16 = get(lh2_fig16,'Position'); %[left bottom width height]
get_pos_lh2title_width_fig16 = get_pos_lh2andtitle_fig16(3);
get_pos_lh2title_height_fig16 = get_pos_lh2andtitle_fig16(4) - get_pos_lh2_fig16(4); 

lh2_fig16.EdgeColor = figProps.Legend.ColorAndStyling.EdgeColor;

% Set background color of second legend to white (instead of grey by default)
set(lh2_fig16,'Color',get(tileax2a_fig16,'Color'));

% Adjust horizontal spacing for legend lh2
lh2_fig16.ItemTokenSize(1) = 0; %figProps.Legend.Fonts.ItemTokenSize

% Set position relative to first legend
lh2titleTop_fig16 = annotation('line',[0.5, 0.75], [0.70, 0.70], 'Color', 'black');
drawnow
set(lh2titleTop_fig16,'Units','centimeters');
drawnow
get_pos_lh2titleTop_fig16 = get(lh2titleTop_fig16,'Position');

lh2titleBottom_fig16 = annotation('line',[0.5, 0.75], [0.65, 0.65], 'Color', 'black');
drawnow
set(lh2titleBottom_fig16,'Units','centimeters');
drawnow
get_pos_lh2titleBottom_fig16 = get(lh2titleBottom_fig16,'Position');

lh2titleBottom2_fig16 = annotation('line',[0.5, 0.75], [0.60, 0.60], 'Color', 'black');
drawnow
set(lh2titleBottom2_fig16,'Units','centimeters');
drawnow
get_pos_lh2titleBottom2_fig16 = get(lh2titleBottom2_fig16,'Position');

%% Third legend [W/D] [TESTING]

% Copy the axes and plot the third legend
tileax2c_fig16 = copyobj(tileax2b_fig16,gcf); %Copies ah1 object (but duplicates data points)
delete(get(tileax2c_fig16,'Children')) %Deletes duplicate data points

% Remove third axis visibility 
set(tileax2c_fig16, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

%---Legends strings
% zLabelStrRounded_lh3_fig16 = {'3.25','3.66','2.44','3.25','4.14','3.66','3.25','4.14','3.66','4.14'};
zLabelStrRounded_lh3_fig16 = B2KFormatValuesToStrOneSigFigMaxUnc(paramStruct.wD_VParameter.ValueSorted,paramStruct.wD_VParameter.UncertaintySorted);

% Plot third legend
hold on
lh3dump1_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2c_fig16);
lh3dump2_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2c_fig16);
lh3dump3_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2c_fig16);
lh3dump4_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2c_fig16);
lh3dump5_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2c_fig16);
lh3dump6_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2c_fig16);
lh3dump7_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2c_fig16);
lh3dump8_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2c_fig16);
lh3dump9_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2c_fig16);
lh3dump10_fig16 = plot(NaN,'Marker','none','LineStyle','none','Parent',tileax2c_fig16);
hold off

lh3_fig16 = legend(tileax2c_fig16,[lh3dump1_fig16 lh3dump2_fig16 lh3dump3_fig16 lh3dump4_fig16 lh3dump5_fig16 lh3dump6_fig16 lh3dump7_fig16 lh3dump8_fig16 lh3dump9_fig16 lh3dump10_fig16],zLabelStrRounded_lh3_fig16,'FontSize',figProps.Legend.Fonts.FontSize,'NumColumns',1,'Interpreter','latex');
lh3_fig16.Units = 'centimeters';
lh3_fig16.Location = 'northeast';
drawnow
get_pos_lh3_fig16 = get(lh3_fig16,'Position');
% lh3title_fig16 = title(lh3_fig16,"$W/d_\mathrm{max}$",'FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
lh3title_fig16 = title(lh3_fig16,"$W/d_\mathrm{V}$",'FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
drawnow
get_pos_lh3andtitle_fig16 = get(lh3_fig16,'Position'); %[left bottom width height]
get_pos_lh3title_width_fig16 = get_pos_lh3andtitle_fig16(3);
get_pos_lh3title_height_fig16 = get_pos_lh3andtitle_fig16(4) - get_pos_lh3_fig16(4); 

lh3_fig16.EdgeColor = figProps.Legend.ColorAndStyling.EdgeColor;

% Set background color of third legend to white (instead of grey by default)
set(lh3_fig16,'Color',get(tileax2a_fig16,'Color'));

% Adjust horizontal spacing for legend lh3
lh3_fig16.ItemTokenSize(1) = 0; %figProps.Legend.Fonts.ItemTokenSize

% Set position relative to second legend
lh3titleTop_fig16 = annotation('line',[0.5, 0.75], [0.55, 0.55], 'Color', 'black');
drawnow
set(lh3titleTop_fig16,'Units','centimeters');
drawnow
get_pos_lh3titleTop_fig16 = get(lh3titleTop_fig16,'Position');

lh3titleBottom_fig16 = annotation('line',[0.5, 0.75], [0.50, 0.50], 'Color', 'black');
drawnow
set(lh3titleBottom_fig16,'Units','centimeters');
drawnow
get_pos_lh3titleBottom_fig16 = get(lh3titleBottom_fig16,'Position');

lh3titleBottom2_fig16 = annotation('line',[0.5, 0.75], [0.45, 0.45], 'Color', 'black');
drawnow
set(lh3titleBottom2_fig16,'Units','centimeters');
drawnow
get_pos_lh3titleBottom2_fig16 = get(lh3titleBottom2_fig16,'Position');

%% nexttile - legend(2)
% Third legend tile (spans 3 x 4)
tileax3_fig16 = nexttile(mainTiledLayout_fig16,Tile3Size_fig16);
set(tileax3_fig16,'Units','centimeters')
hold on

axis off
box off
tileax3_fig16.Visible = 'off';
set(tileax3_fig16,'xtick',[]);
set(tileax3_fig16,'XColor','none','YColor','none');

%% Get initial position
% Get axis position ('centimeters')
set(tileax3_fig16,'Units','centimeters')
drawnow
get_ax3_fig16_pos = get(tileax3_fig16,'Position'); %[figtom width height]
drawnow
top_pos_centimeter_ax3_fig16 = get_ax3_fig16_pos(2) + get_ax3_fig16_pos(4);
drawnow
left_pos_centimeter_ax3_fig16 = get_ax3_fig16_pos(1);
drawnow
right_pos_centimeter_ax3_fig16 = get_ax3_fig16_pos(1) + get_ax3_fig16_pos(3);
drawnow
bottom_pos_centimeter_ax3_fig16 = get_ax3_fig16_pos(2);

%% Fourth legend [Solvent]
% Create array for markers and labels for fourth legend
lh4_fig16b_arrayMarkerdump = ["o" "square" "diamond"];
lh4_fig16b_typeSolvent = cellstr(["Isopropanol","Water","Methanol"]);

% Plot "dummy" data for fourth legend; tried using single variable with loop but does not work
hold on
lh4dump1_fig16 = plot(NaN,'Marker',lh4_fig16b_arrayMarkerdump(1),'DisplayName',lh4_fig16b_typeSolvent{1},'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none','Parent',tileax3_fig16);
lh4dump2_fig16 = plot(NaN,'Marker',lh4_fig16b_arrayMarkerdump(2),'DisplayName',lh4_fig16b_typeSolvent{2},'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none','Parent',tileax3_fig16);
lh4dump3_fig16 = plot(NaN,'Marker',lh4_fig16b_arrayMarkerdump(3),'DisplayName',lh4_fig16b_typeSolvent{3},'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none','Parent',tileax3_fig16);
hold off

lh4_fig16 = legend(tileax3_fig16,[lh4dump1_fig16 lh4dump2_fig16 lh4dump3_fig16],lh4_fig16b_typeSolvent,'FontSize',figProps.Legend.Fonts.FontSize,'Interpreter',figProps.Legend.Interpreter);
lh4_fig16.Units = 'centimeters';
lh4_fig16.Location = 'southwest';
drawnow
get_pos_lh4_fig16 = get(lh4_fig16,'Position');
lh4title_fig16 = title(lh4_fig16,'Solvent','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
drawnow
get_pos_lh4andtitle_fig16 = get(lh4_fig16,'Position'); %[left bottom width height]
get_pos_lh4title_width_fig16 = get_pos_lh4andtitle_fig16(3);
get_pos_lh4title_height_fig16 = get_pos_lh4andtitle_fig16(4) - get_pos_lh4_fig16(4); 

lh4_fig16.EdgeColor = figProps.Legend.ColorAndStyling.EdgeColor;

% Set background color of fourth legend to white (instead of grey by default)
set(lh4_fig16,'Color',get(tileax3_fig16,'Color'));

% Adjust horizontal spacing for legend lh4
lh4_fig16.ItemTokenSize(1) = figProps.Legend.Fonts.ItemTokenSize;

% Set position relative to third legend
lh4titleTop_fig16 = annotation('line',[0.5, 0.75], [0.40, 0.40], 'Color', 'black');
drawnow
set(lh4titleTop_fig16,'Units','centimeters');
drawnow
get_pos_lh4titleTop_fig16 = get(lh4titleTop_fig16,'Position');

lh4titleBottom_fig16 = annotation('line',[0.5, 0.75], [0.35, 0.35], 'Color', 'black');
drawnow
set(lh4titleBottom_fig16,'Units','centimeters');
drawnow
get_pos_lh4titleBottom_fig16 = get(lh4titleBottom_fig16,'Position');

lh4titleBottom2_fig16 = annotation('line',[0.5, 0.75], [0.30, 0.30], 'Color', 'black');
drawnow
set(lh4titleBottom2_fig16,'Units','centimeters');
drawnow
get_pos_lh4titleBottom2_fig16 = get(lh4titleBottom2_fig16,'Position');

%% Fifth legend
% Copy the axes and plot the fifth legend
ax5_fig16 = copyobj(tileax1_fig16,gcf); %Copies ah1 object (but duplicates data points)
delete(get(ax5_fig16,'Children')) %Deletes duplicate data points

% Remove fourth axis visibility 
set(ax5_fig16, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

% Plot "dummy" data for fourth legend; tried using single variable with loop but does not work
hold on
lh5dump1_fig16 = plot(NaN, NaN, 'k-', 'LineWidth', 0.5,'Parent',ax5_fig16);
hold off

% Legend
lh5_fig16 = legend(ax5_fig16,lh5dump1_fig16,{['Rabinovich and' newline 'Kalman (2009a), model']},'Box','off','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Legend.Interpreter);
lh5_fig16.IconColumnWidth = 10;

lh5_fig16.Units = 'centimeters';

lh5_fig16.Position = [1.9660    1    2.8613    1.1504];

%% Evaluate proper size and positioning
drawnow
get_pos_lh1_fig16 = get(lh1_fig16,'Position');
get_pos_lh2_fig16 = get(lh2_fig16,'Position');
get_pos_lh3_fig16 = get(lh3_fig16,'Position');
get_pos_lh4_fig16 = get(lh4_fig16,'Position');

total_linewidth_lh0_fig16 = get_pos_lh1_fig16(3) + get_pos_lh2_fig16(3) + get_pos_lh3_fig16(3) + leg1_xSpacing_centimeters + get_pos_lh4_fig16(3);

%% Reposition First legend [HW/(D_{V})^2, set progressive color order]
%---set legend position---[left bottom width height]
drawnow
refreshdata
set(lh1_fig16,'Position',[get_ax2a_fig16_pos(1)+leg1_axis_offset_x_fig16, get_ax2a_fig16_pos(2)+get_ax2a_fig16_pos(4)-get_pos_lh1_fig16(4)-leg1_axis_offset_y_fig16, get_pos_lh1_fig16(3), get_pos_lh1_fig16(4)]); %Issues with position calculation; need to run set command twice for proper rendering
drawnow
refreshdata
set(lh1_fig16,'Position',[get_ax2a_fig16_pos(1)+leg1_axis_offset_x_fig16, get_ax2a_fig16_pos(2)+get_ax2a_fig16_pos(4)-get_pos_lh1_fig16(4)-leg1_axis_offset_y_fig16, get_pos_lh1_fig16(3), get_pos_lh1_fig16(4)]); %Issues with position calculation; need to run set command twice for proper rendering

%---set annotation line position---[x_begin y_begin dx dy]
set(lh1titleTop_fig16,'Position',[left_pos_centimeter_ax2a_fig16+leg1_axis_offset_x_fig16, top_pos_centimeter_ax2a_fig16-get_pos_an_fig16(4)-leg1_axis_offset_y_fig16, get_pos_lh1title_width_fig16, 0]);
set(lh1titleBottom_fig16,'Position',[left_pos_centimeter_ax2a_fig16+leg1_axis_offset_x_fig16, top_pos_centimeter_ax2a_fig16-get_pos_an_fig16(4)-get_pos_lh1_fig16(4)-leg1_axis_offset_y_fig16, get_pos_lh1title_width_fig16, 0]);
set(lh1titleBottom2_fig16,'Position',[left_pos_centimeter_ax2a_fig16+leg1_axis_offset_x_fig16, top_pos_centimeter_ax2a_fig16-get_pos_an_fig16(4)-get_pos_lh1title_height_fig16-leg1_axis_offset_y_fig16, get_pos_lh1title_width_fig16, 0]);

%% Reposition Second legend [H/D_{V}, dependent color order]
%---set legend position---[left bottom width height]
drawnow
refreshdata
set(lh2_fig16,'Position',[left_pos_centimeter_ax2a_fig16 + get_pos_lh1_fig16(3) + leg1_axis_offset_x_fig16, top_pos_centimeter_ax2a_fig16 - get_pos_an_fig16(4) - get_pos_lh2_fig16(4) - leg1_axis_offset_y_fig16, get_pos_lh2_fig16(3), get_pos_lh2_fig16(4)]);

%---set annotation line position---[x_begin y_begin dx dy]
set(lh2titleTop_fig16,'Position',[left_pos_centimeter_ax2a_fig16 + get_pos_lh1title_width_fig16 + leg1_axis_offset_x_fig16, top_pos_centimeter_ax2a_fig16 - get_pos_an_fig16(4) - leg1_axis_offset_y_fig16, get_pos_lh2_fig16(3), 0]);
set(lh2titleBottom_fig16,'Position',[left_pos_centimeter_ax2a_fig16 + get_pos_lh1title_width_fig16 + leg1_axis_offset_x_fig16, top_pos_centimeter_ax2a_fig16 - get_pos_an_fig16(4) - get_pos_lh2_fig16(4) - leg1_axis_offset_y_fig16, get_pos_lh2_fig16(3), 0]);
set(lh2titleBottom2_fig16,'Position',[left_pos_centimeter_ax2a_fig16 + get_pos_lh1title_width_fig16 + leg1_axis_offset_x_fig16, top_pos_centimeter_ax2a_fig16 - get_pos_an_fig16(4) - get_pos_lh2title_height_fig16 - leg1_axis_offset_y_fig16, get_pos_lh2_fig16(3), 0]);

%% Reposition Third legend [W/D_{V}]
%---set legend position---[left bottom width height]
drawnow
refreshdata
set(lh3_fig16,'Position',[left_pos_centimeter_ax2a_fig16 + get_pos_lh1_fig16(3) + get_pos_lh2_fig16(3) + leg1_axis_offset_x_fig16, top_pos_centimeter_ax2a_fig16 - get_pos_an_fig16(4) - get_pos_lh3_fig16(4) - leg1_axis_offset_y_fig16, get_pos_lh3_fig16(3), get_pos_lh3_fig16(4)]);

%---set annotation line position---[x_begin y_begin dx dy]
set(lh3titleTop_fig16,'Position',[left_pos_centimeter_ax2a_fig16 + get_pos_lh1title_width_fig16 + get_pos_lh2_fig16(3) + leg1_axis_offset_x_fig16, top_pos_centimeter_ax2a_fig16 - get_pos_an_fig16(4) - leg1_axis_offset_y_fig16, get_pos_lh3_fig16(3), 0]);
set(lh3titleBottom_fig16,'Position',[left_pos_centimeter_ax2a_fig16 + get_pos_lh1title_width_fig16 + get_pos_lh2_fig16(3) + leg1_axis_offset_x_fig16, top_pos_centimeter_ax2a_fig16 - get_pos_an_fig16(4) - get_pos_lh3_fig16(4) - leg1_axis_offset_y_fig16, get_pos_lh3_fig16(3), 0]);
set(lh3titleBottom2_fig16,'Position',[left_pos_centimeter_ax2a_fig16 + get_pos_lh1title_width_fig16 + get_pos_lh2_fig16(3) + leg1_axis_offset_x_fig16, top_pos_centimeter_ax2a_fig16 - get_pos_an_fig16(4) - get_pos_lh3title_height_fig16 - leg1_axis_offset_y_fig16, get_pos_lh3_fig16(3), 0]);

%% Reposition Fourth legend [Solvent]
%---set legend position---[left bottom width height]
drawnow
refreshdata
set(lh4_fig16,'Position',[left_pos_centimeter_ax3_fig16+leg1_axis_offset_x_fig16, top_pos_centimeter_ax3_fig16-leg1_axis_offset_y_fig16, get_pos_lh4_fig16(3), get_pos_lh4_fig16(4)]);
drawnow
refreshdata
set(lh4_fig16,'Position',[left_pos_centimeter_ax3_fig16+leg1_axis_offset_x_fig16, top_pos_centimeter_ax3_fig16-leg1_axis_offset_y_fig16, get_pos_lh4_fig16(3), get_pos_lh4_fig16(4)]);

%---set annotation line position---[x_begin y_begin dx dy]
set(lh4titleTop_fig16,'Position',[left_pos_centimeter_ax3_fig16+leg1_axis_offset_x_fig16, top_pos_centimeter_ax3_fig16+get_pos_lh4_fig16(4)-leg1_axis_offset_y_fig16, get_pos_lh4_fig16(3), 0]);
set(lh4titleBottom_fig16,'Position',[left_pos_centimeter_ax3_fig16+leg1_axis_offset_x_fig16, top_pos_centimeter_ax3_fig16-get_pos_lh4_fig16(4)+get_pos_lh4_fig16(4)-leg1_axis_offset_y_fig16, get_pos_lh4_fig16(3), 0]);
set(lh4titleBottom2_fig16,'Position',[left_pos_centimeter_ax3_fig16+leg1_axis_offset_x_fig16, top_pos_centimeter_ax3_fig16-get_pos_lh4title_height_fig16+get_pos_lh4_fig16(4)-leg1_axis_offset_y_fig16, get_pos_lh4_fig16(3), 0]);

%% Figure 17) ModifiedArchimedesNumberUsingD_V [mParameter] vs. ModifiedParticleReynoldsNumberUsingU_avgD_VCRIT [rParameter] with non-linear regression
fig17 = figure(17);
ax1_fig17 = gca;

%---[Get Figure Properties]---
figProps = getDefaultFigProperties('2Column',1,1,1,1);

%---[Set Figure Size]--- 
set(fig17,'Units','centimeters','Position',[figProps.Figure.Position.Left figProps.Figure.Position.Bottom figProps.Figure.Position.Width figProps.Figure.Position.Width*figProps.Figure.Position.HWRatio])
pos = get(fig17,'Position');
set(fig17,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%% Drawing original points as scatter markers
for m = 1:refSize
    pOriginal(m) = scatter(paramStruct.mParameter.ValueSorted(m),paramStruct.rParameter.ValueSorted(m),figProps.Scatter.Markers.SizeData);
    set(pOriginal(m),{'Marker','MarkerFaceColor','MarkerFaceAlpha'},{lineSpec(m,1),lineSpec(m,2),figProps.Scatter.Markers.MarkerFaceAlpha})
    hold on
end
set(pOriginal,'MarkerEdgeColor',figProps.Scatter.Markers.MarkerEdgeColor)

%% Set figure properties

%----[Markers]----
set(p,'MarkerEdgeColor',figProps.Line.Markers.MarkerEdgeColor)
% set(p,'LineStyle','none') %Can't change 'LineStyle' because plotting each individual point in a loop (rather than providing access to whole array of data to connect points)

%---[Axis properties]---
%----[Font]----
ax1_fig17.FontSize = figProps.Axis.Font.FontSize;
%----[Ticks]----
ax1_fig17.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;
%----[Rulers]----
xlim([10 3000])
xticks(logspace(-3,4,8));
ylim([3 100])
yticks(logspace(0,3,4));

ax1_fig17.XAxis.Exponent = 0; %Scientific notation
xtickformat('%.0f') %Scientific notation
ax1_fig17.YAxis.Exponent = 0; %Scientific notation
ytickformat('%.0f') %Scientific notation
set(gca,'XScale','log');
set(gca,'YScale','log');
%----[Grids]----
grid on
ax1_fig17.GridLineStyle = figProps.Axis.Grids.GridLineStyle;
ax1_fig17.GridLineWidth = figProps.Axis.Grids.GridLineWidth;
ax1_fig17.GridColor = figProps.Axis.Grids.GridColor;
ax1_fig17.MinorGridLineStyle = figProps.Axis.Grids.MinorGridLineStyle;
ax1_fig17.MinorGridLineWidth = figProps.Axis.Grids.MinorGridLineWidth;
ax1_fig17.MinorGridColor = figProps.Axis.Grids.MinorGridColor;
%----[Labels]----
xlabel(paramStruct.mParameter.NewStr,'Interpreter',figProps.Axis.Labels.LabelInterpreter)

ax1_fig17.XLabel.FontSize = figProps.Axis.Labels.LabelFontSize;
ax1_fig17.XLabel.FontWeight = figProps.Axis.Labels.LabelFontWeight;

ylabel({'Modified critical';'particle Reynolds number, $\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}^*$'},'Interpreter',figProps.Axis.Labels.LabelInterpreter)

ax1_fig17.YLabel.FontSize = figProps.Axis.Labels.LabelFontSize;
ax1_fig17.YLabel.FontWeight = figProps.Axis.Labels.LabelFontWeight;

%----[Box Styling]----
set(ax1_fig17, 'Box', 'Off')

%% Fit c and m in power law model relationship at each HW/(D_{V})^2
dielectricConstant = [78.39,19.92,32.7]; %water, isopropanol, methanol
ETrio = repmat(dielectricConstant,10,1);

HW_DVSqTrio = zeros(10,3); %Preallocate
H_DVTrio = zeros(10,3); %Preallocate
W_DVTrio = zeros(10,3); %Preallocate
DV_Trio = zeros(10,3); %Preallocate

KalmanArTrio = zeros(10,3); % Preallocate
KalmanReTrio = zeros(10,3); % Preallocate
KalmanArTrioLog = zeros(10,3); % Preallocate
KalmanReTrioLog = zeros(10,3); % Preallocate

pFit = zeros(10,3); % Preallocate for second-order polynomial
yFit = zeros(10,3); % Preallocate for second-order polynomial
yFitLog = zeros(10,3); % Preallocate for second-order polynomial

for i = 1:10  
    HW_DVSqTrio(i,:) = paramStruct.zParameterHW_D_maxSq.ValueSorted(i:10:end);
    H_DVTrio(i,:) = paramStruct.hD_maxParameter.ValueSorted(i:10:end);
    W_DVTrio(i,:) = paramStruct.wD_maxParameter.ValueSorted(i:10:end);

    KalmanArTrio(i,:) = paramStruct.mParameter.ValueSorted(i:10:end);
    KalmanArTrioLog(i,:) = log10(KalmanArTrio(i,:));

    KalmanReTrio(i,:) = paramStruct.rParameter.ValueSorted(i:10:end);
    KalmanReTrioLog(i,:) = log10(KalmanReTrio(i,:));
end

ETrio(:, [1, 2]) = ETrio(:, [2, 1]);
KalmanArTrio(:, [1, 2]) = KalmanArTrio(:, [2, 1]);
KalmanReTrio(:, [1, 2]) = KalmanReTrio(:, [2, 1]);
HW_DVSqTrio(:, [1, 2]) = HW_DVSqTrio(:, [2, 1]);
H_DVTrio(:, [1, 2]) = H_DVTrio(:, [2, 1]);
W_DVTrio(:, [1, 2]) = W_DVTrio(:, [2, 1]);

KalmanArTrioLog(:, [1, 2]) = KalmanArTrioLog(:, [2, 1]);
KalmanReTrioLog(:, [1, 2]) = KalmanReTrioLog(:, [2, 1]);

%% Finding modifier for Ar
modelChoice = 'modelFinal'; %adapted from modelFour
% modelChoice = 'modelOne';
% modelChoice = 'modelTwo';
% modelChoice = 'modelThree';
% modelChoice = 'modelFour'; %CURRENT B2K MODEL
% modelChoice = 'modelFive';
% modelChoice = 'modelSix';
% modelChoice = 'modelSeven';
% modelChoice = 'modelEight';
% modelChoice = 'modelNine';
% modelChoice = 'modelTen';

%% ---[Models]---
% Model F)  [4-Parameters] 2.7 .* ( a.*(Ar_star,d_V).^(b) .* (d_V/W).^(c) .* (E).^(d) ).^(3/7)
% Model 1)  [8-Parameters] 2.7 .* ( a.*(Ar_star,d_V).^(b) .* c.*(H/W).^(d) .* e.*(d_V/W).^(f) .* g.*(E).^(h) ).^(3/7)
% Model 2)  [4-Parameters] 2.7 .* ( (Ar_star,d_V).^(a) .* (H/W).^(b) .* (d_V/W).^(c) .* (E).^(d) ).^(3/7)
% Model 3)  [3-Parameters] 2.7 .* ( (Ar_star,d_V).^(a) .* (d_V/W).^(b) .* (E).^(c) ).^(3/7)
% Model 4)  [4-Parameters] 2.7 .* ( a.*(Ar_star,d_V).^(b) .* (d_V/W).^(c) .* (E).^(d) ).^(3/7)
% Model 5)  [5-Parameters] 2.7 .* ( a.*(Ar_star,d_V).^(b) .* (H/W).^(c) .* (d_V/W).^(d) .* (E).^(e) ).^(3/7)
% Model 6)  [5-Parameters] 2.7 .* ( (Ar_star,d_V).^(a) .* b.*(H/W).^(c) .* (d_V/W).^(d) .* (E).^(e) ).^(3/7)
% Model 7)  [5-Parameters] 2.7 .* ( (Ar_star,d_V).^(a) .* (H/W).^(b) .* c.*(d_V/W).^(d) .* (E).^(e) ).^(3/7)
% Model 8)  [5-Parameters] 2.7 .* ( (Ar_star,d_V).^(a) .* (H/W).^(b) .* (d_V/W).^(c) .* d.*(E).^(e) ).^(3/7)
% Model 9)  [4-Parameters] 2.7 .* ( a.*(Ar_star,d_V).^(b) .* (H/W).^(c) .* (d_V/W).^(d) ).^(3/7)
% Model 10) [6-Parameters] a .* ( (Ar_star,d_V).^(b) .* (H/W).^(c) .* (d_V/W).^(d) .* (E).^(e) ).^(f)

% Assemble data
Ar_star = reshape(KalmanArTrio, [], 1);
H_DV = reshape(H_DVTrio, [], 1);
W_DV = reshape(W_DVTrio, [], 1);
E = reshape(ETrio, [], 1);

% Convert HW_D2Trio to a single column array
HW_DVSqTrio_reshaped = reshape(HW_DVSqTrio, [], 1);
% Round HW_D2Trio array 1 decimal place
HW_DVSqTrio_reshaped_rounded = round(HW_DVSqTrio_reshaped, 1);

Re = reshape(KalmanReTrio, [], 1);

% Model via anonymous function
switch modelChoice
    case 'modelFinal'
        data = [Ar_star 1./W_DV E];
        model = @(beta, data)2.7 .*(beta(1) .* data(:,1).^(beta(2)) .* data(:,2).^(beta(3)) .* data(:,3).^(beta(4))).^(3/7);
        % Initial guesses for parameters: [a, b, c, d];
        beta0 = [1 1 1 1];        
    case 'modelOne'
        data = [Ar_star H_DV./W_DV 1./W_DV E];
        model = @(beta, data)2.7 .*(beta(1).* data(:,1).^(beta(2)) .* beta(3).*data(:,2).^(beta(4)) .* beta(5).*data(:,3).^(beta(6)) .* beta(7).*data(:,4).^(beta(8)) ).^(3/7);
        % Initial guesses for parameters: [a, b, c, d, e, f, g, h];
        beta0 = [1 1 1 1 1 1 1 1];
    case 'modelTwo'
        data = [Ar_star H_DV./W_DV 1./W_DV E];
        model = @(beta, data)2.7 .*(data(:,1).^(beta(1)) .* data(:,2).^(beta(2)) .* data(:,3).^(beta(3)) .* data(:,4).^(beta(4)) ).^(3/7);
        % Initial guesses for parameters: [a, b, c, d];
        beta0 = [1 1 1 1];    
    case 'modelThree'
        data = [Ar_star 1./W_DV E];
        model = @(beta, data)2.7 .*(data(:,1).^(beta(1)) .* data(:,2).^(beta(2)) .* data(:,3).^(beta(3))).^(3/7);
        % Initial guesses for parameters: [a, b, c];
        beta0 = [1 1 1];
    case 'modelFour'
        data = [Ar_star 1./W_DV E];
        model = @(beta, data)2.7 .*(beta(1) .* data(:,1).^(beta(2)) .* data(:,2).^(beta(3)) .* data(:,3).^(beta(4))).^(3/7);
        % Initial guesses for parameters: [a, b, c, d];
        beta0 = [1 1 1 1];
    case 'modelFive'
        data = [Ar_star H_DV./W_DV 1./W_DV E];
        model = @(beta, data)2.7 .*(beta(1).* data(:,1).^(beta(2)) .* data(:,2).^(beta(3)) .* data(:,3).^(beta(4)) .* data(:,4).^(beta(5)) ).^(3/7);
        % Initial guesses for parameters: [a, b, c, d, e];
        beta0 = [1 1 1 1 1];
    case 'modelSix'
        data = [Ar_star H_DV./W_DV 1./W_DV E];
        model = @(beta, data)2.7 .*(data(:,1).^(beta(1)) .* beta(2).*data(:,2).^(beta(3)) .* data(:,3).^(beta(4)) .* data(:,4).^(beta(5)) ).^(3/7);
        % Initial guesses for parameters: [a, b, c, d, e];
        beta0 = [1 1 1 1 1];
    case 'modelSeven'
        data = [Ar_star H_DV./W_DV 1./W_DV E];
        model = @(beta, data)2.7 .*(data(:,1).^(beta(1)) .* data(:,2).^(beta(2)) .* beta(3).*data(:,3).^(beta(4)) .* data(:,4).^(beta(5)) ).^(3/7);
        % Initial guesses for parameters: [a, b, c, d, e];
        beta0 = [1 1 1 1 1];
    case 'modelEight'
        data = [Ar_star H_DV./W_DV 1./W_DV E];
        model = @(beta, data)2.7 .*(data(:,1).^(beta(1)) .* data(:,2).^(beta(2)) .* data(:,3).^(beta(3)) .* beta(4).*data(:,4).^(beta(5)) ).^(3/7);
        % Initial guesses for parameters: [a, b, c, d, e];
        beta0 = [1 1 1 1 1];
    case 'modelNine'
        data = [Ar_star H_DV./W_DV 1./W_DV];
        model = @(beta, data)2.7 .*(beta(1).*data(:,1).^(beta(2)) .* data(:,2).^(beta(3)) .* data(:,3).^(beta(4)) ).^(3/7);
        % Initial guesses for parameters: [a, b, c, d];
        beta0 = [1 1 1 1];
    case 'modelTen'
        data = [Ar_star H_DV./W_DV 1./W_DV E];
        model = @(beta, data)beta(1) .*(data(:,1).^(beta(2)) .* data(:,2).^(beta(3)) .* data(:,3).^(beta(4)) .* data(:,4).^(beta(5)) ).^(beta(6));
        % Initial guesses for parameters: [a, b, c, d, e, f];
        beta0 = [1 1 1 1 1 1];
end

[beta, res, J, Sigma, mse] = nlinfit(data, Re, model, beta0);

disp('beta')
disp(beta)
disp('res')
disp(res)
disp('J')
disp(J)
disp('Sigma')
disp(Sigma)
disp('mse')
disp(mse)

% Calculate standard errors
SE = sqrt(diag(Sigma));

% Define confidence level (e.g., 95% confidence interval)
alpha = 0.05;

% Compute confidence intervals
par_ci = nlparci(beta, res, 'Jacobian', J, 'Alpha', alpha);

% Compute prediction intervals
[ypred, delta] = nlpredci(model, data, beta, res, 'Jacobian', J, 'Alpha', alpha);

% Define the degrees of freedom (number of observations minus number of parameters)
dof = length(data) - length(beta);

% Critical t-value for the specified confidence level
t_crit = tinv(1 - alpha/2, dof);

% Calculate standard errors of the predictions
SE_pred = delta / t_crit;
% Calculate average standard error
average_SE_pred = mean(SE_pred);
disp('average_SE_pred')
disp(average_SE_pred)

% ypred: Predicted responses
% delta: Half-widths of the prediction intervals

% Number of parameters
num_params = length(beta);

% Display parameter estimates
fprintf('Parameter Estimates and Uncertainties:\n');
for i = 1:num_params
    fprintf('Parameter %d:\n', i);
    fprintf('  Estimate       = %.4f\n', beta(i));
    fprintf('  Standard Error = %.4f\n', SE(i));
    fprintf('  95%% CI         = [%.4f, %.4f]\n', par_ci(i, 1), par_ci(i, 2));
end

% Combine data into a table
results_table = table(data, Re, ypred, ypred - delta, ypred + delta, ...
    'VariableNames', {'Data', 'Observed', 'Predicted', 'Lower_CI', 'Upper_CI'});

% Display the table
disp(results_table);

%%
% Define the range of Ar values for the plot
Arplot = (133:1:1971)'; % Ensure Arplot is a column vector

% Unique E and Ar values
E_values = [19.92; 78.39; 32.70];
Ar_values = [133.5510; 793.3536; 1970.9811];

% Ensure Ar_values is sorted in increasing order
[Ar_values_sorted, sort_idx] = sort(Ar_values);
E_values_sorted = E_values(sort_idx);

% Interpolate E for the range of Ar values
E_interp = interp1(Ar_values_sorted, E_values_sorted, Arplot, 'linear', 'extrap');

data_plot = cell(10,1);
Re_plot = cell(10,1);

hold on
hold on

switch modelChoice
    case 'modelFinal'
        for I_fixed = 1:10
            % Use the interpolated E values
            E_rep = E_interp; % Vector of E values corresponding to Arplot

            % Construct the data matrix for plotting
            data_plot{I_fixed,1} = [Arplot, repmat(1./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), E_rep];

            % Evaluate the model with the fitted parameters
            Re_plot{I_fixed,1} = model(beta, data_plot{I_fixed,1});

            plot(Arplot, Re_plot{I_fixed,1}, 'Color', paramStruct.plotData.arrayColor(I_fixed), 'LineWidth', 2);
        end
    case 'modelOne'
        for I_fixed = 1:10
            % Use the interpolated E values
            E_rep = E_interp; % Vector of E values corresponding to Arplot

            % Construct the data matrix for plotting
            data_plot{I_fixed,1} = [Arplot, repmat((paramStruct.hD_maxParameter.ValueSorted(I_fixed))./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), repmat(1./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), E_rep];

            % Evaluate the model with the fitted parameters
            Re_plot{I_fixed,1} = model(beta, data_plot{I_fixed,1});

            plot(Arplot, Re_plot{I_fixed,1}, 'Color', paramStruct.plotData.arrayColor(I_fixed), 'LineWidth', 2);
        end
    case 'modelTwo'
        for I_fixed = 1:10
            % Use the interpolated E values
            E_rep = E_interp; % Vector of E values corresponding to Arplot

            % Construct the data matrix for plotting
            data_plot{I_fixed,1} = [Arplot, repmat((paramStruct.hD_maxParameter.ValueSorted(I_fixed))./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), repmat(1./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), E_rep];

            % Evaluate the model with the fitted parameters
            Re_plot{I_fixed,1} = model(beta, data_plot{I_fixed,1});

            plot(Arplot, Re_plot{I_fixed,1}, 'Color', paramStruct.plotData.arrayColor(I_fixed), 'LineWidth', 2);        
        end
    case 'modelThree'
        for I_fixed = 1:10
            % Use the interpolated E values
            E_rep = E_interp; % Vector of E values corresponding to Arplot

            % Construct the data matrix for plotting
            data_plot{I_fixed,1} = [Arplot, repmat(1./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), E_rep];

            % Evaluate the model with the fitted parameters
            Re_plot{I_fixed,1} = model(beta, data_plot{I_fixed,1});

            plot(Arplot, Re_plot{I_fixed,1}, 'Color', paramStruct.plotData.arrayColor(I_fixed), 'LineWidth', 2);
        end
    case 'modelFour'
        for I_fixed = 1:10
            % Use the interpolated E values
            E_rep = E_interp; % Vector of E values corresponding to Arplot

            % Construct the data matrix for plotting
            data_plot{I_fixed,1} = [Arplot, repmat(1./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), E_rep];

            % Evaluate the model with the fitted parameters
            Re_plot{I_fixed,1} = model(beta, data_plot{I_fixed,1});

            plot(Arplot, Re_plot{I_fixed,1}, 'Color', paramStruct.plotData.arrayColor(I_fixed), 'LineWidth', 2);
        end
    case 'modelFive'
        for I_fixed = 1:10
            % Use the interpolated E values
            E_rep = E_interp; % Vector of E values corresponding to Arplot

            % Construct the data matrix for plotting
            data_plot{I_fixed,1} = [Arplot, repmat((paramStruct.hD_maxParameter.ValueSorted(I_fixed))./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), repmat(1./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), E_rep];

            % Evaluate the model with the fitted parameters
            Re_plot{I_fixed,1} = model(beta, data_plot{I_fixed,1});

            plot(Arplot, Re_plot{I_fixed,1}, 'Color', paramStruct.plotData.arrayColor(I_fixed), 'LineWidth', 2);
        end
    case 'modelSix'
        for I_fixed = 1:10
            % Use the interpolated E values
            E_rep = E_interp; % Vector of E values corresponding to Arplot

            % Construct the data matrix for plotting
            data_plot{I_fixed,1} = [Arplot, repmat((paramStruct.hD_maxParameter.ValueSorted(I_fixed))./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), repmat(1./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), E_rep];

            % Evaluate the model with the fitted parameters
            Re_plot{I_fixed,1} = model(beta, data_plot{I_fixed,1});

            plot(Arplot, Re_plot{I_fixed,1}, 'Color', paramStruct.plotData.arrayColor(I_fixed), 'LineWidth', 2);
        end
    case 'modelSeven'
        for I_fixed = 1:10
            % Use the interpolated E values
            E_rep = E_interp; % Vector of E values corresponding to Arplot

            % Construct the data matrix for plotting
            data_plot{I_fixed,1} = [Arplot, repmat((paramStruct.hD_maxParameter.ValueSorted(I_fixed))./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), repmat(1./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), E_rep];

            % Evaluate the model with the fitted parameters
            Re_plot{I_fixed,1} = model(beta, data_plot{I_fixed,1});

            plot(Arplot, Re_plot{I_fixed,1}, 'Color', paramStruct.plotData.arrayColor(I_fixed), 'LineWidth', 2);
        end
    case 'modelEight'
        for I_fixed = 1:10
            % Use the interpolated E values
            E_rep = E_interp; % Vector of E values corresponding to Arplot

            % Construct the data matrix for plotting
            data_plot{I_fixed,1} = [Arplot, repmat((paramStruct.hD_maxParameter.ValueSorted(I_fixed))./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), repmat(1./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), E_rep];

            % Evaluate the model with the fitted parameters
            Re_plot{I_fixed,1} = model(beta, data_plot{I_fixed,1});

            plot(Arplot, Re_plot{I_fixed,1}, 'Color', paramStruct.plotData.arrayColor(I_fixed), 'LineWidth', 2);
        end
    case 'modelNine'
        for I_fixed = 1:10
            % Construct the data matrix for plotting
            data_plot{I_fixed,1} = [Arplot, repmat((paramStruct.hD_maxParameter.ValueSorted(I_fixed))./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), repmat(1./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1)];

            % Evaluate the model with the fitted parameters
            Re_plot{I_fixed,1} = model(beta, data_plot{I_fixed,1});

            plot(Arplot, Re_plot{I_fixed,1}, 'Color', paramStruct.plotData.arrayColor(I_fixed), 'LineWidth', 2);
        end
    case 'modelTen'
        for I_fixed = 1:10
            % Use the interpolated E values
            E_rep = E_interp; % Vector of E values corresponding to Arplot

            % Construct the data matrix for plotting
            data_plot{I_fixed,1} = [Arplot, repmat((paramStruct.hD_maxParameter.ValueSorted(I_fixed))./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), repmat(1./(paramStruct.wD_maxParameter.ValueSorted(I_fixed)), length(Arplot), 1), E_rep];

            % Evaluate the model with the fitted parameters
            Re_plot{I_fixed,1} = model(beta, data_plot{I_fixed,1});

            plot(Arplot, Re_plot{I_fixed,1}, 'Color', paramStruct.plotData.arrayColor(I_fixed), 'LineWidth', 2);
        end
end

%% Zero legend [Title: Flat-plate, experimental]
yShiftAnnotation_fig17 = 0.03;
anPositionX_fig17 = 0.1446;
anPositionY_fig17 = 0.9017-yShiftAnnotation_fig17;
anPositionWidth_fig17 = 0.4122;
anPositionHeight_fig17 = 0.05;
an_fig17 = annotation('textbox','String','Flat-plate, experimental','FontSize',figProps.Legend.Labels.Title.FontSize,'HorizontalAlignment','center','BackgroundColor',[1 1 1],'EdgeColor',[1 1 1],'FaceAlpha',1,'Interpreter','latex','Position',[anPositionX_fig17 anPositionY_fig17 anPositionWidth_fig17 anPositionHeight_fig17]); %[0.1446 0.4602 0.3 0.05] %0.3872
% anTop_fig17 = annotation('line', [anPositionX_fig17, anPositionX_fig17 + anPositionWidth_fig17], [anPositionY_fig17 + anPositionHeight_fig17, anPositionY_fig17 + anPositionHeight_fig17], 'Color', 'black');
anTop_fig17 = annotation('line', [anPositionX_fig17, anPositionX_fig17 + anPositionWidth_fig17 + 0.0213], [anPositionY_fig17 + anPositionHeight_fig17, anPositionY_fig17 + anPositionHeight_fig17], 'Color', 'black');
anBottom_fig17 = annotation('line', [anPositionX_fig17, anPositionX_fig17 + anPositionWidth_fig17], [anPositionY_fig17,anPositionY_fig17], 'Color', 'black');

%% First legend [HW/(d_{V})^2, set progressive color order]
%---[Legends strings - Dynamic]---
zLabelStrRounded = B2KFormatValuesToStrOneSigFigMaxUnc(paramStruct.zParameterHW_D_VSq.ValueSorted,paramStruct.zParameterHW_D_VSq.UncertaintySorted);

% Create dummy plot to generate first legend since default legend for scatter cannot modify legend marker size
for m = 1:refSize
    pDummy(m) = plot(NaN,NaN,'MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none');
    set(pDummy(m),{'Marker','MarkerFaceColor'},{lineSpec(m,1),[lineSpec(m,2)]})
    hold on
end
set(pDummy,'MarkerEdgeColor',figProps.Line.Markers.MarkerEdgeColor)

% Plot first legend
lh1_fig17 = legend(gca,pDummy(1+refSizeEach:refSizeEach*2),zLabelStrRounded,'FontSize',figProps.Legend.Fonts.FontSize,'NumColumns',1,'Interpreter',figProps.Legend.Interpreter); %circle marker
lh1_fig17.Units = 'normalized';
lh1_fig17.Location = 'northwest';
lh1_fig17_title = title(lh1_fig17,paramStruct.zParameterHW_D_VSq.LegendStr,'FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
lh1_fig17.EdgeColor = figProps.Legend.ColorAndStyling.EdgeColor;

% Adjust horizontal spacing for legend lh1
lh1_fig17.ItemTokenSize(1) = figProps.Legend.Fonts.ItemTokenSize;

% Set position of first legend
lh1_OrigPos_fig17 = get(lh1_fig17,'Position');
set(lh1_fig17,'Position',[lh1_OrigPos_fig17(1), lh1_OrigPos_fig17(2)-yShiftAnnotation_fig17, lh1_OrigPos_fig17(3), lh1_OrigPos_fig17(4)]);
lh1titleTop_fig17 = annotation('line', [lh1_OrigPos_fig17(1), lh1_OrigPos_fig17(1) + lh1_OrigPos_fig17(3)], [lh1_OrigPos_fig17(2)-yShiftAnnotation_fig17*2.6 + lh1_OrigPos_fig17(4), lh1_OrigPos_fig17(2)-yShiftAnnotation_fig17*2.6 + lh1_OrigPos_fig17(4)], 'Color', 'black');
lh1titleBottom_fig17 = annotation('line', [lh1_OrigPos_fig17(1), lh1_OrigPos_fig17(1) + lh1_OrigPos_fig17(3)], [lh1_OrigPos_fig17(2)-yShiftAnnotation_fig17,lh1_OrigPos_fig17(2)-yShiftAnnotation_fig17], 'Color', 'black');
lh1titleBottom2_fig17 = annotation('line', [lh1_OrigPos_fig17(1), lh1_OrigPos_fig17(1) + lh1_OrigPos_fig17(3)], [lh1_OrigPos_fig17(2)-yShiftAnnotation_fig17 + lh1_OrigPos_fig17(4), lh1_OrigPos_fig17(2)-yShiftAnnotation_fig17 + lh1_OrigPos_fig17(4)], 'Color', 'black');

%% Second legend [H/d_{V}, dependent color order] [TESTING]
% Copy the axes and plot the second legend
ax2_fig17 = copyobj(ax1_fig17,gcf); %Copies ah1 object (but duplicates data points)
delete(get(ax2_fig17,'Children')) %Deletes duplicate data points

% Remove second axis visibility 
set(ax2_fig17, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

%---Legends strings
zLabelStrRounded_lh2_fig17 = B2KFormatValuesToStrOneSigFigMaxUnc(paramStruct.hD_VParameter.ValueSorted,paramStruct.hD_VParameter.UncertaintySorted);

% Plot second legend
hold on
lh2dump1_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig17);
lh2dump2_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig17);
lh2dump3_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig17);
lh2dump4_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig17);
lh2dump5_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig17);
lh2dump6_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig17);
lh2dump7_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig17);
lh2dump8_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig17);
lh2dump9_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig17);
lh2dump10_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig17);
hold off

lh2_fig17 = legend(ax2_fig17,[lh2dump1_fig17 lh2dump2_fig17 lh2dump3_fig17 lh2dump4_fig17 lh2dump5_fig17 lh2dump6_fig17 lh2dump7_fig17 lh2dump8_fig17 lh2dump9_fig17 lh2dump10_fig17],zLabelStrRounded_lh2_fig17,'FontSize',figProps.Legend.Fonts.FontSize,'NumColumns',1,'Interpreter',figProps.Legend.Interpreter);
lh2_fig17.Units = 'normalized';
lh2_fig17.Location = 'northwest';
lh2title_fig17 = title(lh2_fig17,paramStruct.hD_VParameter.LegendStr,'FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
lh2_fig17.EdgeColor = figProps.Legend.ColorAndStyling.EdgeColor;

% Set background color of second legend to white (instead of grey by default)
set(lh2_fig17,'Color',get(ax1_fig2,'Color'));

% Adjust horizontal spacing for legend lh2
lh2_fig17.ItemTokenSize(1) = 0; %figProps.Legend.Fonts.ItemTokenSize

% Set position relative to first legend
lh1_OrigPos_fig17 = get(lh1_fig17,'Position');
lh2_OrigPos_fig17 = get(lh2_fig17,'Position');
set(lh2_fig17,'Position',[lh2_OrigPos_fig17(1)+lh1_OrigPos_fig17(3), lh2_OrigPos_fig17(2)-yShiftAnnotation_fig17, lh2_OrigPos_fig17(3), lh2_OrigPos_fig17(4)]);
lh2titleTop_fig17 = annotation('line', [lh2_OrigPos_fig17(1)+lh1_OrigPos_fig17(3), lh2_OrigPos_fig17(1)+lh1_OrigPos_fig17(3) + lh2_OrigPos_fig17(3)], [lh2_OrigPos_fig17(2)-yShiftAnnotation_fig17*2.6 + lh2_OrigPos_fig17(4), lh2_OrigPos_fig17(2)-yShiftAnnotation_fig17*2.6 + lh2_OrigPos_fig17(4)], 'Color', 'black');
lh2titleBottom_fig17 = annotation('line', [lh2_OrigPos_fig17(1)+lh1_OrigPos_fig17(3), lh2_OrigPos_fig17(1)+lh1_OrigPos_fig17(3) + lh2_OrigPos_fig17(3)], [lh2_OrigPos_fig17(2)-yShiftAnnotation_fig17,lh2_OrigPos_fig17(2)-yShiftAnnotation_fig17], 'Color', 'black');
lh2titleBottom2_fig17 = annotation('line', [lh2_OrigPos_fig17(1)+lh1_OrigPos_fig17(3), lh2_OrigPos_fig17(1)+lh1_OrigPos_fig17(3) + lh2_OrigPos_fig17(3)], [lh2_OrigPos_fig17(2)-yShiftAnnotation_fig17 + lh2_OrigPos_fig17(4), lh2_OrigPos_fig17(2)-yShiftAnnotation_fig17 + lh2_OrigPos_fig17(4)], 'Color', 'black');

%% Third legend [AR] [TESTING]
% Copy the axes and plot the third legend
ax3_fig17 = copyobj(ax2_fig17,gcf); %Copies ah1 object (but duplicates data points)
delete(get(ax3_fig17,'Children')) %Deletes duplicate data points

% Remove third axis visibility 
set(ax3_fig17, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

%---Legends strings
% zLabelStrRounded_lh3_fig17 = {'0.67','0.67','1.50','1.00','0.67','1.00','1.33','1.00','1.33','1.33'};
zLabelStrRounded_lh3_fig17 = B2KFormatValuesToStrOneSigFigMaxUnc(paramStruct.wD_maxParameter.ValueSorted,paramStruct.wD_maxParameter.UncertaintySorted);

% Plot third legend
hold on
lh3dump1_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig17);
lh3dump2_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig17);
lh3dump3_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig17);
lh3dump4_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig17);
lh3dump5_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig17);
lh3dump6_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig17);
lh3dump7_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig17);
lh3dump8_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig17);
lh3dump9_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig17);
lh3dump10_fig17 = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig17);
hold off

lh3_fig17 = legend(ax3_fig17,[lh3dump1_fig17 lh3dump2_fig17 lh3dump3_fig17 lh3dump4_fig17 lh3dump5_fig17 lh3dump6_fig17 lh3dump7_fig17 lh3dump8_fig17 lh3dump9_fig17 lh3dump10_fig17],zLabelStrRounded_lh3_fig17,'FontSize',figProps.Legend.Fonts.FontSize,'NumColumns',1,'Interpreter','latex');
lh3_fig17.Units = 'normalized';
lh3_fig17.Location = 'northwest';
lh3title_fig17 = title(lh3_fig17,paramStruct.aParameter.LegendStr,'FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
lh3_fig17.EdgeColor = figProps.Legend.ColorAndStyling.EdgeColor;

% Set background color of third legend to white (instead of grey by default)
set(lh3_fig17,'Color',get(ax1_fig17,'Color'));

% Adjust horizontal spacing for legend lh3
lh3_fig17.ItemTokenSize(1) = 0; %figProps.Legend.Fonts.ItemTokenSize

% Set position relative to second legend
lh2_OrigPos_fig17 = get(lh2_fig17,'Position');
lh3_OrigPos_fig17 = get(lh3_fig17,'Position');
set(lh3_fig17,'Position',[lh3_OrigPos_fig17(1)+lh2_OrigPos_fig17(3)+lh1_OrigPos_fig17(3), lh3_OrigPos_fig17(2)-yShiftAnnotation_fig17, lh3_OrigPos_fig17(3), lh3_OrigPos_fig17(4)]);
lh3titleTop_fig17 = annotation('line', [lh3_OrigPos_fig17(1)+lh2_OrigPos_fig17(3)+lh1_OrigPos_fig17(3), lh3_OrigPos_fig17(1)+lh2_OrigPos_fig17(3)+lh1_OrigPos_fig17(3) + lh3_OrigPos_fig17(3)], [lh3_OrigPos_fig17(2)-yShiftAnnotation_fig17*2.6 + lh3_OrigPos_fig17(4), lh3_OrigPos_fig17(2)-yShiftAnnotation_fig17*2.6 + lh3_OrigPos_fig17(4)], 'Color', 'black');
lh3titleBottom_fig17 = annotation('line', [lh3_OrigPos_fig17(1)+lh2_OrigPos_fig17(3)+lh1_OrigPos_fig17(3), lh3_OrigPos_fig17(1)+lh2_OrigPos_fig17(3)+lh1_OrigPos_fig17(3) + lh3_OrigPos_fig17(3)], [lh3_OrigPos_fig17(2)-yShiftAnnotation_fig17,lh3_OrigPos_fig17(2)-yShiftAnnotation_fig17], 'Color', 'black');
lh2titleBottom2_fig17 = annotation('line', [lh3_OrigPos_fig17(1)+lh2_OrigPos_fig17(3)+lh1_OrigPos_fig17(3), lh3_OrigPos_fig17(1)+lh2_OrigPos_fig17(3)+lh1_OrigPos_fig17(3) + lh3_OrigPos_fig17(3)], [lh3_OrigPos_fig17(2)-yShiftAnnotation_fig17 + lh3_OrigPos_fig17(4), lh3_OrigPos_fig17(2)-yShiftAnnotation_fig17 + lh3_OrigPos_fig17(4)], 'Color', 'black');

%% Fourth legend [Solvent]
% Create array for markers and labels for fourth legend
lh4_fig17_arrayMarkerdump = ["o" "square" "diamond"];
lh4_fig17_typeSolvent = cellstr(["Isopropanol","Water","Methanol"]);

% Copy the axes and plot the fourth legend
ax4_fig17 = copyobj(ax3_fig17,gcf); %Copies ah1 object (but duplicates data points)
delete(get(ax4_fig17,'Children')) %Deletes duplicate data points

% Remove fourth axis visibility 
set(ax4_fig17, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

% Plot "dummy" data for fourth legend; tried using single variable with loop but does not work
hold on
lh4dump1_fig17 = plot(NaN,'Marker',lh4_fig17_arrayMarkerdump(1),'DisplayName',lh4_fig17_typeSolvent{1},'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none','Parent',ax4_fig17);
lh4dump2_fig17 = plot(NaN,'Marker',lh4_fig17_arrayMarkerdump(2),'DisplayName',lh4_fig17_typeSolvent{2},'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none','Parent',ax4_fig17);
lh4dump3_fig17 = plot(NaN,'Marker',lh4_fig17_arrayMarkerdump(3),'DisplayName',lh4_fig17_typeSolvent{3},'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none','Parent',ax4_fig17);
hold off

lh4_fig17 = legend(ax4_fig17,[lh4dump1_fig17 lh4dump2_fig17 lh4dump3_fig17],lh4_fig17_typeSolvent,'FontSize',figProps.Legend.Fonts.FontSize,'Interpreter',figProps.Legend.Interpreter);
lh4_fig17.Units = 'normalized';
lh4_fig17.Location = 'northwest';
lh4title_fig17 = title(lh4_fig17,'Solvent','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
lh4_fig17.EdgeColor = [1 1 1];

% Set background color of fourth legend to white (instead of grey by default)
set(lh4_fig17,'Color',get(ax1_fig2,'Color'));

% Adjust horizontal spacing for legend lh4
lh4_fig17.ItemTokenSize(1) = figProps.Legend.Fonts.ItemTokenSize;

% Set position relative to third legend
lh3_OrigPos_fig17 = get(lh3_fig17,'Position');
lh4_OrigPos_fig17 = get(lh4_fig17,'Position');
lh4_addLeftSpace_fig17 = 0.025;
set(lh4_fig17,'Position',[lh4_addLeftSpace_fig17+lh4_OrigPos_fig17(1)+lh3_OrigPos_fig17(3)+lh2_OrigPos_fig17(3)+lh1_OrigPos_fig17(3), lh4_OrigPos_fig17(2)-yShiftAnnotation_fig17, lh4_OrigPos_fig17(3), lh4_OrigPos_fig17(4)]);
lh4titleTop_fig17 = annotation('line', [lh4_addLeftSpace_fig17+lh4_OrigPos_fig17(1)+lh3_OrigPos_fig17(3)+lh2_OrigPos_fig17(3)+lh1_OrigPos_fig17(3), lh4_addLeftSpace_fig17+lh4_OrigPos_fig17(1)+lh3_OrigPos_fig17(3)+lh2_OrigPos_fig17(3)+lh1_OrigPos_fig17(3) + lh4_OrigPos_fig17(3)], [lh4_OrigPos_fig17(2)-yShiftAnnotation_fig17*2.6 + lh4_OrigPos_fig17(4), lh4_OrigPos_fig17(2)-yShiftAnnotation_fig17*2.6 + lh4_OrigPos_fig17(4)], 'Color', 'black');
lh4titleBottom_fig17 = annotation('line', [lh4_addLeftSpace_fig17+lh4_OrigPos_fig17(1)+lh3_OrigPos_fig17(3)+lh2_OrigPos_fig17(3)+lh1_OrigPos_fig17(3), lh4_addLeftSpace_fig17+lh4_OrigPos_fig17(1)+lh3_OrigPos_fig17(3)+lh2_OrigPos_fig17(3)+lh1_OrigPos_fig17(3) + lh4_OrigPos_fig17(3)], [lh4_OrigPos_fig17(2)-yShiftAnnotation_fig17,lh4_OrigPos_fig17(2)-yShiftAnnotation_fig17], 'Color', 'black');
lh4titleBottom2_fig17 = annotation('line', [lh4_addLeftSpace_fig17+lh4_OrigPos_fig17(1)+lh3_OrigPos_fig17(3)+lh2_OrigPos_fig17(3)+lh1_OrigPos_fig17(3), lh4_addLeftSpace_fig17+lh4_OrigPos_fig17(1)+lh3_OrigPos_fig17(3)+lh2_OrigPos_fig17(3)+lh1_OrigPos_fig17(3) + lh4_OrigPos_fig17(3)], [lh4_OrigPos_fig17(2)-yShiftAnnotation_fig17 + lh4_OrigPos_fig17(4), lh4_OrigPos_fig17(2)-yShiftAnnotation_fig17 + lh4_OrigPos_fig17(4)], 'Color', 'black');

%%
uistack(pOriginal,'up',30)

%% Q-Q Plot of Residuals
fig63 = figure(63);
qqplot(res);
title('Q-Q Plot of Residuals');

%% Figure 18) - Boxplot for nlinfit ignoring group-specific effects
fig18 = figure(18);
ax1_fig18 = gca;

%---[Get Figure Properties]---
figProps = getDefaultFigProperties('2Column',1,1,1,1);

%---[Set Figure Size]---
set(fig18,'Units','centimeters','Position',[3 3 figProps.Figure.Position.Width figProps.Figure.Position.Width*figProps.Figure.Position.HWRatio])
pos = get(fig18,'Position');
set(fig18,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

% Create dummy grouping variables for boxplot; need to create because using grouping variables directly causes issues due to grouping of same values
arrayGroup = repmat([1:10]',3,1);

% Initialize an empty matrix for RGB colors
arrayColorRGB = zeros(length(paramStruct.plotData.arrayColorBase), 3);

% Loop through each hex color and convert to RGB
for i = 1:length(paramStruct.plotData.arrayColorBase)
    arrayColorInd = paramStruct.plotData.arrayColorBase(i);
    arrayColorRGB(i, :) = sscanf(arrayColorInd{1}(2:end), '%2x%2x%2x', [1 3]) / 255;
end

% boxplot1_fig18 = boxplot(res,HW_D2Trio_reshaped_rounded,'colors',arrayColorRGB,'symbol','o');
boxplot1_fig18 = boxplot(real(res),arrayGroup,'colors',arrayColorRGB,'symbol','o');
set(boxplot1_fig18(~isnan(boxplot1_fig18)),'LineWidth',2)
hold on
% boxplot1_fig18 = boxplot(res,HW_D2Trio_reshaped_rounded,'colors','k','symbol','ko');
boxplot2_fig18 = boxplot(real(res),arrayGroup,'colors','k','symbol','ko');
set(boxplot2_fig18(~isnan(boxplot2_fig18)),'LineWidth',1)

set(gca,'XTickLabel',zLabelStrRounded)
grid on
xlabel('Group')
ylabel('Residual')
hold off

%% Set figure properties

%---[Axis properties]---
%----[Font]----
ax1_fig18.FontSize = figProps.Axis.Font.FontSize;
%----[Ticks]----
ax1_fig18.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;

%---[Legend properties]---
%----[Position and Layout]----
%----[Labels]----
xlabel(paramStruct.zParameterHW_D_maxSq.LegendStr,'Interpreter',figProps.Axis.Labels.LabelInterpreter)
ax1_fig18.XLabel.FontSize = figProps.Axis.Labels.LabelFontSize;
ax1_fig18.XLabel.FontWeight = figProps.Axis.Labels.LabelFontWeight;

ylabel('Residual','Interpreter',figProps.Axis.Labels.LabelInterpreter)
ax1_fig18.YLabel.FontSize = figProps.Axis.Labels.LabelFontSize;
ax1_fig18.YLabel.FontWeight = figProps.Axis.Labels.LabelFontWeight;

%% Figure 19 - Generalized master curve 

% Define the equation parameters
a_Kalman = 2.7; % coefficient
fricCoeff = 0.70;

% Define the bounds for x
x_min_Kalman = 1e-3;
x_max_Kalman = 1e8;
x_Kalman = [x_min_Kalman x_max_Kalman];

% Compute the corresponding y values
y_Kalman = a_Kalman .* (x_Kalman).^(3/7);

hold on

% Plot +30% and -30% (as specified by Kalman 2005; 90% of measured points bounded within limits of +_ 30%)
yplus_Kalman = 1.3 .* a_Kalman .* (x_Kalman).^(3/7);

yminus_Kalman = 0.7 .* a_Kalman .* (x_Kalman).^(3/7);

%% Compute Modified Ar**

%Real values only from beta
betaReal = real(beta);

switch modelChoice
    case 'modelFinal'
        Ar_starstar = betaReal(1).*Ar_star.^(betaReal(2)) .* (1./W_DV).^(betaReal(3)) .* (E).^(betaReal(4));
        Re_starstar = 2.7 .* (Ar_starstar).^(3/7);
    case 'modelOne'
        Ar_starstar = betaReal(1).*Ar_star.^(betaReal(2)) .* betaReal(3).*(H_DV./W_DV).^(betaReal(4)) .* betaReal(5).*(1./W_DV).^(betaReal(6)) .* betaReal(7).*(E).^(betaReal(8));
        Re_starstar = 2.7 .* (Ar_starstar).^(3/7);
    case 'modelTwo'
        Ar_starstar = Ar_star.^(betaReal(1)) .* (H_DV./W_DV).^(betaReal(2)) .* (1./W_DV).^(betaReal(3)) .* (E).^(betaReal(4));
        Re_starstar = 2.7 .* (Ar_starstar).^(3/7);
    case 'modelThree'
        Ar_starstar = Ar_star.^(betaReal(1)) .* (1./W_DV).^(betaReal(2)) .* (E).^(betaReal(3));
        Re_starstar = 2.7 .* (Ar_starstar).^(3/7);
    case 'modelFour'
        Ar_starstar = betaReal(1).*Ar_star.^(betaReal(2)) .* (1./W_DV).^(betaReal(3)) .* (E).^(betaReal(4));
        Re_starstar = 2.7 .* (Ar_starstar).^(3/7);
    case 'modelFive'
        Ar_starstar = betaReal(1).*Ar_star.^(betaReal(2)) .* (H_DV./W_DV).^(betaReal(3)) .* (1./W_DV).^(betaReal(4)) .* (E).^(betaReal(5));
        Re_starstar = 2.7 .* (Ar_starstar).^(3/7);
    case 'modelSix'
        Ar_starstar = Ar_star.^(betaReal(1)) .* betaReal(2).*(H_DV./W_DV).^(betaReal(3)) .* (1./W_DV).^(betaReal(4)) .* (E).^(betaReal(5));
        Re_starstar = 2.7 .* (Ar_starstar).^(3/7);
    case 'modelSeven'
        Ar_starstar = Ar_star.^(betaReal(1)) .* (H_DV./W_DV).^(betaReal(2)) .* betaReal(3).*(1./W_DV).^(betaReal(4)) .* (E).^(betaReal(5));
        Re_starstar = 2.7 .* (Ar_starstar).^(3/7);
    case 'modelEight'
        Ar_starstar = Ar_star.^(betaReal(1)) .* (H_DV./W_DV).^(betaReal(2)) .* (1./W_DV).^(betaReal(3)) .* betaReal(4).*(E).^(betaReal(5));
        Re_starstar = 2.7 .* (Ar_starstar).^(3/7);
    case 'modelNine'
        Ar_starstar = betaReal(1).*Ar_star.^(betaReal(2)) .* (H_DV./W_DV).^(betaReal(3)) .* (1./W_DV).^(betaReal(4));
        Re_starstar = 2.7 .* (Ar_starstar).^(3/7);
    case 'modelTen'
        Ar_starstar = Ar_star.^(betaReal(2)) .* (H_DV./W_DV).^(betaReal(3)) .* (1./W_DV).^(betaReal(4)) .* (E).^(betaReal(5));
        Re_starstar = betaReal(1) .* (Ar_starstar).^(betaReal(6));
end

SortedShapes = lineSpec(:,1);
SortedShapes(11:20,1) = lineSpec(1:10,1);
SortedShapes(1:10,1) = lineSpec(11:20,1);

%% Figure 20 - Generalized master curve 
fig20 = figure(20);
ax1_fig20 = gca;

AxisFontSizeMultiplier_fig20 = 6/4.9; %1
AxisLabelFontSizeMultiplier_fig20 = 9/7; %1
LegendLabelFontSizeMultiplier_fig20 = 6/4.9; %1
LegendFontSizeMultiplier_fig20 = 6/4.9; %1

% Custom figure size scale factors (to account for limitations in exported fig size using exportgraphics MATLAB R2024b)
expGraphicsScaleFactor_fig20 = (14/12.3472).*(14/13.8289).*(14/14.0406).*(14/13.8289); %13.8994

%---[Get Figure Properties]---
figProps = getDefaultFigProperties('1.5Column',AxisFontSizeMultiplier_fig20,AxisLabelFontSizeMultiplier_fig20,LegendLabelFontSizeMultiplier_fig20,LegendFontSizeMultiplier_fig20);

%---[Set Figure Size]---
set(fig20,'Units','centimeters','Position',[figProps.Figure.Position.Left figProps.Figure.Position.Bottom figProps.Figure.Position.Width.*expGraphicsScaleFactor_fig20 figProps.Figure.Position.Width*figProps.Figure.Position.HWRatio.*expGraphicsScaleFactor_fig20])
pos = get(fig20,'Position');
set(fig20,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%% Plot Kalman model line
switch modelChoice
    case 'modelFinal'
        % Plot generalized master curve for pickup from a layer of particles in liquid
        % Define the equation parameters
        a_Kalman = 2.7; % coefficient
        % fricCoeff = 0.65;
        fricCoeff = 0.70;

        % Define the bounds for x
        x_min_Kalman = 1e-3;
        x_max_Kalman = 1e8;
        x_Kalman = [x_min_Kalman x_max_Kalman];

        % Compute the corresponding y values
        y_Kalman = a_Kalman .* (x_Kalman).^(3/7);
        KalmanPlotLineObj = plot(x_Kalman, y_Kalman, 'k-', 'LineWidth', 0.5);

        hold on

        % Plot +30% and -30% (as specified by Kalman 2005; 90% of measured points bounded within limits of +_ 30%)
        yplus_Kalman = 1.3 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yplus_Kalman, 'k--', 'LineWidth', 0.25);

        yminus_Kalman = 0.7 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yminus_Kalman, 'k--', 'LineWidth', 0.25);
    case 'modelOne'
        % Plot generalized master curve for pickup from a layer of particles in liquid
        % Define the equation parameters
        a_Kalman = 2.7; % coefficient
        % fricCoeff = 0.65;
        fricCoeff = 0.70;

        % Define the bounds for x
        x_min_Kalman = 1e-3;
        x_max_Kalman = 1e8;
        x_Kalman = [x_min_Kalman x_max_Kalman];

        % Compute the corresponding y values
        y_Kalman = a_Kalman .* (x_Kalman).^(3/7);
        KalmanPlotLineObj = plot(x_Kalman, y_Kalman, 'k-', 'LineWidth', 0.5);

        hold on

        % Plot +30% and -30% (as specified by Kalman 2005; 90% of measured points bounded within limits of +_ 30%)
        yplus_Kalman = 1.3 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yplus_Kalman, 'k--', 'LineWidth', 0.25);

        yminus_Kalman = 0.7 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yminus_Kalman, 'k--', 'LineWidth', 0.25);
    case 'modelTwo'
        % Plot generalized master curve for pickup from a layer of particles in liquid
        % Define the equation parameters
        a_Kalman = 2.7; % coefficient
        fricCoeff = 0.70;

        % Define the bounds for x
        x_min_Kalman = 1e-3;
        x_max_Kalman = 1e8;
        x_Kalman = [x_min_Kalman x_max_Kalman];

        % Compute the corresponding y values
        y_Kalman = a_Kalman .* (x_Kalman).^(3/7);
        KalmanPlotLineObj = plot(x_Kalman, y_Kalman, 'k-', 'LineWidth', 0.5);

        hold on

        % Plot +30% and -30% (as specified by Kalman 2005; 90% of measured points bounded within limits of +_ 30%)
        yplus_Kalman = 1.3 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yplus_Kalman, 'k--', 'LineWidth', 0.25);

        yminus_Kalman = 0.7 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yminus_Kalman, 'k--', 'LineWidth', 0.25);
    case 'modelThree'
        % Plot generalized master curve for pickup from a layer of particles in liquid
        % Define the equation parameters
        a_Kalman = 2.7; % coefficient
        fricCoeff = 0.70;

        % Define the bounds for x
        x_min_Kalman = 1e-3;
        x_max_Kalman = 1e8;
        x_Kalman = [x_min_Kalman x_max_Kalman];

        % Compute the corresponding y values
        y_Kalman = a_Kalman .* (x_Kalman).^(3/7);
        KalmanPlotLineObj = plot(x_Kalman, y_Kalman, 'k-', 'LineWidth', 0.5);

        hold on

        % Plot +30% and -30% (as specified by Kalman 2005; 90% of measured points bounded within limits of +_ 30%)
        yplus_Kalman = 1.3 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yplus_Kalman, 'k--', 'LineWidth', 0.25);

        yminus_Kalman = 0.7 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yminus_Kalman, 'k--', 'LineWidth', 0.25);
    case 'modelFour'
        % Plot generalized master curve for pickup from a layer of particles in liquid
        % Define the equation parameters
        a_Kalman = 2.7; % coefficient
        fricCoeff = 0.70;

        % Define the bounds for x
        x_min_Kalman = 1e-3;
        x_max_Kalman = 1e8;
        x_Kalman = [x_min_Kalman x_max_Kalman];

        % Compute the corresponding y values
        y_Kalman = a_Kalman .* (x_Kalman).^(3/7);
        KalmanPlotLineObj = plot(x_Kalman, y_Kalman, 'k-', 'LineWidth', 0.5);

        hold on

        % Plot +30% and -30% (as specified by Kalman 2005; 90% of measured points bounded within limits of +_ 30%)
        yplus_Kalman = 1.3 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yplus_Kalman, 'k--', 'LineWidth', 0.25);

        yminus_Kalman = 0.7 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yminus_Kalman, 'k--', 'LineWidth', 0.25);
    case 'modelFive'
        % Plot generalized master curve for pickup from a layer of particles in liquid
        % Define the equation parameters
        a_Kalman = 2.7; % coefficient
        fricCoeff = 0.70;

        % Define the bounds for x
        x_min_Kalman = 1e-3;
        x_max_Kalman = 1e8;
        x_Kalman = [x_min_Kalman x_max_Kalman];

        % Compute the corresponding y values
        y_Kalman = a_Kalman .* (x_Kalman).^(3/7);
        KalmanPlotLineObj = plot(x_Kalman, y_Kalman, 'k-', 'LineWidth', 0.5);

        hold on

        % Plot +30% and -30% (as specified by Kalman 2005; 90% of measured points bounded within limits of +_ 30%)
        yplus_Kalman = 1.3 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yplus_Kalman, 'k--', 'LineWidth', 0.25);

        yminus_Kalman = 0.7 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yminus_Kalman, 'k--', 'LineWidth', 0.25);
    case 'modelSix'
        % Plot generalized master curve for pickup from a layer of particles in liquid
        % Define the equation parameters
        a_Kalman = 2.7; % coefficient
        fricCoeff = 0.70;

        % Define the bounds for x
        x_min_Kalman = 1e-3;
        x_max_Kalman = 1e8;
        x_Kalman = [x_min_Kalman x_max_Kalman];

        % Compute the corresponding y values
        y_Kalman = a_Kalman .* (x_Kalman).^(3/7);
        KalmanPlotLineObj = plot(x_Kalman, y_Kalman, 'k-', 'LineWidth', 0.5);

        hold on

        % Plot +30% and -30% (as specified by Kalman 2005; 90% of measured points bounded within limits of +_ 30%)
        yplus_Kalman = 1.3 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yplus_Kalman, 'k--', 'LineWidth', 0.25);

        yminus_Kalman = 0.7 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yminus_Kalman, 'k--', 'LineWidth', 0.25);
    case 'modelSeven'
        % Plot generalized master curve for pickup from a layer of particles in liquid
        % Define the equation parameters
        a_Kalman = 2.7; % coefficient
        fricCoeff = 0.70;

        % Define the bounds for x
        x_min_Kalman = 1e-3;
        x_max_Kalman = 1e8;
        x_Kalman = [x_min_Kalman x_max_Kalman];

        % Compute the corresponding y values
        y_Kalman = a_Kalman .* (x_Kalman).^(3/7);
        KalmanPlotLineObj = plot(x_Kalman, y_Kalman, 'k-', 'LineWidth', 0.5);

        hold on

        % Plot +30% and -30% (as specified by Kalman 2005; 90% of measured points bounded within limits of +_ 30%)
        yplus_Kalman = 1.3 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yplus_Kalman, 'k--', 'LineWidth', 0.25);

        yminus_Kalman = 0.7 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yminus_Kalman, 'k--', 'LineWidth', 0.25);
    case 'modelEight'
        % Plot generalized master curve for pickup from a layer of particles in liquid
        % Define the equation parameters
        a_Kalman = 2.7; % coefficient
        fricCoeff = 0.70;

        % Define the bounds for x
        x_min_Kalman = 1e-3;
        x_max_Kalman = 1e8;
        x_Kalman = [x_min_Kalman x_max_Kalman];

        % Compute the corresponding y values
        y_Kalman = a_Kalman .* (x_Kalman).^(3/7);
        KalmanPlotLineObj = plot(x_Kalman, y_Kalman, 'k-', 'LineWidth', 0.5);

        hold on

        % Plot +30% and -30% (as specified by Kalman 2005; 90% of measured points bounded within limits of +_ 30%)
        yplus_Kalman = 1.3 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yplus_Kalman, 'k--', 'LineWidth', 0.25);

        yminus_Kalman = 0.7 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yminus_Kalman, 'k--', 'LineWidth', 0.25);
    case 'modelNine'
        % Plot generalized master curve for pickup from a layer of particles in liquid
        % Define the equation parameters
        a_Kalman = 2.7; % coefficient
        % fricCoeff = 0.65;
        fricCoeff = 0.70;

        % Define the bounds for x
        x_min_Kalman = 1e-3;
        x_max_Kalman = 1e8;
        x_Kalman = [x_min_Kalman x_max_Kalman];

        % Compute the corresponding y values
        y_Kalman = a_Kalman .* (x_Kalman).^(3/7);
        KalmanPlotLineObj = plot(x_Kalman, y_Kalman, 'k-', 'LineWidth', 0.5);

        hold on

        % Plot +30% and -30% (as specified by Kalman 2005; 90% of measured points bounded within limits of +_ 30%)
        yplus_Kalman = 1.3 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yplus_Kalman, 'k--', 'LineWidth', 0.25);

        yminus_Kalman = 0.7 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yminus_Kalman, 'k--', 'LineWidth', 0.25);
    case 'modelTen'
        % Plot generalized master curve for pickup from a layer of particles in liquid
        % Define the equation parameters
        a_Kalman = 2.7; % coefficient
        fricCoeff = 0.70;

        % Define the bounds for x
        x_min_Kalman = 1e-3;
        x_max_Kalman = 1e8;
        x_Kalman = [x_min_Kalman x_max_Kalman];

        % Compute the corresponding y values
        y_Kalman = a_Kalman .* (x_Kalman).^(3/7);
        KalmanPlotLineObj = plot(x_Kalman, y_Kalman, 'r-', 'LineWidth', 0.5);

        hold on

        % Plot +30% and -30% (as specified by Kalman 2005; 90% of measured points bounded within limits of +_ 30%)
        yplus_Kalman = 1.3 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yplus_Kalman, 'r--', 'LineWidth', 0.25);

        yminus_Kalman = 0.7 .* a_Kalman .* (x_Kalman).^(3/7);
        plot(x_Kalman, yminus_Kalman, 'r--', 'LineWidth', 0.25);

        % Plot generalized master curve for pickup from a layer of particles in liquid
        % Define the equation parameters
        omega_multiplier = beta(1);
        zeta_exponent = beta(6);

        % Define the bounds for x
        x_min_Kalman = 1e-3;
        x_max_Kalman = 1e8;
        x_Kalman = [x_min_Kalman x_max_Kalman];

        % Compute the corresponding y values
        y_Kalman = omega_multiplier .* (x_Kalman).^(zeta_exponent);
        KalmanPlotLineObj = plot(x_Kalman, y_Kalman, 'k-', 'LineWidth', 0.5);

        hold on

        % Plot +30% and -30%
        yplus_Kalman = 1.3 .* omega_multiplier .* (x_Kalman).^(zeta_exponent);
        plot(x_Kalman, yplus_Kalman, 'k--', 'LineWidth', 0.25);

        yminus_Kalman = 0.7 .* omega_multiplier .* (x_Kalman).^(zeta_exponent);
        plot(x_Kalman, yminus_Kalman, 'k--', 'LineWidth', 0.25); 
end

%% Compute Modified Ar**

%Real values only from beta
betaReal = real(beta);

switch modelChoice
    case 'modelFinal'
        Ar_starstar = betaReal(1).*Ar_star.^(betaReal(2)) .* (1./W_DV).^(betaReal(3)) .* (E).^(betaReal(4));
    case 'modelOne'
        Ar_starstar = betaReal(1).*Ar_star.^(betaReal(2)) .* betaReal(3).*(H_DV./W_DV).^(betaReal(4)) .* betaReal(5).*(1./W_DV).^(betaReal(6)) .* betaReal(7).*(E).^(betaReal(8));
    case 'modelTwo'
        Ar_starstar = Ar_star.^(betaReal(1)) .* (H_DV./W_DV).^(betaReal(2)) .* (1./W_DV).^(betaReal(3)) .* (E).^(betaReal(4));
    case 'modelThree'
        Ar_starstar = Ar_star.^(betaReal(1)) .* (1./W_DV).^(betaReal(2)) .* (E).^(betaReal(3));
    case 'modelFour'
        Ar_starstar = betaReal(1).*Ar_star.^(betaReal(2)) .* (1./W_DV).^(betaReal(3)) .* (E).^(betaReal(4));
    case 'modelFive'
        Ar_starstar = betaReal(1).*Ar_star.^(betaReal(2)) .* (H_DV./W_DV).^(betaReal(3)) .* (1./W_DV).^(betaReal(4)) .* (E).^(betaReal(5));
    case 'modelSix'
        Ar_starstar = Ar_star.^(betaReal(1)) .* betaReal(2).*(H_DV./W_DV).^(betaReal(3)) .* (1./W_DV).^(betaReal(4)) .* (E).^(betaReal(5));
    case 'modelSeven'
        Ar_starstar = Ar_star.^(betaReal(1)) .* (H_DV./W_DV).^(betaReal(2)) .* betaReal(3).*(1./W_DV).^(betaReal(4)) .* (E).^(betaReal(5));
    case 'modelEight'
        Ar_starstar = Ar_star.^(betaReal(1)) .* (H_DV./W_DV).^(betaReal(2)) .* (1./W_DV).^(betaReal(3)) .* betaReal(4).*(E).^(betaReal(5));
    case 'modelNine'
        Ar_starstar = betaReal(1).*Ar_star.^(betaReal(2)) .* (H_DV./W_DV).^(betaReal(3)) .* (1./W_DV).^(betaReal(4));
    case 'modelTen'
        Ar_starstar = Ar_star.^(betaReal(2)) .* (H_DV./W_DV).^(betaReal(3)) .* (1./W_DV).^(betaReal(4)) .* (E).^(betaReal(5));
end

SortedShapes = lineSpec(:,1);
SortedShapes(11:20,1) = lineSpec(1:10,1);
SortedShapes(1:10,1) = lineSpec(11:20,1);

%% Drawing modified points as scatter markers
hold on
for mm = 1:length(Re_starstar)
    pstarstar(mm) = scatter(Ar_starstar(mm),Re_starstar(mm)+real(res(mm)),figProps.Scatter.Markers.SizeData);
    % pstarstar(mm) = scatter(Ar_starstar(mm),Re_starstar(mm)+real(res(mm)),figProps.Scatter.Markers.SizeData);
    set(pstarstar(mm),{'Marker','MarkerFaceColor','MarkerFaceAlpha'},{SortedShapes(mm,1),lineSpec(mm,2),figProps.Scatter.Markers.MarkerFaceAlpha})
end
set(pstarstar,'MarkerEdgeColor',figProps.Scatter.Markers.MarkerEdgeColor)

%%
uistack(pOriginal,'up',30)

%% Set figure properties
%---[Line properties]---
%----[Line]----
%----[Markers]----
set(pstarstar,'MarkerEdgeColor',figProps.Line.Markers.MarkerEdgeColor)

%---[Axis properties]---
%----[Font]----
ax1_fig20.FontSize = figProps.Axis.Font.FontSize;
%----[Ticks]----
ax1_fig20.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;
%----[Rulers]----
switch modelChoice
    case 'modelFinal'
        xlim([1 3000])
        xticks(logspace(0,4,5));
        ylim([1 1000])
        yticks(logspace(0,4,5));
    case {'modelOne','modelTwo','modelThree','modelFour','modelFive','modelSix','modelSeven','modelEight','modelNine', 'modelTen'}
        xlim([1 3000])
        xticks(logspace(0,4,5));
        ylim([1 3000])
        yticks(logspace(0,4,5));
end

ax1_fig20.XAxis.Exponent = 0; %Scientific notation
xtickformat('%.0f') %Scientific notation
ax1_fig20.YAxis.Exponent = 0; %Scientific notation
ytickformat('%.0f') %Scientific notation
set(gca,'XScale','log');
set(gca,'YScale','log');
%----[Grids]----
grid on
ax1_fig20.GridLineStyle = figProps.Axis.Grids.GridLineStyle;
ax1_fig20.GridLineWidth = figProps.Axis.Grids.GridLineWidth;
ax1_fig20.GridColor = figProps.Axis.Grids.GridColor;
ax1_fig20.MinorGridLineStyle = figProps.Axis.Grids.MinorGridLineStyle;
ax1_fig20.MinorGridLineWidth = figProps.Axis.Grids.MinorGridLineWidth;
ax1_fig20.MinorGridColor = figProps.Axis.Grids.MinorGridColor;
%----[Labels]----
switch modelChoice
    case 'modelFinal'
        xlabel("Remodified Archimedes number, $\mathrm{Ar}_{d_\mathrm{V}}^{**}$",'Interpreter',figProps.Axis.Labels.LabelInterpreter)
    case 'modelOne'
        xlabel("Remodified Archimedes number, $\mathrm{Ar}_{d_\mathrm{V}}^{[LO-1]}$",'Interpreter',figProps.Axis.Labels.LabelInterpreter)
    case 'modelTwo'
        xlabel("Remodified Archimedes number, $\mathrm{Ar}_{d_\mathrm{V}}^{[LO-2]}$",'Interpreter',figProps.Axis.Labels.LabelInterpreter)
    case 'modelThree'
        xlabel("Remodified Archimedes number, $\mathrm{Ar}_{d_\mathrm{V}}^{[LO-3]}$",'Interpreter',figProps.Axis.Labels.LabelInterpreter)
    case 'modelFour'
        xlabel("Remodified Archimedes number, $\mathrm{Ar}_{d_\mathrm{V}}^{[LO-4]}$",'Interpreter',figProps.Axis.Labels.LabelInterpreter)
    case 'modelFive'
        xlabel("Remodified Archimedes number, $\mathrm{Ar}_{d_\mathrm{V}}^{[LO-5]}$",'Interpreter',figProps.Axis.Labels.LabelInterpreter)
    case 'modelSix'
        xlabel("Remodified Archimedes number, $\mathrm{Ar}_{d_\mathrm{V}}^{[LO-6]}$",'Interpreter',figProps.Axis.Labels.LabelInterpreter)
    case 'modelSeven'
        xlabel("Remodified Archimedes number, $\mathrm{Ar}_{d_\mathrm{V}}^{[LO-7]}$",'Interpreter',figProps.Axis.Labels.LabelInterpreter)
    case 'modelEight'
        xlabel("Remodified Archimedes number, $\mathrm{Ar}_{d_\mathrm{V}}^{[LO-8]}$",'Interpreter',figProps.Axis.Labels.LabelInterpreter)
    case 'modelNine'
        xlabel("Remodified Archimedes number, $\mathrm{Ar}_{d_\mathrm{V}}^{[LO-9]}$",'Interpreter',figProps.Axis.Labels.LabelInterpreter)
    case 'modelTen'
        xlabel("Remodified Archimedes number, $\mathrm{Ar}_{d_\mathrm{V}}^{[LO-10]}$",'Interpreter',figProps.Axis.Labels.LabelInterpreter)
end

ax1_fig20.XLabel.FontSize = figProps.Axis.Labels.LabelFontSize;
ax1_fig20.XLabel.FontWeight = figProps.Axis.Labels.LabelFontWeight;

ylabel({'Modified critical';'particle Reynolds number, $\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}^*$'},'Interpreter',figProps.Axis.Labels.LabelInterpreter)

ax1_fig20.YLabel.FontSize = figProps.Axis.Labels.LabelFontSize;
ax1_fig20.YLabel.FontWeight = figProps.Axis.Labels.LabelFontWeight;
%----[Multiple Plots]----
%----[Color and Transparency]----
%----[Box Styling]----
set(ax1_fig20, 'Box', 'Off')

%% Add text annotation parallel to the KalmanPlotLineObj line
switch modelChoice
    case 'modelFinal'
        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_main = 200;
        y_annot_main = 70;
        %---[Add the text annotation]---
        str_main = '${\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}}^* = 2.7 (\overbrace{\alpha \left( \mathrm{Ar}_{d_\mathrm{V}}^{*}\right)^{\beta} \left( {d_\mathrm{V}/W}\right)^{\theta} \left( {\epsilon_\mathrm{r}}\right)^{\psi}}^{\displaystyle{\mathrm{Ar}_{d_\mathrm{V}}^{**}}} )^{3/7}$';
    case 'modelOne'
        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_main = 200;
        y_annot_main = 70;
        str_main = '${\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}}^* = 2.7 (\overbrace{\alpha_{1} \left( \mathrm{Ar}_{d_\mathrm{V}}^{*}\right)^{\beta} \alpha_{2} \left( {H/W}\right)^{\xi} \alpha_{3} \left( {d_\mathrm{V}/W}\right)^{\theta} \alpha_{4} \left( {\epsilon_\mathrm{r}}\right)^{\psi}}^{\displaystyle{\mathrm{Ar}_{d_\mathrm{V}}^{[LO-1]}}} )^{3/7}$';
    case 'modelTwo'
        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_main = 200;
        y_annot_main = 70;
        str_main = '${\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}}^* = 2.7 (\overbrace{\left( \mathrm{Ar}_{d_\mathrm{V}}^{*}\right)^{\beta} \left( {H/W}\right)^{\xi} \left( {d_\mathrm{V}/W}\right)^{\theta} \left( {\epsilon_\mathrm{r}}\right)^{\psi}}^{\displaystyle{\mathrm{Ar}_{d_\mathrm{V}}^{[LO-2]}}} )^{3/7}$';
    case 'modelThree'
        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_main = 500;
        y_annot_main = 90;
        %---[Add the text annotation]---
        str_main = '${\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}}^* = 2.7 (\overbrace{\left( \mathrm{Ar}_{d_\mathrm{V}}^{*}\right)^{\beta} \left( {d_\mathrm{V}/W}\right)^{\theta} \left( {\epsilon_\mathrm{r}}\right)^{\psi}}^{\displaystyle{\mathrm{Ar}_{d_\mathrm{V}}^{[LO-3]}}} )^{3/7}$';
    case 'modelFour'
        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_main = 500;
        y_annot_main = 90;
        %---[Add the text annotation]---
        str_main = '${\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}}^* = 2.7 (\overbrace{\alpha \left( \mathrm{Ar}_{d_\mathrm{V}}^{*}\right)^{\beta} \left( {d_\mathrm{V}/W}\right)^{\theta} \left( {\epsilon_\mathrm{r}}\right)^{\psi}}^{\displaystyle{\mathrm{Ar}_{d_\mathrm{V}}^{[LO-4]}}} )^{3/7}$';
    case 'modelFive'
        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_main = 200;
        y_annot_main = 70;
        %---[Add the text annotation]---
        str_main = '${\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}}^* = 2.7 (\overbrace{\alpha_{1} \left( \mathrm{Ar}_{d_\mathrm{V}}^{*}\right)^{\beta} \left( {H/W}\right)^{\xi} \left( {d_\mathrm{V}/W}\right)^{\theta} \left( {\epsilon_\mathrm{r}}\right)^{\psi}}^{\displaystyle{\mathrm{Ar}_{d_\mathrm{V}}^{[LO-5]}}} )^{3/7}$';
    case 'modelSix'
        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_main = 200;
        y_annot_main = 70;
        %---[Add the text annotation]---
        str_main = '${\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}}^* = 2.7 (\overbrace{\left( \mathrm{Ar}_{d_\mathrm{V}}^{*}\right)^{\beta} \alpha_{2} \left( {H/W}\right)^{\xi} \left( {d_\mathrm{V}/W}\right)^{\theta} \left( {\epsilon_\mathrm{r}}\right)^{\psi}}^{\displaystyle{\mathrm{Ar}_{d_\mathrm{V}}^{[LO-6]}}} )^{3/7}$';
    case 'modelSeven'
        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_main = 200;
        y_annot_main = 70;
        %---[Add the text annotation]---
        str_main = '${\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}}^* = 2.7 (\overbrace{\left( \mathrm{Ar}_{d_\mathrm{V}}^{*}\right)^{\beta} \left( {H/W}\right)^{\xi} \alpha_{3} \left( {d_\mathrm{V}/W}\right)^{\theta} \left( {\epsilon_\mathrm{r}}\right)^{\psi}}^{\displaystyle{\mathrm{Ar}_{d_\mathrm{V}}^{[LO-7]}}} )^{3/7}$';
    case 'modelEight'
        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_main = 200;
        y_annot_main = 70;
        %---[Add the text annotation]---
        str_main = '${\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}}^* = 2.7 (\overbrace{\left( \mathrm{Ar}_{d_\mathrm{V}}^{*}\right)^{\beta} \left( {H/W}\right)^{\xi} \left( {d_\mathrm{V}/W}\right)^{\theta} \alpha_{4} \left( {\epsilon_\mathrm{r}}\right)^{\psi}}^{\displaystyle{\mathrm{Ar}_{d_\mathrm{V}}^{[LO-8]}}} )^{3/7}$';
    case 'modelNine'
        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_main = 200;
        y_annot_main = 70;
        %---[Add the text annotation]---
        str_main = '${\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}}^* = 2.7 (\overbrace{\alpha_{1} \left( \mathrm{Ar}_{d_\mathrm{V}}^{*}\right)^{\beta} \left( {H/W}\right)^{\xi} \left( {d_\mathrm{V}/W}\right)^{\theta} )^{\psi}}^{\displaystyle{\mathrm{Ar}_{d_\mathrm{V}}^{[LO-9]}}} )^{3/7}$';
    case 'modelTen'
        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_main1 = 70;
        y_annot_main1 = 30;
        str_main1 = '${\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}}^* = 2.7 \left( \mathrm{Ar}_{d_\mathrm{V}}^{*}\right)^{3/7}$';

        x_annot_main2 = 500;
        y_annot_main2 = 4;        
        str_main2 = '${\mathrm{\check{Re}}_{\mathrm{p},d_\mathrm{V}}}^* = \omega (\overbrace{\left( \mathrm{Ar}_{d_\mathrm{V}}^{*}\right)^{\beta} \left( {H/W}\right)^{\xi} \left( {d_\mathrm{V}/W}\right)^{\theta} \left( {\epsilon_\mathrm{r}}\right)^{\psi}}^{\displaystyle{\mathrm{Ar}_{d_\mathrm{V}}^{[LO-10]}}} )^{\zeta}$';
end

switch modelChoice
    case {'modelFinal','modelOne','modelTwo','modelThree','modelFour','modelFive','modelSix','modelSeven','modelEight','modelNine'}
        eqnTextAnnot_fig20 = B2KTextAlignedToLineLogLogPlot(gca,x_annot_main,y_annot_main,str_main,3/7,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter','latex','Color','k');

        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_plus = 1.5;
        y_annot_plus = 4.5;
        %---[Add the text annotation]---
        str_plus = '$\mathrm{+}$30\%';
        plusTextAnnot = B2KTextAlignedToLineLogLogPlot(gca,x_annot_plus,y_annot_plus,str_plus,3/7,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter','latex');

        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_minus = 1.5;
        y_annot_minus = 1.4;
        %---[Add the text annotation]---
        str_minus = '$\mathrm{-}$30\%';
        minusTextAnnot = B2KTextAlignedToLineLogLogPlot(gca,x_annot_minus,y_annot_minus,str_minus,3/7,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter','latex');
    case 'modelTen'
        eqnTextAnnot1_fig20 = B2KTextAlignedToLineLogLogPlot(gca,x_annot_main1,y_annot_main1,str_main1,3/7,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter','latex','Color','r');

        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_plus = 1.5;
        y_annot_plus = 4.5;
        %---[Add the text annotation]---
        str_plus = '$\mathrm{+}$30\%';
        plusTextAnnot = B2KTextAlignedToLineLogLogPlot(gca,x_annot_plus,y_annot_plus,str_plus,3/7,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter','latex','Color','r');

        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_minus = 1.5;
        y_annot_minus = 1.4;
        %---[Add the text annotation]---
        str_minus = '$\mathrm{-}$30\%';
        minusTextAnnot = B2KTextAlignedToLineLogLogPlot(gca,x_annot_minus,y_annot_minus,str_minus,3/7,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter','latex','Color','r');       

        eqnTextAnnot2_fig20 = B2KTextAlignedToLineLogLogPlot(gca,x_annot_main2,y_annot_main2,str_main2,beta(6),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter','latex','Color','k');

        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_plus = 30;
        y_annot_plus = 4.25;
        %---[Add the text annotation]---
        str_plus = '$\mathrm{+}$30\%';
        plusTextAnnot = B2KTextAlignedToLineLogLogPlot(gca,x_annot_plus,y_annot_plus,str_plus,beta(6),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter','latex');

        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_minus = 30;
        y_annot_minus = 1.25;
        %---[Add the text annotation]---
        str_minus = '$\mathrm{-}$30\%';
        minusTextAnnot = B2KTextAlignedToLineLogLogPlot(gca,x_annot_minus,y_annot_minus,str_minus,beta(6),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter','latex');
end

%% Text annotation for assumed friction coefficient
switch modelChoice
    case {'modelFinal','modelOne','modelTwo','modelThree','modelFour','modelFive','modelSix','modelSeven','modelEight','modelNine'}
        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_fric = 32.5; %120
        y_annot_fric = 2.4;
        %---[Add the text annotation]---
        str_fric = sprintf('Assumed $\\mu_{s} = %.2f$', fricCoeff);
        text(x_annot_fric,y_annot_fric,str_fric,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter','latex');
    case 'modelTen'
        %---[Choose a point (x, y) on the line to place the annotation]---
        x_annot_fric = 3.25;
        y_annot_fric = 30;
        %---[Add the text annotation]---
        str_fric = sprintf('Assumed $\\mu_{s} = %.2f$', fricCoeff);
        text(x_annot_fric,y_annot_fric,str_fric,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter','latex');
end

%% Text annotation for regressed parameters
% Pick symbols
switch modelChoice
    case 'modelFinal'
        % Anchor & vertical spacing size
        x_annot  = 600;
        y_annot  = 15;
        deltaVal = 0.175;  % spacing in decades (if log) or units (if linear)
        syms = {'\alpha','\beta','\theta','\psi'};
    case 'modelOne'
        % Anchor & vertical spacing size
        x_annot  = 600;
        y_annot  = 15;
        deltaVal = 0.15;  % spacing in decades (if log) or units (if linear)
        syms = {'\alpha_1','\beta','\alpha_2','\xi','\alpha_3','\theta','\alpha_4','\psi'};
    case 'modelTwo'
        % Anchor & vertical spacing size
        x_annot  = 600;
        y_annot  = 15;
        deltaVal = 0.15;  % spacing in decades (if log) or units (if linear)
        syms = {'\beta','\xi','\theta','\psi'};
    case 'modelThree'
        % Anchor & vertical spacing size
        x_annot  = 600;
        y_annot  = 15;
        deltaVal = 0.15;  % spacing in decades (if log) or units (if linear)
        syms = {'\beta','\theta','\psi'};
    case 'modelFour'
        % Anchor & vertical spacing size
        x_annot  = 600;
        y_annot  = 15;
        deltaVal = 0.15;  % spacing in decades (if log) or units (if linear)
        syms = {'\alpha','\beta','\theta','\psi'};
    case 'modelFive'
        % Anchor & vertical spacing size
        x_annot  = 600;
        y_annot  = 15;
        deltaVal = 0.15;  % spacing in decades (if log) or units (if linear)
        syms = {'\alpha_1','\beta','\xi','\theta','\psi'};
    case 'modelSix'
        % Anchor & vertical spacing size
        x_annot  = 600;
        y_annot  = 15;
        deltaVal = 0.15;  % spacing in decades (if log) or units (if linear)
        syms = {'\beta','\alpha_2','\xi','\theta','\psi'};
    case 'modelSeven'
        % Anchor & vertical spacing size
        x_annot  = 600;
        y_annot  = 15;
        deltaVal = 0.15;  % spacing in decades (if log) or units (if linear)
        syms = {'\beta','\xi','\alpha_3','\theta','\psi'};
    case 'modelEight'
        % Anchor & vertical spacing size
        x_annot  = 600;
        y_annot  = 15;
        deltaVal = 0.15;  % spacing in decades (if log) or units (if linear)
        syms = {'\beta','\xi','\theta','\alpha_4','\psi'};
    case 'modelNine'
        % Anchor & vertical spacing size
        x_annot  = 600;
        y_annot  = 15;
        deltaVal = 0.15;  % spacing in decades (if log) or units (if linear)
        syms = {'\alpha_1','\beta','\xi','\theta'};
    case 'modelTen'
        % Anchor & vertical spacing size
        x_annot  = 600;
        y_annot  = 600;
        deltaVal = 0.15;  % spacing in decades (if log) or units (if linear)
        syms = {'\omega','\beta','\xi','\theta','\psi','\zeta'};
    otherwise
        error('Unknown modelChoice: %s', modelChoice);
end

% No need to specify 'YSpacing' unless overriding auto-detection
h = B2KAnnotateModelParams( ...
  gca, syms, beta, SE, ...
  x_annot, y_annot, deltaVal, ...
  'Interpreter', 'latex', ...
  'FontSize',    figProps.Legend.Labels.Title.FontSize, ...
  'AlignEqual',  true, ...
  'AlignPM',     true, ... % <- auto-detects 'log' vs 'linear'
  'Color', 'k');  

%% Zero legend [Title: Flat-plate, experimental]
yShiftAnnotation_fig20b = -0.0225;

%% First legend [HW/(D_{V})^2, set progressive color order]
%---[Legends strings - Dynamic]---
zLabelStrRounded = B2KFormatValuesToStrOneSigFigMaxUnc(paramStruct.zParameterHW_D_VSq.ValueSorted,paramStruct.zParameterHW_D_VSq.UncertaintySorted);

% Create dummy plot to generate first legend since default legend for scatter cannot modify legend marker size
for m = 1:refSize
    pDummy(m) = plot(NaN,NaN,'MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none');
    set(pDummy(m),{'Marker','MarkerFaceColor'},{lineSpec(m,1),[lineSpec(m,2)]})
    hold on
end
set(pDummy,'MarkerEdgeColor',figProps.Line.Markers.MarkerEdgeColor)

% Plot first legend
hold on
lh1_fig20b = legend(gca,pDummy(1+refSizeEach:refSizeEach*2),zLabelStrRounded,'FontSize',figProps.Legend.Fonts.FontSize,'NumColumns',1,'Interpreter',figProps.Legend.Interpreter); %circle marker
lh1_fig20b.Units = 'normalized';
lh1_fig20b.Location = 'northwest';
lh1_fig20b_title = title(lh1_fig20b,paramStruct.zParameterHW_D_VSq.LegendStr,'FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
lh1_fig20b.EdgeColor = figProps.Legend.ColorAndStyling.EdgeColor;

% Adjust horizontal spacing for legend lh1
lh1_fig20b.ItemTokenSize(1) = figProps.Legend.Fonts.ItemTokenSize;

% Set position of first legend
lh1_OrigPos_fig20b = get(lh1_fig20b,'Position');
set(lh1_fig20b,'Position',[lh1_OrigPos_fig20b(1), lh1_OrigPos_fig20b(2)-yShiftAnnotation_fig20b, lh1_OrigPos_fig20b(3), lh1_OrigPos_fig20b(4)]);
lh1titleTop_fig20b = annotation('line', [lh1_OrigPos_fig20b(1), lh1_OrigPos_fig20b(1) + lh1_OrigPos_fig20b(3)], [lh1_OrigPos_fig20b(2)-yShiftAnnotation_fig20b + lh1_OrigPos_fig20b(4), lh1_OrigPos_fig20b(2)-yShiftAnnotation_fig20b + lh1_OrigPos_fig20b(4)], 'Color', 'black');
lh1titleBottom_fig20b = annotation('line', [lh1_OrigPos_fig20b(1), lh1_OrigPos_fig20b(1) + lh1_OrigPos_fig20b(3)], [lh1_OrigPos_fig20b(2)-yShiftAnnotation_fig20b,lh1_OrigPos_fig20b(2)-yShiftAnnotation_fig20b], 'Color', 'black');
lh1titleBottom2_fig20b = annotation('line', [lh1_OrigPos_fig20b(1), lh1_OrigPos_fig20b(1) + lh1_OrigPos_fig20b(3)], [lh1_OrigPos_fig20b(2)+yShiftAnnotation_fig20b*1.2 + lh1_OrigPos_fig20b(4), lh1_OrigPos_fig20b(2)+yShiftAnnotation_fig20b*1.2 + lh1_OrigPos_fig20b(4)], 'Color', 'black');

%% Second legend [H/D_{V}, dependent color order] [TESTING]

% Copy the axes and plot the second legend
ax2_fig20 = copyobj(ax1_fig20,gcf); %Copies ah1 object (but duplicates data points)
delete(get(ax2_fig20,'Children')) %Deletes duplicate data points

% Remove second axis visibility 
set(ax2_fig20, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

%---Legends strings
zLabelStrRounded_lh2_fig20b = B2KFormatValuesToStrOneSigFigMaxUnc(paramStruct.hD_VParameter.ValueSorted,paramStruct.hD_VParameter.UncertaintySorted);

% Plot second legend
hold on
lh2dump1_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig20);
lh2dump2_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig20);
lh2dump3_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig20);
lh2dump4_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig20);
lh2dump5_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig20);
lh2dump6_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig20);
lh2dump7_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig20);
lh2dump8_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig20);
lh2dump9_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig20);
lh2dump10_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax2_fig20);
hold off

lh2_fig20b = legend(ax2_fig20,[lh2dump1_fig20b lh2dump2_fig20b lh2dump3_fig20b lh2dump4_fig20b lh2dump5_fig20b lh2dump6_fig20b lh2dump7_fig20b lh2dump8_fig20b lh2dump9_fig20b lh2dump10_fig20b],zLabelStrRounded_lh2_fig20b,'FontSize',figProps.Legend.Fonts.FontSize,'NumColumns',1,'Interpreter',figProps.Legend.Interpreter);
lh2_fig20b.Units = 'normalized';
lh2_fig20b.Location = 'northwest';
lh2title_fig20b = title(lh2_fig20b,"$H/d_\mathrm{V}$",'FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
lh2_fig20b.EdgeColor = figProps.Legend.ColorAndStyling.EdgeColor;

% Set background color of second legend to white (instead of grey by default)
set(lh2_fig20b,'Color',get(ax1_fig2,'Color'));

% Adjust horizontal spacing for legend lh2
lh2_fig20b.ItemTokenSize(1) = 0; %figProps.Legend.Fonts.ItemTokenSize

% Set position relative to first legend
lh1_OrigPos_fig20b = get(lh1_fig20b,'Position');
lh2_OrigPos_fig20b = get(lh2_fig20b,'Position');
set(lh2_fig20b,'Position',[lh2_OrigPos_fig20b(1)+lh1_OrigPos_fig20b(3), lh2_OrigPos_fig20b(2)-yShiftAnnotation_fig20b, lh2_OrigPos_fig20b(3), lh2_OrigPos_fig20b(4)]);
lh2titleTop_fig20b = annotation('line', [lh2_OrigPos_fig20b(1)+lh1_OrigPos_fig20b(3), lh2_OrigPos_fig20b(1)+lh1_OrigPos_fig20b(3) + lh2_OrigPos_fig20b(3)], [lh2_OrigPos_fig20b(2)-yShiftAnnotation_fig20b + lh2_OrigPos_fig20b(4), lh2_OrigPos_fig20b(2)-yShiftAnnotation_fig20b + lh2_OrigPos_fig20b(4)], 'Color', 'black');
lh2titleBottom_fig20b = annotation('line', [lh2_OrigPos_fig20b(1)+lh1_OrigPos_fig20b(3), lh2_OrigPos_fig20b(1)+lh1_OrigPos_fig20b(3) + lh2_OrigPos_fig20b(3)], [lh2_OrigPos_fig20b(2)-yShiftAnnotation_fig20b,lh2_OrigPos_fig20b(2)-yShiftAnnotation_fig20b], 'Color', 'black');
lh2titleBottom2_fig20b = annotation('line', [lh2_OrigPos_fig20b(1)+lh1_OrigPos_fig20b(3), lh2_OrigPos_fig20b(1)+lh1_OrigPos_fig20b(3) + lh2_OrigPos_fig20b(3)], [lh2_OrigPos_fig20b(2)+yShiftAnnotation_fig20b*1.2 + lh2_OrigPos_fig20b(4), lh2_OrigPos_fig20b(2)+yShiftAnnotation_fig20b*1.2 + lh2_OrigPos_fig20b(4)], 'Color', 'black');

%% Third legend [W/D_{V}] [TESTING]

% Copy the axes and plot the third legend
ax3_fig20 = copyobj(ax2_fig20,gcf); %Copies ah1 object (but duplicates data points)
delete(get(ax3_fig20,'Children')) %Deletes duplicate data points

% Remove third axis visibility 
set(ax3_fig20, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

%---Legends strings
zLabelStrRounded_lh3_fig20b = B2KFormatValuesToStrOneSigFigMaxUnc(paramStruct.wD_VParameter.ValueSorted,paramStruct.wD_VParameter.UncertaintySorted);

% Plot third legend
hold on
lh3dump1_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig20);
lh3dump2_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig20);
lh3dump3_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig20);
lh3dump4_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig20);
lh3dump5_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig20);
lh3dump6_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig20);
lh3dump7_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig20);
lh3dump8_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig20);
lh3dump9_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig20);
lh3dump10_fig20b = plot(NaN,'Marker','none','LineStyle','none','Parent',ax3_fig20);
hold off

lh3_fig20b = legend(ax3_fig20,[lh3dump1_fig20b lh3dump2_fig20b lh3dump3_fig20b lh3dump4_fig20b lh3dump5_fig20b lh3dump6_fig20b lh3dump7_fig20b lh3dump8_fig20b lh3dump9_fig20b lh3dump10_fig20b],zLabelStrRounded_lh3_fig20b,'FontSize',figProps.Legend.Fonts.FontSize,'NumColumns',1,'Interpreter','latex');
lh3_fig20b.Units = 'normalized';
lh3_fig20b.Location = 'northwest';
lh3title_fig20b = title(lh3_fig20b,"$W/d_\mathrm{V}$",'FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
lh3_fig20b.EdgeColor = figProps.Legend.ColorAndStyling.EdgeColor;

% Set background color of third legend to white (instead of grey by default)
set(lh3_fig20b,'Color',get(tileax1_fig16,'Color'));

% Adjust horizontal spacing for legend lh3
lh3_fig20b.ItemTokenSize(1) = 0; %figProps.Legend.Fonts.ItemTokenSize

% Set position relative to second legend
lh2_OrigPos_fig20b = get(lh2_fig20b,'Position');
lh3_OrigPos_fig20b = get(lh3_fig20b,'Position');
set(lh3_fig20b,'Position',[lh3_OrigPos_fig20b(1)+lh2_OrigPos_fig20b(3)+lh1_OrigPos_fig20b(3), lh3_OrigPos_fig20b(2)-yShiftAnnotation_fig20b, lh3_OrigPos_fig20b(3), lh3_OrigPos_fig20b(4)]);
lh3titleTop_fig20b = annotation('line', [lh3_OrigPos_fig20b(1)+lh2_OrigPos_fig20b(3)+lh1_OrigPos_fig20b(3), lh3_OrigPos_fig20b(1)+lh2_OrigPos_fig20b(3)+lh1_OrigPos_fig20b(3) + lh3_OrigPos_fig20b(3)], [lh3_OrigPos_fig20b(2)-yShiftAnnotation_fig20b + lh3_OrigPos_fig20b(4), lh3_OrigPos_fig20b(2)-yShiftAnnotation_fig20b + lh3_OrigPos_fig20b(4)], 'Color', 'black');
lh3titleBottom_fig20b = annotation('line', [lh3_OrigPos_fig20b(1)+lh2_OrigPos_fig20b(3)+lh1_OrigPos_fig20b(3), lh3_OrigPos_fig20b(1)+lh2_OrigPos_fig20b(3)+lh1_OrigPos_fig20b(3) + lh3_OrigPos_fig20b(3)], [lh3_OrigPos_fig20b(2)-yShiftAnnotation_fig20b,lh3_OrigPos_fig20b(2)-yShiftAnnotation_fig20b], 'Color', 'black');
lh2titleBottom2_fig20b = annotation('line', [lh3_OrigPos_fig20b(1)+lh2_OrigPos_fig20b(3)+lh1_OrigPos_fig20b(3), lh3_OrigPos_fig20b(1)+lh2_OrigPos_fig20b(3)+lh1_OrigPos_fig20b(3) + lh3_OrigPos_fig20b(3)], [lh3_OrigPos_fig20b(2)+yShiftAnnotation_fig20b*1.2 + lh3_OrigPos_fig20b(4), lh3_OrigPos_fig20b(2)+yShiftAnnotation_fig20b*1.2 + lh3_OrigPos_fig20b(4)], 'Color', 'black');

%% Fourth legend [Solvent]

% Create array for markers and labels for fourth legend
lh4_fig20b_arrayMarkerdump = ["o" "square" "diamond"];
lh4_fig20b_typeSolvent = cellstr(["Isopropanol","Water","Methanol"]);

% Copy the axes and plot the fourth legend
ax4_fig20 = copyobj(ax3_fig20,gcf); %Copies ah1 object (but duplicates data points)
delete(get(ax4_fig20,'Children')) %Deletes duplicate data points

% Remove fourth axis visibility 
set(ax4_fig20, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

% Plot "dummy" data for fourth legend; tried using single variable with loop but does not work
hold on
lh4dump1_fig20b = plot(NaN,'Marker',lh4_fig20b_arrayMarkerdump(1),'DisplayName',lh4_fig20b_typeSolvent{1},'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none','Parent',ax4_fig20);
lh4dump2_fig20b = plot(NaN,'Marker',lh4_fig20b_arrayMarkerdump(2),'DisplayName',lh4_fig20b_typeSolvent{2},'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none','Parent',ax4_fig20);
lh4dump3_fig20b = plot(NaN,'Marker',lh4_fig20b_arrayMarkerdump(3),'DisplayName',lh4_fig20b_typeSolvent{3},'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',sqrt(figProps.Scatter.Markers.SizeData),'LineStyle','none','Parent',ax4_fig20);
hold off

lh4_fig20b = legend(ax4_fig20,[lh4dump1_fig20b lh4dump2_fig20b lh4dump3_fig20b],lh4_fig20b_typeSolvent,'FontSize',figProps.Legend.Fonts.FontSize,'Interpreter',figProps.Legend.Interpreter);
lh4_fig20b.Units = 'normalized';
lh4_fig20b.Location = 'northwest';
lh4title_fig20b = title(lh4_fig20b,'Solvent','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Text.Interpreter);
lh4_fig20b.EdgeColor = [1 1 1];

% Set background color of fourth legend to white (instead of grey by default)
set(lh4_fig20b,'Color',get(ax1_fig2,'Color'));

% Adjust horizontal spacing for legend lh4
lh4_fig20b.ItemTokenSize(1) = figProps.Legend.Fonts.ItemTokenSize;

% Set position relative to third legend
lh3_OrigPos_fig20b = get(lh3_fig20b,'Position');
lh4_OrigPos_fig20b = get(lh4_fig20b,'Position');
lh4_addLeftSpace_fig20b = 0.025;
set(lh4_fig20b,'Position',[lh4_addLeftSpace_fig20b+lh4_OrigPos_fig20b(1)+lh3_OrigPos_fig20b(3)+lh2_OrigPos_fig20b(3)+lh1_OrigPos_fig20b(3), lh4_OrigPos_fig20b(2)-yShiftAnnotation_fig20b, lh4_OrigPos_fig20b(3), lh4_OrigPos_fig20b(4)]);
lh4titleTop_fig20b = annotation('line', [lh4_addLeftSpace_fig20b+lh4_OrigPos_fig20b(1)+lh3_OrigPos_fig20b(3)+lh2_OrigPos_fig20b(3)+lh1_OrigPos_fig20b(3), lh4_addLeftSpace_fig20b+lh4_OrigPos_fig20b(1)+lh3_OrigPos_fig20b(3)+lh2_OrigPos_fig20b(3)+lh1_OrigPos_fig20b(3) + lh4_OrigPos_fig20b(3)], [lh4_OrigPos_fig20b(2)-yShiftAnnotation_fig20b + lh4_OrigPos_fig20b(4), lh4_OrigPos_fig20b(2)-yShiftAnnotation_fig20b + lh4_OrigPos_fig20b(4)], 'Color', 'black');
lh4titleBottom_fig20b = annotation('line', [lh4_addLeftSpace_fig20b+lh4_OrigPos_fig20b(1)+lh3_OrigPos_fig20b(3)+lh2_OrigPos_fig20b(3)+lh1_OrigPos_fig20b(3), lh4_addLeftSpace_fig20b+lh4_OrigPos_fig20b(1)+lh3_OrigPos_fig20b(3)+lh2_OrigPos_fig20b(3)+lh1_OrigPos_fig20b(3) + lh4_OrigPos_fig20b(3)], [lh4_OrigPos_fig20b(2)-yShiftAnnotation_fig20b,lh4_OrigPos_fig20b(2)-yShiftAnnotation_fig20b], 'Color', 'black');
lh4titleBottom2_fig20b = annotation('line', [lh4_addLeftSpace_fig20b+lh4_OrigPos_fig20b(1)+lh3_OrigPos_fig20b(3)+lh2_OrigPos_fig20b(3)+lh1_OrigPos_fig20b(3), lh4_addLeftSpace_fig20b+lh4_OrigPos_fig20b(1)+lh3_OrigPos_fig20b(3)+lh2_OrigPos_fig20b(3)+lh1_OrigPos_fig20b(3) + lh4_OrigPos_fig20b(3)], [lh4_OrigPos_fig20b(2)+yShiftAnnotation_fig20b*1.2 + lh4_OrigPos_fig20b(4), lh4_OrigPos_fig20b(2)+yShiftAnnotation_fig20b*1.2 + lh4_OrigPos_fig20b(4)], 'Color', 'black');

%% Fifth legend

switch modelChoice
    case {'modelFinal','modelOne','modelTwo','modelThree','modelFour','modelFive','modelSix','modelSeven','modelEight','modelNine'}
        % Copy the axes and plot the fifth legend
        ax5_fig20 = copyobj(ax4_fig20,gcf); %Copies ah1 object (but duplicates data points)
        delete(get(ax5_fig20,'Children')) %Deletes duplicate data points

        % Remove fourth axis visibility 
        set(ax5_fig20, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

        % Plot "dummy" data for fourth legend; tried using single variable with loop but does not work
        hold on
        lh5dump1_fig20b = plot(NaN, NaN, 'k-', 'LineWidth', 0.5,'Parent',ax5_fig20);
        hold off

        % Legend
        lh5_fig20 = legend(ax5_fig20,lh5dump1_fig20b,'Rabinovich and Kalman (2009a), model','Box','off','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Legend.Interpreter);
        lh5_fig20.Units = 'normalized';
        lh5_fig20.Location = 'south';
    case 'modelTen'
% Copy the axes and plot the fifth legend
        ax5_fig20 = copyobj(ax4_fig20,gcf); %Copies ah1 object (but duplicates data points)
        delete(get(ax5_fig20,'Children')) %Deletes duplicate data points

        % Remove fourth axis visibility 
        set(ax5_fig20, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

        % Plot "dummy" data for fourth legend; tried using single variable with loop but does not work
        hold on
        lh5dump1_fig20b = plot(NaN, NaN, 'r-', 'LineWidth', 0.5,'Parent',ax5_fig20);
        hold off

        % Legend
        lh5_fig20 = legend(ax5_fig20,lh5dump1_fig20b,'Rabinovich and Kalman (2009a), model','Box','off','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Legend.Interpreter);
        lh5_fig20.Units = 'normalized';
        lh5_fig20.Location = 'northeast';
        lh5_fig20.IconColumnWidth = 15; %modify legend line width, default: 30
        lh5_fig20.TextColor = 'r';


        % Copy the axes and plot the sixth legend
        ax6_fig20 = copyobj(ax5_fig20,gcf); %Copies ah1 object (but duplicates data points)
        delete(get(ax6_fig20,'Children')) %Deletes duplicate data points

        % Remove fourth axis visibility 
        set(ax6_fig20, 'Color', 'none', 'XTick', [], 'Box', 'Off', 'Visible', 'off')

        % Plot "dummy" data for fourth legend; tried using single variable with loop but does not work
        hold on
        lh6dump1_fig20b = plot(NaN, NaN, 'k-', 'LineWidth', 0.5,'Parent',ax6_fig20);
        hold off
        lh6_fig20 = legend(ax6_fig20,lh6dump1_fig20b,'B2K, model','Box','off','FontSize',figProps.Legend.Labels.Title.FontSize,'Interpreter',figProps.Legend.Interpreter);
        lh6_fig20.Units = 'normalized';
        lh6_fig20.Location = 'northeast';
        lh6_fig20.IconColumnWidth = 15; %modify legend line width, default: 30

        lh5_OrigPos_fig20 = get(lh5_fig20,'Position');
        lh6_OrigPos_fig20 = get(lh6_fig20,'Position');
        lh6_addYSpace = 0.05;
        set(lh6_fig20,'Position',[lh5_OrigPos_fig20(1), lh5_OrigPos_fig20(2)-lh6_addYSpace, lh6_OrigPos_fig20(3), lh6_OrigPos_fig20(4)]);
end

%% Export all figures
%---[Find all figure handles]---
figHandles = findall(groot, 'Type', 'figure');

% % Loop through each figure handle and export graphics in PDF format
for i = 1:length(figHandles)
    % Construct the full path for each figure
    fullFilePath = fullfile(fullBaseExportPath, [scriptName '_pdf_fig' num2str(figHandles(i).Number) '.pdf']);

    % Display the full path
    disp(['Exporting figure to: ' fullFilePath]);

    % Export the figure
    exportgraphics(figHandles(i), fullFilePath,'ColorSpace','rgb','ContentType', 'vector');
end

% % % Loop through each figure handle and export graphics in PDF format
% for i = 1:length(figHandles)
%     % Construct the full path for each figure
%     fullFilePath = fullfile(fullBaseExportPath, [scriptName '_pdf_fig' num2str(figHandles(i).Number) '.pdf']);
% 
%     % Display the full path
%     disp(['Exporting figure to: ' fullFilePath]);
% 
%     % Export the figure
%     print(figHandles(i), fullFilePath,'-vector','-dpdf');
% end

% % % Loop through each figure handle and export graphics in PDF format
% for i = 1:length(figHandles)
%     % Construct the full path for each figure
%     fullFilePath = fullfile(fullBaseExportPath, [scriptName '_pdf_fig' num2str(figHandles(i).Number) '.pdf']);
% 
%     % Display the full path
%     disp(['Exporting figure to: ' fullFilePath]);
% 
%     % Export the figure
%     saveas(figHandles(i), fullFilePath);
%     export_fig fullFilePath
% end

% Loop through each figure handle and export graphics in JPG format at 600 dpi
% for i = 1:length(figHandles)
%     % Construct the full path for each figure
%     fullFilePath = fullfile(fullBaseExportPath, [scriptName '_pdf_fig' num2str(figHandles(i).Number) '.jpg']);
% 
%     % Display the full path
%     disp(['Exporting figure to: ' fullFilePath]);
% 
%     % Export the figure
%     exportgraphics(figHandles(i), fullFilePath,'ColorSpace','rgb','ContentType', 'image','Resolution',600);
% end

%% Save used functions to log .txt file

currScriptNameWithExt = append(currScriptName,'.m');
[fList,pList] = matlab.codetools.requiredFilesAndProducts(currScriptNameWithExt);

logFileName = append(currScriptNameAddDate,'_log','.txt');
logFileFullPath = fullfile(fullBaseExportPath,logFileName);
writelines([newline,'Function run (function_YYMMDD_instance):'],logFileFullPath);
writelines(currScriptNameAddDate,logFileFullPath,"WriteMode","append");

writelines([newline,'Main function:'],logFileFullPath,"WriteMode","append");
writelines(currScriptNameWithExt,logFileFullPath,"WriteMode","append");

writelines([newline,'Functions used:'],logFileFullPath,"WriteMode","append");
for idxfList = 1:length(fList)
    fListFullPath = fList{idxfList};
    [fListFilePath,fListName,FListExt] = fileparts(fListFullPath);
    writelines(fListName,logFileFullPath,"WriteMode","append");
end

writelines([newline,'Required program files:'],logFileFullPath,"WriteMode","append");

for idxpList = 1:length(pList)
    pListName = pList(idxpList).Name;
    pListVer = pList(idxpList).Version;
    pListNameVer = append(pListName,"_",pListVer);
    writelines(pListNameVer,logFileFullPath,"WriteMode","append")
end

%% ---[LOCAL FUNCTIONS]---

%% Function: evaluateParameters - Evaluate values, uncertainties, labels, and legend strings for any parameters
function [paramStruct,refSize,refSizeEach] = evaluateParameters(flagCombine, Parameters, ParametersCombined, dataCase, paramStruct, paramLists, sortingParamName)

%% 0. Choose the data cell array
if flagCombine == 0
    DataArray = Parameters;
    refSize = length(Parameters);
    refSizeEach = refSize / 3;
else
    DataArray = ParametersCombined;
    refSize = length(ParametersCombined);
    refSizeEach = refSize / 3;
end

%% 1. Update the NewStr field for all parameter fields that have a 1x1 Str
allParamFields = fieldnames(paramStruct);
for i = 1:length(allParamFields)
    fld = allParamFields{i};
    if any(strcmp(fld, {'plotData','zSorted'}))
        continue;
    end
    if isfield(paramStruct.(fld), 'Str')
        if numel(paramStruct.(fld).Str) == 1
            paramStruct.(fld) = updateParameterStr(paramStruct.(fld), paramLists);
            idx = find(paramLists.propertyStr == paramStruct.(fld).Str, 1);
            if ~isempty(idx)
                paramStruct.(fld).LegendStr = paramLists.abbrev(idx);
            end
        end
    end
end

%% 2. Process the designated sorting parameter
if isfield(paramStruct.(sortingParamName), 'Str') && (numel(paramStruct.(sortingParamName).Str) > 1)
    [processedParam, zDefStr] = processZParameter(paramStruct.(sortingParamName), paramLists);
    paramStruct.(sortingParamName) = processedParam;
    [zVal, zUncert] = computeZParameterWithUncertainty(DataArray, zDefStr);
    paramStruct.(sortingParamName).Value       = zVal;
    paramStruct.(sortingParamName).Uncertainty = zUncert;
else
    [val, err] = extractParameterData(DataArray, paramStruct.(sortingParamName));
    paramStruct.(sortingParamName).Value       = val;
    paramStruct.(sortingParamName).Uncertainty = err;
end

%% 3. Extract Value & Uncertainty for all other parameters
for i = 1:numel(allParamFields)
    fld = allParamFields{i};
    if any(strcmp(fld, {'plotData','zSorted'})),     continue;  end
    if ~isfield(paramStruct.(fld), 'Str'),           continue;  end
    if strcmp(fld, sortingParamName),                continue;  end

    if numel(paramStruct.(fld).Str) > 1
        [processedParam, zDefStr] = processZParameter(paramStruct.(fld), paramLists);
        paramStruct.(fld)     = processedParam;
        [val, err]           = computeZParameterWithUncertainty(DataArray, zDefStr);
        paramStruct.(fld).Value       = val;
        paramStruct.(fld).Uncertainty = err;
        continue;
    end

    [val, err] = extractParameterData(DataArray, paramStruct.(fld));
    paramStruct.(fld).Value       = val;
    paramStruct.(fld).Uncertainty = err;
end

%% 4. Compute sorted indices for the designated sorting parameter
sortedIdxCell = computeSegmentSortedIndices(paramStruct.(sortingParamName).Value, 3);
paramStruct.zSorted = sortedIdxCell;

%% 5. For every parameter, compute sorted arrays
for i = 1:length(allParamFields)
    fld = allParamFields{i};
    if any(strcmp(fld, {'plotData','zSorted'}))
        continue;
    end
    if isfield(paramStruct.(fld), 'Value') && isnumeric(paramStruct.(fld).Value)
        paramStruct.(fld).ValueSorted = sortSegments(paramStruct.(fld).Value, sortedIdxCell);
        if isfield(paramStruct.(fld), 'Uncertainty') && ~isempty(paramStruct.(fld).Uncertainty)
            paramStruct.(fld).UncertaintySorted = sortSegments(paramStruct.(fld).Uncertainty, sortedIdxCell);
        end
    end
end

%% 6. Set up plotting markers and colors (and sort colors)
[arrayMarker, arrayColor, arrayColorBase, reshapeSize] = computePlottingMarkers(dataCase, flagCombine);
paramStruct.plotData.arrayMarker = arrayMarker;
paramStruct.plotData.arrayColor  = arrayColor;
paramStruct.plotData.arrayColorBase  = arrayColorBase;
paramStruct.plotData.reshapeSize = reshapeSize;
paramStruct.plotData.arrayColorSorted = sortSegments(arrayColor, sortedIdxCell);

end

%% Local Helper Functions

function param = updateParameterStr(param, paramLists)
idx = find(paramLists.propertyStr == param.Str, 1);
if ~isempty(idx)
    if ~isempty(paramLists.dim(idx)) && ~strcmp(paramLists.dim(idx), "")
        param.NewStr = append(paramLists.str(idx), ", ", paramLists.abbrev(idx), " ", paramLists.dim(idx));
    else
        param.NewStr = append(paramLists.str(idx), ", ", paramLists.abbrev(idx));
    end
end
end

function [value, uncert] = extractParameterData(DataArray, param)
n = length(DataArray);
value  = NaN(n, 1);
uncert = NaN(n, 1);
for i = 1:n
    if strcmp(class(DataArray{i}.(param.Str)), 'DimVar')
        value(i)  = DataArray{i}.(param.Str).Value.Value;
        uncert(i) = DataArray{i}.(param.Str).Value.Err;
    else
        value(i)  = DataArray{i}.(param.Str).Value;
        uncert(i) = DataArray{i}.(param.Str).Err;
    end
end
end

function [zParam, zDefStr] = processZParameter(zParam, paramLists)
if numel(zParam.Str) == 1
    zParam.NewStr = "";
    zDefStr = {};
else
    nZ = numel(zParam.Str);
    zIdx = NaN(nZ,1);
    unitDef = cell(nZ,1);
    for zz = 1:nZ
        zIdx(zz) = find(paramLists.propertyStr == zParam.Str(zz), 1);
        unitDef{zz} = paramLists.unitArray.(zParam.Str(zz));
    end
    zUnitMultiplied = cell(nZ,1);
    for zz = 1:nZ
        zUnitMultiplied{zz} = zParam.Unit(zz) * unitDef{zz};
    end
    zUnitCat = cat(3, zUnitMultiplied{:});
    zUnitSum = sum(zUnitCat, 3);
    nzIdx = find(zUnitSum);
    
    zDefStr = cell(length(nzIdx), 4);
    for i = 1:length(nzIdx)
        zDefStr{i,1} = paramLists.unitDef(nzIdx(i));
        zDefStr{i,2} = zUnitSum(nzIdx(i));
        zDefStr{i,3} = nzIdx(i);
        zDefStr{i,4} = paramLists.unitPropertyStr{nzIdx(i)};
    end
    
    legendNum = "";
    legendDen = "";
    for i = 1:size(zDefStr,1)
        exponent = zDefStr{i,2};
        unitCore = extractBetween(zDefStr{i,1}, '$', '$');
        if exponent > 0
            if exponent == 1
                legendNum = append(legendNum, unitCore);
            else
                legendNum = append(legendNum, unitCore, "^", num2str(exponent));
            end
        elseif exponent < 0
            if exponent == -1
                legendDen = append(legendDen, unitCore);
            else
                legendDen = append(legendDen, unitCore, "^", num2str(abs(exponent)));
            end
        end
    end
    zParam.LegendNumStr = legendNum;
    zParam.LegendDenStr = legendDen;
    zParam.LegendStr    = append("$", legendNum, "/", legendDen, "$");
end
end

function [zVal, zUncert] = computeZParameterWithUncertainty(DataArray, zDefStr)
n = length(DataArray);
zVal = NaN(n,1);
zUncert = NaN(n,1);
for j = 1:n
    prodVal = 1;
    sumRelSq = 0;
    for k = 1:size(zDefStr,1)
        unitField = zDefStr{k,4};
        exponent  = zDefStr{k,2};
        factorObj = DataArray{j}.(unitField);
        
        if isa(factorObj, 'DimVar')
            numFactor = factorObj.Value.Value;
            numFactorErr = factorObj.Value.Err;
        elseif isa(factorObj, 'UC')
            numFactor = double(factorObj);
            numFactorErr = factorObj.Err;
        else
            numFactor = factorObj;
            numFactorErr = 0;
        end
        
        prodVal = prodVal * (numFactor ^ exponent);
        if numFactor ~= 0
            relUncert = (exponent * (numFactorErr / numFactor))^2;
            sumRelSq = sumRelSq + relUncert;
        end
    end
    zVal(j) = prodVal;
    zUncert(j) = prodVal * sqrt(sumRelSq);
end
end

function sortedIdxCell = computeSegmentSortedIndices(data, nSegments)
n = length(data);
segLength = n / nSegments;
sortedIdxCell = cell(1, nSegments);
for seg = 1:nSegments
    idx = ((seg-1)*segLength + 1):(seg*segLength);
    [~, sortIdx] = sort(data(idx));
    sortedIdxCell{seg} = sortIdx;
end
end

function sortedData = sortSegments(data, sortedIdxCell)
n = length(data) / numel(sortedIdxCell);
sortedData = [];
for i = 1:numel(sortedIdxCell)
    seg = data((i-1)*n + 1 : i*n);
    sortedData = [sortedData; seg(sortedIdxCell{i})];
end
end

function [arrayMarker, arrayColor, arrayColorBase, reshapeSize] = computePlottingMarkers(dataCase, flagCombine, baseColor)
if nargin < 3 || isempty(baseColor)
    if flagCombine == 0
        baseColor = ...
            ["#9e0142"; %Dark red 1
            "#d53e4f"; %Light red 2
            "#f46d43"; %Dark orange 3
            "#fdae61"; %Orange 4
            "#fee08b"; %Light orange 5
            "#e6f598"; %Lime 6
            "#abdda4"; %Green 7
            "#66c2a5"; %Dark Green 8
            "#3288bd"; %Blue 9
            "#f46d43"; %Dark orange 3
            "#5e4fa2"]; %Purple 11
    else
        baseColor = ...
            ["#9e0142"; %Dark red 1
            "#d53e4f"; %Light red 2
            "#f46d43"; %Dark orange 3
            "#fdae61"; %Orange 4
            "#fee08b"; %Light orange 5
            "#e6f598"; %Lime 6
            "#abdda4"; %Green 7
            "#66c2a5"; %Dark Green 8
            "#3288bd"; %Blue 9
            "#5e4fa2"]; %Purple 10
    end
end
if flagCombine == 0
    fullCount = 11;
else
    fullCount = 10;
end
switch dataCase
    case 'All'
        rowsToKeep = 1:fullCount;
    case 'H>=W'
        if flagCombine == 0
            rowsToKeep = [2,3,5,7,8,9,10,11];
        else
            rowsToKeep = [2,3,5,7,8,9,10];
        end
    case 'H>W'
        rowsToKeep = [2,7,8,9];
    case 'H<=W'
        if flagCombine == 0
            rowsToKeep = [1,3,4,5,6,10,11];
        else
            rowsToKeep = [1,3,4,5,6,10];
        end
    case 'H<W'
        rowsToKeep = [1,4,6];
    otherwise  % 'H=W'
        if flagCombine == 0
            rowsToKeep = [3,5,10,11];
        else
            rowsToKeep = [3,5,10];
        end
end
countRows = numel(rowsToKeep);
arrayMarker = [repmat("square", countRows, 1); repmat("o", countRows, 1); repmat("diamond", countRows, 1)];
arrayColorBase = baseColor(rowsToKeep);
arrayColor = repmat(arrayColorBase, 3, 1);
reshapeSize = [countRows, 3];
end

%% Function: getForce00_Data - Extracts all the appropriate data based on selected case for evaluated force plots
%channel(10)
%-surface(6)
%--botClear(5)
%---mesh(11)
%----solvent(3)
%-----totalDrag
%-----pressureDrag
%-----viscousDrag
%-----totalLift
%-----pressureLift
%-----viscousLift

function [concatenated_values_cell_00, selectedxAxisData, xAxisLabel] = getForce00_Data(Force00_CaseType,channel_00,eval_data_00)

switch Force00_CaseType
    %% xDistInlet
    case "channelALL10_surface1_botClear4_mesh6_solvent1_downstreamsplitChannelAR"
        % Define indices for each level of the nested structure
        indices_A_channel_00 = 1:10; % Index for A
        indices_B_surface_00 = 1; % Index for B
        indices_C_botClear_00 = 4; % Index for C
        indices_D_mesh_00 = 6; % Index for D
        indices_E_solvent_00 = 1; % Indices for E

        % xDistInletObj_array = [1E-4 9E-4 0.0024 0.00301 0.023]; %[m]
        xDistInletObj_array = [0.10 0.90 2.40 3.01 23.00]; %[mm]
        selectedxAxisData = xDistInletObj_array;

        % Initialize cell array to store concatenated values for each A
        concatenated_values_cell_00 = cell(length(indices_A_channel_00), 1);

        % Labels
        % xAxisLabel = 'X-Distance between particle and inlet (m)';
        xAxisLabel = 'X-distance between particle and inlet (mm)';
        
    case "channelALL10_surface1_botClear4_mesh6_solvent2_downstreamsplitChannelAR"
        % Define indices for each level of the nested structure
        indices_A_channel_00 = 1:10; % Index for A
        indices_B_surface_00 = 1; % Index for B
        indices_C_botClear_00 = 4; % Index for C
        indices_D_mesh_00 = 6; % Index for D
        indices_E_solvent_00 = 2; % Indices for E

        % xDistInletObj_array = [1E-4 9E-4 0.0024 0.00301 0.023]; %[m]
        xDistInletObj_array = [0.10 0.90 2.40 3.01 23.00]; %[mm]
        selectedxAxisData = xDistInletObj_array;

        % Initialize cell array to store concatenated values for each A
        concatenated_values_cell_00 = cell(length(indices_A_channel_00), 1);

        % Labels
        % xAxisLabel = 'X-Distance between particle and inlet (m)';
        xAxisLabel = 'X-distance between particle and inlet (mm)';

    case "channelALL10_surface1_botClear4_mesh6_solvent3_downstreamsplitChannelAR"
        % Define indices for each level of the nested structure
        indices_A_channel_00 = 1:10; % Index for A
        indices_B_surface_00 = 1; % Index for B
        indices_C_botClear_00 = 4; % Index for C
        indices_D_mesh_00 = 6; % Index for D
        indices_E_solvent_00 = 3; % Indices for E

        % xDistInletObj_array = [1E-4 9E-4 0.0024 0.00301 0.023]; %[m]
        xDistInletObj_array = [0.10 0.90 2.40 3.01 23.00]; %[mm]
        selectedxAxisData = xDistInletObj_array;

        % Initialize cell array to store concatenated values for each A
        concatenated_values_cell_00 = cell(length(indices_A_channel_00), 1);

        % Labels
        % xAxisLabel = 'X-Distance between particle and inlet (m)';
        xAxisLabel = 'X-distance between particle and inlet (mm)';
    otherwise
        error('Invalid case type')
end

%---[Prepare data]---
% Loop over indices and collect values
for i_A_channel_00 = indices_A_channel_00
    % Initialize concatenated array for this A
    concatenated_values_COMSOL_00 = [];
    for i_B_surface_00 = indices_B_surface_00
        for i_C_botClear_00 = indices_C_botClear_00
            for i_D_xDistInletObj_00 = indices_D_mesh_00
                for i_E_solvent_00 = indices_E_solvent_00
                    concatenated_values_COMSOL_00 = [concatenated_values_COMSOL_00, channel_00(i_A_channel_00).surface_01_11(i_B_surface_00).botClear_01_11(i_C_botClear_00).xDistInletObj_01_11(i_D_xDistInletObj_00).solvent_01_11(i_E_solvent_00).(eval_data_00)];
                end
            end
        end
    end
    % Store concatenated values for this A
    concatenated_values_cell_00{i_A_channel_00} = concatenated_values_COMSOL_00;
end

end

%% Function: getForce01_11Data - Extracts all the appropriate data based on selected case for evaluated force plots
%channel(10)
%-surface(6)
%--botClear(3)
%---xDistInletObj(5)
%----solvent(3)
%-----totalDrag
%-----pressureDrag
%-----viscousDrag
%-----totalLift
%-----pressureLift
%-----viscousLift

function [concatenated_values_cell_01_11, selectedxAxisData, xAxisLabel] = getForce01_11Data(Force01_11CaseType,channel_01_11,eval_data_01_11)

switch Force01_11CaseType
    %% xDistInlet
    case "channelALL10_surface1_botClear1_xDistInletObjALL5_solvent1_splitChannelAR"
        % Define indices for each level of the nested structure
        indices_A_channel_01_11 = 1:10; % Index for A
        indices_B_surface_01_11 = 1; % Index for B
        indices_C_botClear_01_11 = 2; % Index for C
        indices_D_xDistInletObj_01_11 = 1:5; % Index for D
        indices_E_solvent_01_11 = 1; % Indices for E

        % xDistInletObj_array = [1E-4 9E-4 0.0024 0.00301 0.023]; %[m]
        xDistInletObj_array = [0.10 0.90 2.40 3.01 23.00]; %[mm]
        selectedxAxisData = xDistInletObj_array;

        % Initialize cell array to store concatenated values for each A
        concatenated_values_cell_01_11 = cell(length(indices_A_channel_01_11), 1);

        % Labels
        % xAxisLabel = 'X-Distance between particle and inlet (m)';
        xAxisLabel = 'X-distance between particle and inlet (mm)';

    case "channelALL10_surface1_botClear1_xDistInletObjALL5_solvent2_splitChannelAR"
        % Define indices for each level of the nested structure
        indices_A_channel_01_11 = 1:10; % Index for A
        indices_B_surface_01_11 = 1; % Index for B
        indices_C_botClear_01_11 = 2; % Index for C
        indices_D_xDistInletObj_01_11 = 1:5; % Index for D
        indices_E_solvent_01_11 = 2; % Indices for E

        % xDistInletObj_array = [1E-4 9E-4 0.0024 0.00301 0.023]; %[m]
        xDistInletObj_array = [0.10 0.90 2.40 3.01 23.00]; %[mm]
        selectedxAxisData = xDistInletObj_array;

        % Initialize cell array to store concatenated values for each A
        concatenated_values_cell_01_11 = cell(length(indices_A_channel_01_11), 1);

        % Labels
        % xAxisLabel = 'X-Distance between particle and inlet (m)';
        xAxisLabel = 'X-distance between particle and inlet (mm)';

    case "channelALL10_surface1_botClear1_xDistInletObjALL5_solvent3_splitChannelAR"
        % Define indices for each level of the nested structure
        indices_A_channel_01_11 = 1:10; % Index for A
        indices_B_surface_01_11 = 1; % Index for B
        indices_C_botClear_01_11 = 2; % Index for C
        indices_D_xDistInletObj_01_11 = 1:5; % Index for D
        indices_E_solvent_01_11 = 3; % Indices for E

        % xDistInletObj_array = [1E-4 9E-4 0.0024 0.00301 0.023]; %[m]
        xDistInletObj_array = [0.10 0.90 2.40 3.01 23.00]; %[mm]
        selectedxAxisData = xDistInletObj_array;

        % Initialize cell array to store concatenated values for each A
        concatenated_values_cell_01_11 = cell(length(indices_A_channel_01_11), 1);

        % Labels
        % xAxisLabel = 'X-Distance between particle and inlet (m)';
        xAxisLabel = 'X-distance between particle and inlet (mm)';

    %% botClear
    case "channelALL10_surface1_botClearALL3_xDistInletObj1_solvent1_splitChannelAR"
        % Define indices for each level of the nested structure
        indices_A_channel_01_11 = 1:10; % Index for A
        indices_B_surface_01_11 = 1; % Index for B
        indices_C_botClear_01_11 = 1:3; % Index for C
        indices_D_xDistInletObj_01_11 = 1; % Index for D
        indices_E_solvent_01_11 = 1; % Indices for E

        % botClear_array = [1E-5 5E-5 1E-4]; %[m]
        botClear_array = [0.01    0.05    0.10]; %[mm]
        selectedxAxisData = botClear_array;

        % Initialize cell array to store concatenated values for each A
        concatenated_values_cell_01_11 = cell(length(indices_A_channel_01_11), 1);

        % Labels
        % xAxisLabel = 'Bottom clearance between bottom channel wall and bottom surface of particle (m)';
        xAxisLabel = 'Bottom clearance between bottom channel wall and bottom surface of particle (mm)';

    case "channelALL10_surface1_botClearALL3_xDistInletObj1_solvent2_splitChannelAR"
        % Define indices for each level of the nested structure
        indices_A_channel_01_11 = 1:10; % Index for A
        indices_B_surface_01_11 = 1; % Index for B
        indices_C_botClear_01_11 = 1:3; % Index for C
        indices_D_xDistInletObj_01_11 = 1; % Index for D
        indices_E_solvent_01_11 = 2; % Indices for E

        botClear_array = [0.01    0.05    0.10]; %[mm]
        selectedxAxisData = botClear_array;

        % Initialize cell array to store concatenated values for each A
        concatenated_values_cell_01_11 = cell(length(indices_A_channel_01_11), 1);

        % Labels
        xAxisLabel = 'Clearance distance (mm)';

    case "channelALL10_surface1_botClearALL3_xDistInletObj1_solvent3_splitChannelAR"
        % Define indices for each level of the nested structure
        indices_A_channel_01_11 = 1:10; % Index for A
        indices_B_surface_01_11 = 1; % Index for B
        indices_C_botClear_01_11 = 1:3; % Index for C
        indices_D_xDistInletObj_01_11 = 1; % Index for D
        indices_E_solvent_01_11 = 3; % Indices for E

        botClear_array = [0.01    0.05    0.10]; %[mm]
        selectedxAxisData = botClear_array;

        % Initialize cell array to store concatenated values for each A
        concatenated_values_cell_01_11 = cell(length(indices_A_channel_01_11), 1);

        % Labels
        xAxisLabel = 'Clearance distance (mm)';

    otherwise
        error('Invalid case type')
end

%---[Prepare data]---
% Loop over indices and collect values
for i_A_channel_01_11 = indices_A_channel_01_11
    % Initialize concatenated array for this A
    concatenated_values_COMSOL_01_11 = [];
    for i_B_surface_01_11 = indices_B_surface_01_11
        for i_C_botClear_01_11 = indices_C_botClear_01_11
            for i_D_xDistInletObj_01_11 = indices_D_xDistInletObj_01_11
                for i_E_solvent_01_11 = indices_E_solvent_01_11
                    concatenated_values_COMSOL_01_11 = [concatenated_values_COMSOL_01_11, channel_01_11(i_A_channel_01_11).surface_01_11(i_B_surface_01_11).botClear_01_11(i_C_botClear_01_11).xDistInletObj_01_11(i_D_xDistInletObj_01_11).solvent_01_11(i_E_solvent_01_11).(eval_data_01_11)];
                end
            end
        end
    end
    % Store concatenated values for this A
    concatenated_values_cell_01_11{i_A_channel_01_11} = concatenated_values_COMSOL_01_11;
end

end

%% Function: getFigyLim - Loops through all data to determine the maximum y-value to determine maximum yLim to set a consistent y-axis range
function max_yLim = getFigyLim(casesStrALL_01_11, channel_01_11, eval_data_01_11)

concatenated_values_cell_Count = 0;
for idxCasesStrALL_01_11 = 1:length(casesStrALL_01_11)
    concatenated_values_cell_Count = concatenated_values_cell_Count + 1;
    [concatenated_values_cell_01_11{concatenated_values_cell_Count},~,~] = getForce01_11Data(casesStrALL_01_11(idxCasesStrALL_01_11),channel_01_11,eval_data_01_11);
    
    % Set unified ylim based on maximum value of y among all data
    allDataMat = [];
    for idxconcatenated_values_cell_01_11 = 1:length(concatenated_values_cell_01_11)
        dataMat = cell2mat(concatenated_values_cell_01_11{idxconcatenated_values_cell_01_11});
        allDataMat = [allDataMat dataMat];
    end 
end

max_yLim = max(allDataMat,[],"all");
end

%% Function: genForcePlot - Generate force plot
function genForcePlot(xParameterUsed, channel_01_11, eval_data_01_11, xTileIdxs_Cell_ChannelAR_ALL_01_11, colorOrderRatioHWandDmaxSq, arrayColor, colorOrderRegular, figProps)
    %---[Outer, main tiledlayout - to allow for common axis labels]---
    mainTile = tiledlayout(1,1); % Outer layout
    %---[Inner, nested plot Tiledlayout]---
    nestedPlotTile1 = tiledlayout(mainTile,3,3); % Inner layout
    nestedPlotTile1.TileSpacing = 'tight'; 
    nestedPlotTile1.Padding = 'tight';

    %---[Formulate selected plot data into case string]---
    if strcmp(xParameterUsed,'xDistInletObj') 
        casesStrALL_01_11 = ["channelALL10_surface1_botClear1_xDistInletObjALL5_solvent1_splitChannelAR", ...
                             "channelALL10_surface1_botClear1_xDistInletObjALL5_solvent2_splitChannelAR", ...
                             "channelALL10_surface1_botClear1_xDistInletObjALL5_solvent3_splitChannelAR"];
    else % strcmp(xParameterUsed,'botClear') 
        casesStrALL_01_11 = ["channelALL10_surface1_botClearALL3_xDistInletObj1_solvent1_splitChannelAR", ...
                             "channelALL10_surface1_botClearALL3_xDistInletObj1_solvent2_splitChannelAR", ...
                             "channelALL10_surface1_botClearALL3_xDistInletObj1_solvent3_splitChannelAR"];
    end

    %---[Evaluate and get maximum y-value for all sets of plots for selected case towards setting consistent y-lim range]---
    max_yLim = getFigyLim(casesStrALL_01_11, channel_01_11, eval_data_01_11);

    %---[Populate nested tiles with plots, with xParameter as xDistInletObj]---
    if strcmp(xParameterUsed,'xDistInletObj')
        nestedPlotTile1Count = 0;
        concatenated_values_cell_Count = 0;
        for idxCasesStrALL_01_11 = 1:length(casesStrALL_01_11) % Number of solvent cases
            concatenated_values_cell_Count = concatenated_values_cell_Count + 1;
            [concatenated_values_cell_01_11{concatenated_values_cell_Count}, selectedxAxisData, xAxisLabel] = ...
                getForce01_11Data(casesStrALL_01_11(idxCasesStrALL_01_11), channel_01_11, eval_data_01_11); % Get figure data
            for idxChannelARALL_01_11 = 1:length(xTileIdxs_Cell_ChannelAR_ALL_01_11) % Number of different channel aspect ratio ranges
                nestedPlotTile1Count = nestedPlotTile1Count + 1;
                plotCOMSOL01_11_t1 = gobjects(10,1); % Preallocate handles array
                currplotCOMSOL01_ax{nestedPlotTile1Count} = nexttile(nestedPlotTile1, nestedPlotTile1Count, [1 1]);

                %---[Create background tile and axis]---
                breakTiledLayout = tiledlayout(nestedPlotTile1,1,2,'TileSpacing','compact');
                breakTiledLayout.Layout.Tile = nestedPlotTile1Count;
                AxbreakTiledLayout = axes(breakTiledLayout,'XTick',[],'YTick',[],'Box','off');
                AxbreakTiledLayout.Layout.TileSpan = [1 2];
                AxbreakTiledLayout.Visible = 'off';
                set(currplotCOMSOL01_ax{nestedPlotTile1Count},'XTick',[],'YTick',[],'Box','off');

                %---[Get specified indices and indices to plot]---
                specified_indices = xTileIdxs_Cell_ChannelAR_ALL_01_11{idxChannelARALL_01_11};
                indices_to_plot = intersect(colorOrderRatioHWandDmaxSq, specified_indices, 'stable');

                iColor = 0;
                for i_sortedColor_01_11 = colorOrderRatioHWandDmaxSq
                    iColor = iColor + 1;
                    if ismember(i_sortedColor_01_11, indices_to_plot)
                        %---[First interval plotting]---
                        breakAx1{nestedPlotTile1Count} = axes(breakTiledLayout);
                        breakAx1{nestedPlotTile1Count}.Layout.Tile = 1;
                        %---[Plot concatenated values for this set]---
                        plotCOMSOL01_11_t1(i_sortedColor_01_11) = plot(selectedxAxisData, ...
                            concatenated_values_cell_01_11{concatenated_values_cell_Count}{i_sortedColor_01_11}, ...
                            'DisplayName', ['A(', num2str(i_sortedColor_01_11), ')']);

                        %---[Plot settings]---
                        plotCOMSOL01_11_t1(i_sortedColor_01_11).LineWidth = figProps.Line.Line.LineWidth;
                        plotCOMSOL01_11_t1(i_sortedColor_01_11).Marker = 'o';
                        plotCOMSOL01_11_t1(i_sortedColor_01_11).Color = arrayColor(colorOrderRegular(iColor));
                        plotCOMSOL01_11_t1(i_sortedColor_01_11).MarkerFaceColor = arrayColor(colorOrderRegular(iColor));
                        plotCOMSOL01_11_t1(i_sortedColor_01_11).MarkerEdgeColor = 'k';

                        %---[Remove first interval axis]---
                        breakAx1{nestedPlotTile1Count}.Box = 'off';
                        breakAx1{nestedPlotTile1Count}.Visible = 'off';

                        % xlim(breakAx1{nestedPlotTile1Count},[0 0.004]) %[m]
                        xlim(breakAx1{nestedPlotTile1Count},[0 4]) %[mm]
                        breakAx1{nestedPlotTile1Count}.XAxis.Exponent = 0; % Scientific notation
                        xtickformat('%.0f') % Scientific notation %[mm]
                        ylim(breakAx1{nestedPlotTile1Count},[0 max_yLim])
                        breakAx1{nestedPlotTile1Count}.YAxis.Exponent = -6; % Scientific notation

                        xBreakMarkerWidth = 0.2; %[mm]
                        xBreakMarkerHeight = max_yLim/8; %[mm]

                        line([4 - xBreakMarkerWidth, 4 + xBreakMarkerWidth], [-xBreakMarkerHeight, xBreakMarkerHeight],...
                            'Color','k', ...
                            'LineWidth',1.5, ...
                            'Clipping','off');

                        %---[Second interval plotting]---
                        breakAx2{nestedPlotTile1Count} = axes(breakTiledLayout);
                        breakAx2{nestedPlotTile1Count}.Layout.Tile = 2;
                        %---[Plot concatenated values for this A]---
                        plotCOMSOL01_11_t2(i_sortedColor_01_11) = plot(selectedxAxisData, ...
                            concatenated_values_cell_01_11{concatenated_values_cell_Count}{i_sortedColor_01_11}, ...
                            'DisplayName', ['A(', num2str(i_sortedColor_01_11), ')']);

                        %---[Plot settings]---
                        plotCOMSOL01_11_t2(i_sortedColor_01_11).LineWidth = figProps.Line.Line.LineWidth;
                        plotCOMSOL01_11_t2(i_sortedColor_01_11).Marker = 'o';
                        plotCOMSOL01_11_t2(i_sortedColor_01_11).Color = arrayColor(colorOrderRegular(iColor));
                        plotCOMSOL01_11_t2(i_sortedColor_01_11).MarkerFaceColor = arrayColor(colorOrderRegular(iColor));
                        plotCOMSOL01_11_t2(i_sortedColor_01_11).MarkerEdgeColor = 'k';

                        %---[Remove second interval axis]---
                        breakAx2{nestedPlotTile1Count}.Box = 'off';
                        breakAx2{nestedPlotTile1Count}.Visible = 'off';

                        % xlim(breakAx2{nestedPlotTile1Count},[0.022 0.024]) %[m]
                        xlim(breakAx2{nestedPlotTile1Count},[20 24]) %[mm]
                        breakAx2{nestedPlotTile1Count}.XAxis.Exponent = 0; % Scientific notation
                        xtickformat('%.0f') % Scientific notation %[mm]
                        ylim(breakAx2{nestedPlotTile1Count},[0 max_yLim])
                        breakAx2{nestedPlotTile1Count}.YAxis.Exponent = -6; % Scientific notation

                        breakAx2{nestedPlotTile1Count}.YAxis.Visible = 'off'; % Remove extra Y-axis from the second interval that lies in the middle of the combined plot

                        line([20 - xBreakMarkerWidth, 20 + xBreakMarkerWidth], [-xBreakMarkerHeight, xBreakMarkerHeight],...
                            'Color','k', ...
                            'LineWidth',1.5, ...
                            'Clipping','off');

                        %---[Link the axes]---
                        linkaxes([breakAx1{nestedPlotTile1Count}, breakAx2{nestedPlotTile1Count}], 'y')        
                        hold on;

                        %---[Make a single (first) axis visible; all others were previously hidden]---
                        if i_sortedColor_01_11 == indices_to_plot(1)
                            breakAx1{nestedPlotTile1Count}.Visible = 'on'; % Adds default Ax1 ticks and values
                            breakAx2{nestedPlotTile1Count}.Visible = 'on';
                            breakAx1{nestedPlotTile1Count}.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;
                            breakAx2{nestedPlotTile1Count}.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;
                        
                            set(breakAx1{nestedPlotTile1Count},'FontSize',14.6*(8/7.6171));
                            set(breakAx2{nestedPlotTile1Count},'FontSize',14.6*(8/7.6171));

                            set(breakAx1{nestedPlotTile1Count},'TickDir','out');
                            set(breakAx2{nestedPlotTile1Count},'TickDir','out');

                            set(breakAx1{nestedPlotTile1Count},'YTick',0E-05:0.5E-05:2E-05);
                            set(breakAx2{nestedPlotTile1Count},'YTick',0E-05:0.5E-05:2E-05);

                            set(breakAx1{nestedPlotTile1Count},'YMinorTick','on');
                            yAx1 = get(breakAx1{nestedPlotTile1Count},'YAxis');
                            set(yAx1,'MinorTickValues',0.25E-5:0.50E-5:1.75E-05);
                            set(yAx1,'TickLength',[0.1800 0.1950]);

                            xAx1 = get(breakAx1{nestedPlotTile1Count},'XAxis');
                            set(breakAx1{nestedPlotTile1Count},'XMinorTick','on');
                            set(xAx1,'MinorTickValues',[1 3]);
                            set(xAx1,'TickLength',[0.1800 0.1950]);
                            xAx2 = get(breakAx2{nestedPlotTile1Count},'XAxis');
                            set(breakAx2{nestedPlotTile1Count},'XMinorTick','on');
                            set(xAx2,'MinorTickValues',[21 23]);
                            set(xAx2,'TickLength',[0.1800 0.1950]);
                        end
                    end % End of if ismember
                end % End of for i_sortedColor_01_11
            end % End of for idxChannelARALL_01_11
        end % End of for idxCasesStrALL_01_11

    %---[Populate nested tiles with plots, with xParameter as botClear]---
    else % strcmp(xParameterUsed,'botClear') 
        nestedPlotTile1Count = 0;
        concatenated_values_cell_Count = 0;
        for idxCasesStrALL_01_11 = 1:length(casesStrALL_01_11) % Number of solvent cases
            concatenated_values_cell_Count = concatenated_values_cell_Count + 1;
            [concatenated_values_cell_01_11{concatenated_values_cell_Count}, selectedxAxisData, xAxisLabel] = ...
                getForce01_11Data(casesStrALL_01_11(idxCasesStrALL_01_11), channel_01_11, eval_data_01_11);
            for idxChannelARALL_01_11 = 1:length(xTileIdxs_Cell_ChannelAR_ALL_01_11) % Number of different channel aspect ratio ranges
                nestedPlotTile1Count = nestedPlotTile1Count + 1;
                plotCOMSOL01_11_t1 = gobjects(10,1); % Preallocate handles array
                currplotCOMSOL01_ax{nestedPlotTile1Count} = nexttile(nestedPlotTile1, nestedPlotTile1Count, [1 1]);

                %---[Get specified indices and indices to plot]---
                specified_indices = xTileIdxs_Cell_ChannelAR_ALL_01_11{idxChannelARALL_01_11};
                indices_to_plot = intersect(colorOrderRatioHWandDmaxSq, specified_indices, 'stable');

                iColor = 0;
                for i_sortedColor_01_11 = colorOrderRatioHWandDmaxSq
                    iColor = iColor + 1;
                    if ismember(i_sortedColor_01_11, indices_to_plot)
                        % Plot concatenated values for this A
                        plotCOMSOL01_11_t1(i_sortedColor_01_11) = plot(selectedxAxisData, ...
                            concatenated_values_cell_01_11{concatenated_values_cell_Count}{i_sortedColor_01_11}, ...
                            'DisplayName', ['A(', num2str(i_sortedColor_01_11), ')']);
                        plotCOMSOL01_11_t1(i_sortedColor_01_11).LineWidth = figProps.Line.Line.LineWidth;
                        plotCOMSOL01_11_t1(i_sortedColor_01_11).Marker = 'o';
                        plotCOMSOL01_11_t1(i_sortedColor_01_11).Color = arrayColor(colorOrderRegular(iColor));
                        plotCOMSOL01_11_t1(i_sortedColor_01_11).MarkerFaceColor = arrayColor(colorOrderRegular(iColor));
                        plotCOMSOL01_11_t1(i_sortedColor_01_11).MarkerEdgeColor = 'k';
                        hold on;
                    end % End of if ismember
                end % End of for i_sortedColor_01_11
                currplotCOMSOL01_ax{nestedPlotTile1Count}.TickLabelInterpreter = figProps.Axis.Ticks.TickLabelInterpreter;

                set(currplotCOMSOL01_ax{nestedPlotTile1Count},'FontSize',14.6*(8/7.6171));
                
                set(currplotCOMSOL01_ax{nestedPlotTile1Count},'box','off');

                set(currplotCOMSOL01_ax{nestedPlotTile1Count},'TickDir','out');
                
                if strcmp(eval_data_01_11,'totalDrag')
                    set(currplotCOMSOL01_ax{nestedPlotTile1Count},'YTick',0E-05:0.5E-05:2E-05);
                elseif strcmp(eval_data_01_11,'totalLift')
                    set(currplotCOMSOL01_ax{nestedPlotTile1Count},'YTick',0E-05:1.0E-05:2E-05);
                end

                set(currplotCOMSOL01_ax{nestedPlotTile1Count},'YMinorTick','on');
                yAx1 = get(currplotCOMSOL01_ax{nestedPlotTile1Count},'YAxis');
                if strcmp(eval_data_01_11,'totalDrag')
                    set(yAx1,'MinorTickValues',0.25E-5:0.50E-5:1.75E-05);
                elseif strcmp(eval_data_01_11,'totalLift')
                    set(yAx1,'MinorTickValues',0.50E-5:0.50E-5:1.5E-05);
                end
                set(yAx1,'TickLength',[0.06 0.150]);

                set(currplotCOMSOL01_ax{nestedPlotTile1Count},'XMinorTick','on');
                xAx1 = get(currplotCOMSOL01_ax{nestedPlotTile1Count},'XAxis');
                set(xAx1,'MinorTickValues',0.01:0.01:0.09);
                set(xAx1,'TickLength',[0.06 0.150]);
            end % End of for idxChannelARALL_01_11
        end % End of for idxCasesStrALL_01_11

        % Set unified ylim based on maximum value of y among all data
        allDataMat = [];
        for idxconcatenated_values_cell_01_11 = 1:length(concatenated_values_cell_01_11)
            dataMat = cell2mat(concatenated_values_cell_01_11{idxconcatenated_values_cell_01_11});
            allDataMat = [allDataMat, dataMat];
        end

        max_yLim = max(allDataMat,[],"all");

        for idxcurrplotCOMSOL01_ax = 1:length(currplotCOMSOL01_ax)
            ylim(currplotCOMSOL01_ax{idxcurrplotCOMSOL01_ax},[0 max_yLim])
            % ylim(currplotCOMSOL01_ax{idxcurrplotCOMSOL01_ax},[-1 max_yLim])
        end
    end % End of if strcmp(xParameterUsed,'xDistInletObj')

    %---[Add xlabel and ylabel]---
    xlabel(nestedPlotTile1,xAxisLabel,'interpreter',figProps.Axis.Labels.LabelInterpreter,'FontSize',18*(9/9.3908)) %18

    if strcmp(eval_data_01_11,'totalDrag') 
        ylabel(nestedPlotTile1,'Drag force (N)','interpreter',figProps.Axis.Labels.LabelInterpreter,'FontSize',18*(9/9.3908)) %18
    elseif strcmp(eval_data_01_11,'totalLift')
        ylabel(nestedPlotTile1,'Lift force (N)','interpreter',figProps.Axis.Labels.LabelInterpreter,'FontSize',18*(9/9.3908)) %18
    end
end

%% Function: getDefaultFigProperties
function figProps = getDefaultFigProperties(ElsevierFigureColumnSize,AxisFontSizeMultiplier,AxisLabelFontSizeMultiplier,LegendLabelFontSizeMultiplier,LegendFontSizeMultiplier)
    %% [Figure template]
    
    %---[Figure and Axis Handles]---
    
    %---[Prepare data]---
    
    %---[Line properties]---
    %----[Line]----
    %----[Markers]----

    %---[Scatter properties]---
    
    %---[Axis properties]---
    %----[Font]----
    %----[Ticks]----
    %----[Rulers]----
    %----[Grids]----
    %----[Labels]----
    %----[Multiple Plots]----
    %----[Color and Transparency]----
    %----[Box Styling]----
    %----[Position]----
    
    %---[Legend properties]---
    %----[Position and Layout]----
    %----[Labels]----
    %----[Font]----
    %----[Color and Styling]----
    
    %---[Text properties]---
    
    %---[Figure properties]---
    
    %---[Unofficial properties]---

    %% Constant settings (line, axis, legend, figure properties setup)
    %---[Line properties]---
    %----[Line]----
    % def_p_Color = [0 0 0]
    figProps.Line.Line.LineStyle = '-';
    figProps.Line.Line.LineWidth = 1;
    %----[Markers]----
    % def_p_Marker = 'o';
    figProps.Line.Markers.Markersize = 10;
    figProps.Line.Markers.MarkerEdgeColor = 'k';
    % def_p_MarkerFaceColor = 'r';

    %---[Scatter properties]---
    figProps.Scatter.Markers.SizeData = 40; %75, 50, 45, 40, 40, 40
    figProps.Scatter.Markers.MarkerFaceAlpha = 1;
    figProps.Scatter.Markers.MarkerEdgeColor = [0 0 0];

    figProps.Scatter.Markers.SizeData2 = 25; %40, 25, 40, 35, 30, 28, 25
    figProps.Scatter.Markers.LineWidth2 = 2;
    figProps.Scatter.Markers.MarkerEdgeColor2 = [0.6 0.6 0.6];
    
    %---[Axis properties]---
    %----[Font]----
    % all_ax_FontName = 'Arial';
    figProps.Axis.Font.FontWeight = 'bold';
    figProps.Axis.Font.FontSize = 7 * AxisFontSizeMultiplier;
    figProps.Axis.Font.TitleFontWeight = 'bold';
    figProps.Axis.Font.SubtitleFontWeight = 'bold';
    %----[Ticks]----
    % def_ax_XTick = 0:10:100;
    % def_ax_YTick = 0:10:100;
    % def_ax_ZTick = 0:10:100;
    figProps.Axis.Ticks.TickLabelInterpreter = 'latex';
    % def_ax_XMinorTick = 'on';
    % def_ax_YMinorTick = 'on';
    % def_ax_ZMinorTick = 'on';
    % def_ax_TickDir = 'in';
    % def_ax_TickLength = [0.01 0.025];
    %----[Rulers]----
    % def_ax_Xlim = [0 1];
    % def_ax_Ylim = [0 1];
    % def_ax_Zlim = [0 1];
    % def_ax_XScale = 'linear';
    % def_ax_YScale = 'linear';
    % def_ax_ZScale = 'linear';
    %----[Grids]----
    % def_ax_XGrid = 'on';
    % def_ax_YGrid = 'on';
    % def_ax_ZGrid = 'on';
    % def_ax_Layer = 'bottom';
    figProps.Axis.Grids.GridLineStyle = '-';
    figProps.Axis.Grids.GridLineWidth = 1;
    figProps.Axis.Grids.GridColor = [0.6 0.6 0.6]; %[R G B]
    % def_ax_GridAlpha = 0.15;
    % def_ax_XMinorGrid = 'off';
    % def_ax_YMinorGrid = 'off';
    % def_ax_ZMinorGrid = 'off';
    figProps.Axis.Grids.MinorGridLineStyle = ':';
    figProps.Axis.Grids.MinorGridLineWidth = 1;
    figProps.Axis.Grids.MinorGridColor = [0.8 0.8 0.8];
    % def_ax_MinorGridAlpha = 0.25;
    %----[Labels]----
    % def_ax_Title =
    % def_ax_Subtitle =
    % def_ax_TitleHorizontalAlignment = 'center';
    % def_ax_XLabel =
    % def_ax_YLabel =
    % def_ax_ZLabel =
    figProps.Axis.Labels.LabelFontSize = 7 * AxisLabelFontSizeMultiplier;
    figProps.Axis.Labels.LabelFontWeight = 'bold';
    figProps.Axis.Labels.LabelInterpreter = 'latex';
    % def_ax_Legend =
    %----[Multiple Plots]----
    %----[Color and Transparency]----
    %----[Box Styling]----
    % def_ax_Color = [1 1 1];
    % def_ax_LineWidth = 0.5;
    % def_ax_Box = 'off';
    % def_ax_BoxStyle = 'back';
    %----[Position]----
    % def_ax_OuterPosition = [0 0 1 1];
    % def_ax_InnerPosition = [0.1300 0.1100 0.7750 0.8150];
    % def_ax_Position = [0.1300 0.1100 0.7750 0.8150];
    % def_ax_TightInset =
    % def_ax_Units =
    % def_ax_Layout =
    
    %---[Legend properties]---
    %----[Position and Layout]----
    % def_lgd_Location = 'north';
    % def_lgd_Orientation = 'vertical';
    % def_lgd_NumColumns = 1;
    % def_lgd_Position = [0.2 0.6 0.1 0.2];
    % def_lgd_Units = 'normalized';
    % def_lgd_Layout = 'east';
    %----[Labels]----
    % def_lgd_Title.String = 'Example title';
    figProps.Legend.Labels.Title.FontSize = 7 * LegendLabelFontSizeMultiplier;
    figProps.Legend.Interpreter = 'latex';
    % def_lgd_Interpreter = 'latex';
    %----[Font]----
    % def_lgd_FontName = 'Arial';
    figProps.Legend.Fonts.FontSize = 7 * LegendFontSizeMultiplier;
    figProps.Legend.Fonts.FontWeight = 'bold';
    figProps.Legend.Fonts.ItemTokenSize = 10; %Horizontal spacing
    %----[Color and Styling]----
    % main_lgd_TextColor = 'k';
    % main_lgd_Color = [1 1 1];
    figProps.Legend.ColorAndStyling.EdgeColor = [1 1 1];
    % main_lgd_Box = 'on';
    % main_lgd_LineWidth = '1';
    
    %---[Text properties]---
    figProps.Text.Interpreter = 'latex';
    
    %---[Figure properties]---
    % figProps.Figure.Position.Left = 3;
    % figProps.Figure.Position.Bottom = 3;
    % figProps.Figure.Position.Width = 19; % set this parameter and keep it forever %19, 14, 9, 3
    % figProps.Figure.Position.HWRatio = 0.65; % feel free to play with this ratio
    figProps.Figure.Position.Left = 3;
    figProps.Figure.Position.Bottom = 3;
    figProps.Figure.Position.HWRatio = 0.65; % feel free to play with this ratio

    figProps.Line.ErrorLine.LineWidth = 0.5;
    figProps.Line.ErrorLine.CapWidth = 0.125; %0.1, 0.15
    figProps.Line.ErrorLine.Color = [0 0 0];
    figProps.Line.ErrorBar.CapSize = 4;
    
    %---[Unofficial properties]---
    figProps.Unofficial.Axis.XLabelSpacing = 0.0075; %0.4 0.5264; 0.3
    figProps.Unofficial.Axis.YLabelSpacing = 0.0125; %0.4 0.9851; 0.5
    %% Variable settings
    % Define figure size parameters
    switch ElsevierFigureColumnSize
        case '1Column' % Single column (Elsevier 2024 guidelines)
            figProps.Figure.Position.Width = 9; %cm
        case '1.5Column' % 1.5 column (Elsevier 2024 guidelines)
            figProps.Figure.Position.Width = 14; %cm
        case '2Column' % Double column, full width (Elsevier 2024 guidelines)
            figProps.Figure.Position.Width = 19; %cm
        otherwise
            error('RY Error:\nElsevierFigureColumnSize must be 1, 1.5, or 2')
    end
end

%% Function: checkParameters - check if all fields in nParameter structures have assigned values
function missingFields = checkParameters(paramNames)
    missingFields = struct();

    for i = 1:length(paramNames)
        varName = paramNames{i};
        paramStruct = evalin('base', varName);
        fieldsList = fieldnames(paramStruct);
        missing = {};
        
        for j = 1:length(fieldsList)
            currField = fieldsList{j};
            if isempty(paramStruct.(currField))
                missing{end+1} = currField;
            end
        end

        missingFields.(varName) = missing;

        if isempty(missing)
            fprintf('Parameter %s is complete (all fields populated).\n', varName);
        else
            fprintf('Parameter %s has empty fields: %s\n', varName, strjoin(missing, ', '));
        end
    end
end