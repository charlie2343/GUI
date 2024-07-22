function varargout = GPR_GUI(varargin)
% GPR_GUI MATLAB code for GPR_GUI.fig
%      GPR_GUI, by itself, creates a new GPR_GUI or raises the existing
%      singleton*.
%
%      H = GPR_GUI returns the handle to a new GPR_GUI or the handle to
%      the existing singleton*.
%
%      GPR_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GPR_GUI.M with the given input arguments.
%
%      GPR_GUI('Property','Value',...) creates a new GPR_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GPR_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GPR_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GPR_GUI

% Last Modified by GUIDE v2.5 15-Jul-2020 15:32:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GPR_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GPR_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GPR_GUI is made visible.
function GPR_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GPR_GUI (see VARARGIN)

% Choose default command line output for GPR_GUI
global a;
handles.output = hObject;
% 
% % Update handles structure
guidata(hObject, handles);
a = arduino('com3');



% UIWAIT makes GPR_GUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GPR_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
% PUSHBUTTON 1 - SWEEP
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
init_plots(handles)
set(handles.sweep,'string','IN SWEEP');


obj1 = instrfind('Type', 'visa-gpib', 'RsrcName', 'GPIB0::7::INSTR', 'Tag', '');

% Create the VISA-GPIB object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = visa('AGILENT', 'GPIB1::7::INSTR');
else
    fclose(obj1);
    obj1 = obj1(1);
end
set(obj1,'InputBufferSize',5400);

fopen(obj1);
%set channel and format
fprintf(obj1, 'WAVeform:SOUrce CHANnel1');
fprintf(obj1, 'WAVeform:FORMat WORD');
fprintf(obj1, 'WAVeform:BYTeorder LSBfirst');    %want WORD data

type = query(obj1,'WAVeform:TYPE?');
disp("Acquisition type is " + type);
xinc = str2double(query(obj1, 'WAVeform:XINCrement?')); %get x increment
yinc = str2double(query(obj1, 'WAVeform:YINCrement?'));query(obj1,'*OPC?')%get y increment (unused)
yorg = str2double(query(obj1, 'WAVeform:YORigin?'));   query(obj1,'*OPC?')%get y origin
yref = str2double(query(obj1, 'WAVeform:YREFerence?'));query(obj1,'*OPC?')

fprintf(obj1, 'WAVeform:DATA?');
wfmdata = fread(obj1, 2700, 'int16');query(obj1,'*OPC?')

N_averages = handles.N_averages;
if(~isfloat(handles.N_averages))
    N_averages = str2double(get(handles.N_averages,'string'));
end

N_expected_samples = 50;
t = 0:xinc:(length(wfmdata)-1)*xinc;
global large_data;
large_data = zeros(length(wfmdata),N_expected_samples);
global frametimes;
frametimes = zeros(1,N_expected_samples);
global time;
time = t;

axes(handles.axes2);
ylim([-max(t) 0].*1e9);
axes(handles.axes1);
xlim([0 max(t)].*1e9)

%this for loop takes the average of the for loops matrix and makes a matrix
%of N_samples
k = 1;
t_sample = str2double(get(handles.t_trace,'string'));
image_update_period = 10;
N_WINDOW=4;
window = [zeros(N_WINDOW,1);ones(length(wfmdata)-N_WINDOW,1)];
while (strcmp(get(handles.sweep,'string'), 'IN SWEEP'))
%while(1)
    tic
    tic
    %this for loop creates matrix of a few samples to then be averaged
    dataMat = single(zeros(1,length(wfmdata)));
    for j = 1:N_averages        
        fprintf(obj1, ':DIGitize:CHANNEL1'); 
        fprintf(obj1, ':WAVeform:DATA?');
        dataMat = dataMat + ((fread(obj1, 2700, 'int16')-yref).*yinc.*window + yorg)';  %required for non-ADC units, but
        fprintf(obj1,'*CLS');
    end
    
    meanMatrix = dataMat./N_averages;
    meanMatrix = meanMatrix - mean(meanMatrix);
    large_data(:,k) = meanMatrix;
    %guidata(hObject,handles);
    
    %clear and then plot A-Scan on to first Axis
    if(~mod(k,image_update_period))
        axes(handles.axes1); hold on;
        cla(handles.axes1);  
        line(t.*1e9,meanMatrix)
        ytickformat('%.1f')
        %clear and then plot B-Scan on to second Axis
        %cla(handles.axes2);
        axes(handles.axes2); hold on;  colormap(gray(2^12))
        cla(handles.axes2)
        imagesc(1:size(large_data,2),1e9.*fliplr(t-max(t)),single(large_data));
        drawnow;
    end
    
    
    k=k+1;%loop counter
    
    if( k == N_expected_samples)
       temp = large_data;
       [~,c] = size(large_data);
       large_data = zeros(length(wfmdata),c + 20);
       large_data(:,1:c) = temp;
       
       temp = frametimes;
       frametimes = zeros(1,c+20);
       frametimes(1:c) = temp;       
    end
    
    eT = toc;
    pause(t_sample-eT); %normalize each from to t_sample
    frametimes(k) = toc;
end


% --- Executes on button press in pushbutton2.
%PUSHBUTTON 2 - END SWEEP
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushb    utton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.sweep,'string','NOT IN SWEEP');
global large_data;
global frametimes;
global time;
global scanx;

t=time;
LD = [];

frametimes = frametimes(frametimes~=0);

[r,c] = size(large_data);
for k = 1:c
    if(sum(abs(large_data(:,k))) ~= 0)
       LD(:,k) =  large_data(:,k);
    else
        break;
    end
end
[r,c] = size(LD);
large_data = LD;

% temp = frametimes(1:c);
% frametimes = temp;

speed = str2double(get(handles.speed,'string')).*0.0254;
h = str2double(get(handles.height,'string')).*0.0254;
er = str2double(get(handles.rel_perm,'string'));
t_gnd = 2*h/3e8;
d = zeros(1,length(t));
d(t<t_gnd) = 3e8*(t(t<t_gnd) -t_gnd)./2;
d(t>=t_gnd) = 3e8/sqrt(er)*(t(t>=t_gnd)-t_gnd)./2;
d=d.*100;

% dX = mode(frametimes)*speed;
% x = 0:dX:((c-1)*dX);
x = cumsum(frametimes).*speed;

scanx = x;
LD_bkgr = large_data - mean(large_data,2);
% 
% cla(handles.axes3);
% axes(handles.axes3); hold on;
% imagesc(x,fliplr(d),flipud(LD_bkgr));
% set(gca,'YDir','reverse')
% ylim([min(d) max(d)]);
% xlim([0 max(x)]);
% ylabel('Depth, cm');
% xlabel('Along-track distance, cm');

Npics=5;
pic = round(rand*(Npics-1) + 1);
imgfile=char("upper_left_pictures/pic"+num2str(pic)+".PNG");
axes(handles.axes4);
imgread=imread(imgfile);
imshow(imgread);
drawnow;





function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rel_perm_Callback(hObject, eventdata, handles)
% hObject    handle to rel_perm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rel_perm as text
%        str2double(get(hObject,'String')) returns contents of rel_perm as a double

%acquire from GUI and store handle
handles.permittivity = str2double(get(hObject,'String'));

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function rel_perm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rel_perm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%seting default value
handles.permittivity = 5.7;
guidata(hObject,handles);

% --- Executes on button press in pushbutton3.
%PUSHBUTTON 3- SAVE IMAGE  
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global large_data;
fn = get(handles.edit11,'String');
csvwrite(fn+".csv", large_data);


function pushbutton4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double

%acquire from GUI and store handle
handles.N_traces = str2double(get(hObject,'String'));

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%seting default value
handles.N_traces = 10;
guidata(hObject,handles);

function N_averages_Callback(hObject, eventdata, handles)
% hObject    handle to N_averages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N_averages as text
%        str2double(get(hObject,'String')) returns contents of N_averages as a double

%acquire from GUI and store handle
handles.N_averages = str2double(get(hObject,'String'));

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function N_averages_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_averages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%seting default value
%handles.N_averages = 7;
%guidata(hObject,handles);

%% USER MADE FUNCTIONS


function init_plots(handles)
%this function makes the plots look nice on a button press
set(0,'defaultfigurecolor',[1 1 1 ]);
set(0,'defaultAxesFontSize',16);
set(0,'defaultlinelinewidth',2);
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');

% Npics=5;
% pic = round(rand*(Npics-1) + 1);
% imgfile=char("upper_left_pictures/pic"+num2str(pic)+".PNG");
% axes(handles.axes4);
% imgread=imread(imgfile);
% imshow(imgread);
% drawnow;

XLIM_NS = 25; %scan length
YLIM=15; % magnitude voltage
XLIM_SCAN = 1.5;%meters
YLIM_DEPTH = 50;%cm
%grid on for plot1
axes(handles.axes1); grid on; hold on;
xlabel('Time (ns)');
ylabel('Magnitude (V)');
%xlim([0,XLIM_NS]);
%grid on for plot2
axes(handles.axes2); hold on;
xlabel('Scan Number');
ylabel('Time');
%xlim([0 XLIM_SCAN]);
%ylim([-YLIM_DEPTH 0]);
%caxis([0,50]);


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6. csvwrite
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%N_traces = handles.N_traces;
%N_averages = handles.N_averages;

%SCOPEREAD Code for communicating with an instrument.
%
%   This is the machine generated representation of an instrument control
%   session. The instrument control session comprises all the steps you are
%   likely to take when communicating with your instrument. These steps are:
%   
%       1. Create an instrument object
%       2. Connect to the instrument
%       3. Configure properties
%       4. Write and read data
%       5. Disconnect from the instrument
% 
%   To run the instrument control session, type the name of the file,
%   ScopeRead, at the MATLAB command prompt.
% 
%   The file, SCOPEREAD.M must be on your MATLAB PATH. For additional information 
%   on setting your MATLAB PATH, type 'help addpath' at the MATLAB command 
%   prompt.
% 
%   Example:
%       scoperead;
% 
%   See also SERIAL, GPIB, TCPIP, UDP, VISA, BLUETOOTH, I2C, SPI.
% 
 
%   Creation time: 30-Jul-2018 15:47:54

% Find a VISA-GPIB object.
obj1 = instrfind('Type', 'visa-gpib', 'RsrcName', 'GPIB0::7::INSTR', 'Tag', '');

% Create the VISA-GPIB object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = visa('AGILENT', 'GPIB0::7::INSTR');
else
    fclose(obj1);
    obj1 = obj1(1)
end

% Configure instrument object, obj1.
set(obj1, 'ByteOrder', 'bigEndian');

set(obj1, 'InputBufferSize', 100000);

% Configure instrument object, obj1

set(obj1, 'OutputBufferSize', 512);

% Connect to instrument object, obj1.
fopen(obj1);

% Communicating with instrument object, obj1.
fprintf(obj1, 'WAV:FORMAT WORD');
fprintf(obj1, 'WAV:SOURCE CHAN3');
fprintf(obj1, 'WAV:BYTEORDER MSBF');
dataMat = [];
meanMatrix = [];
large_data = [];
%this for loop takes the average of the for loops matrix and makes a matrix
%of N_samples
for i = 1:N_traces
    %this for loop creates matrix of a few samples to then be averaged
    for j = 1:N_averages
        fprintf(obj1, 'WAV:DATA?');
        data1 = binblockread(obj1, 'int16');
        dataMat(:,j) = data1;
    end
    meanMatrix = (sum(dataMat,2) / j); %average the matrix of data and place in new column of new matrix
    
    
    %ploting
    data2 = query(obj1, 'WAV:YINC?');
    data3 = query(obj1, 'WAV:XINC?')
    data4 = query(obj1, 'WAV:YOR?');
    
    %x = [1:1004]* str2num(data3);
    meanMatrix = meanMatrix * str2num(data2) + str2num(data4);
    csvwrite(get(handles.savefilename,'string')+".csv", meanMatrix); %temp for Sam
    
    large_data = cat(2,meanMatrix, large_data);
    
    handles.large_data = large_data;

    pause(handles.N_delay);

    %clear and then plot A-Scan on to first Axis
    cla(handles.axes1);
    axes(handles.axes1); hold on;
    plot(meanMatrix);

    %clear and then plot B-Scan on to second Axis
    cla(handles.axes2);
    axes(handles.axes2); hold on;

    %[y,x] = size(meanMatrix);
    %fprintf('%d',x);
    imagesc(large_data);
    
    %fprintf('%d', meanMatrix);
    
    %this is to calculate and print percent completed
    completion_percent = strcat(num2str((i/N_traces) * 100), '%');
    set(hObject, 'String', completion_percent);
    
    %save handles
    guidata(hObject,handles);
    
end

set(hObject, 'String', 'Begin Sweep');


% Disconnect all objects.
fclose(obj1);

% Clean up all objects.
delete(obj1);

%save handles
guidata(hObject,handles);



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double




% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', '0');



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double

%acquire from GUI and store handle
handles.N_delay = str2double(get(hObject,'String'));

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%seting default value
handles.N_delay = .1;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

handles.checkbox1 = get(hObject,'Value');
guidata(hObject,handles);


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2

handles.checkbox2 = get(hObject,'Value');
guidata(hObject,handles);



function savefilename_Callback(hObject, eventdata, handles)
% hObject    handle to savefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savefilename as text
%        str2double(get(hObject,'String')) returns contents of savefilename as a double


% --- Executes during object creation, after setting all properties.
function savefilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t_trace_Callback(hObject, eventdata, handles)
% hObject    handle to t_trace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_trace as text
%        str2double(get(hObject,'String')) returns contents of t_trace as a double


% --- Executes during object creation, after setting all properties.
function t_trace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_trace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function speed_Callback(hObject, eventdata, handles)
% hObject    handle to speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of speed as text
%        str2double(get(hObject,'String')) returns contents of speed as a double


% --- Executes during object creation, after setting all properties.
function speed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function height_Callback(hObject, eventdata, handles)
% hObject    handle to height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of height as text
%        str2double(get(hObject,'String')) returns contents of height as a double


% --- Executes during object creation, after setting all properties.
function height_CreateFcn(hObject, eventdata, handles)
% hObject    handle to height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles) %save migrated image
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global image;

fn = get(handles.image_name,'String');
csvwrite(fn+".csv", image);







function image_name_Callback(hObject, eventdata, handles)
% hObject    handle to image_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of image_name as text
%        str2double(get(hObject,'String')) returns contents of image_name as a double


% --- Executes during object creation, after setting all properties.
function image_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11. GENERATE IMAGE
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global large_data;
global time;
global scanx;
global image;

h = str2double(get(handles.height,'string')).*0.0254;
er = str2double(get(handles.rel_perm,'string'));
%I = BScan_kirchoff_migration(S,h,er,t,x,xC,zC, progress)
xC = scanx;
zC = 0.01:5e-3:50e-2;


%perform pre-processing
pp_large_data = preprocess(large_data, time,scanx, hObject, eventdata,handles);
I = BScan_kirchoff_migration2(pp_large_data,h,er,time,scanx,xC,zC,0,1);

image = I;

Ipp = 20.*log10(I); %postprocessed image
DR = str2double(get(handles.drange,'string')); %dynamic range
mI = max(max(Ipp));
cla(handles.axes3);
axes(handles.axes3); hold on;
imagesc(xC,zC,Ipp);
set(gca,'YDir','reverse')
colormap(gca,'default');
caxis([mI - DR, mI])
ylim([0 max(zC)]);
xlim([min(xC) max(xC)])
ylabel('Depth, m');
xlabel('Along-track distance, m');

function data = preprocess(data, time, x, hObject, eventdata,handles)

[n,c] = size(data);
envelope_ = get(handles.envelope,'value');
bkgr_ = get(handles.bkgr,'value');
filter_ = get(handles.filter,'value');
jitter_ = get(handles.jitter,'value');

if(jitter_)
   tlims = [0.5e-9, 2.5e-9];
   data = remove_trigger_jitter(data,time,tlims); 
end

if(bkgr_)
    data = data - mean(data,2); %bkgr subtraction
end

if(envelope_)
   data = abs(hilbert(data));  
end

if(filter_)
   flo =  str2double(get(handles.flo,'string')).*1e6;
   fhi =  str2double(get(handles.fhi,'string')).*1e6;
   
   %scrub the data
   tmp = fhi;
   
   fhi = max([fhi,flo,0]);
   flo = min([flo,tmp]);
   if(flo < 0)
       flo = 0;
   end
   Ts = time(2) - time(1);
   Fs = 1/Ts;
   Fx=x(2)-x(1);
   kx = linspace(-Fx/2,Fx/2,numel(x));
   f = linspace(-Fs/2,Fs/2,n);
   
   Pp_Data = fftshift(fft(data,[],1),1);
   
   Pp_Data = fftshift(fft(Pp_Data,[],2),2);

   Pp_Data((abs(f) < flo) | (abs(f) > fhi),:)=0;    
   KX_MAX=20;
   Pp_Data(:,abs(kx)>KX_MAX)=0;
   data = ifft(ifftshift(Pp_Data,1),[],1);
   data = real(ifft(ifftshift(data,2),[],2));
end



% --- Executes on button press in pushbutton12. fake data
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global large_data;
global frametimes;
global time;
dir = "C:\Users\samwa\OneDrive - UC Davis\2019-04-30_GPR_GUI_v1p1\2019-04-03_DPDImages\";
S = csvread(dir+ "gen5_nodpd.csv");
[n,c] = size(S);
large_data = S;
time = linspace(0,12.5e-9,n);
frametimes = 0.25.*ones(1,c);

axes(handles.axes2); hold on;  colormap(gray(2^12))
cla(handles.axes2)
imagesc(1:size(large_data,2),time,flipud(single(large_data)));
xlim([0, size(large_data,2)]);
ylim([0,max(time)]);


function drange_Callback(hObject, eventdata, handles)
% hObject    handle to drange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of drange as text
%        str2double(get(hObject,'String')) returns contents of drange as a double


% --- Executes during object creation, after setting all properties.
function drange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to drange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in envelope.
function envelope_Callback(hObject, eventdata, handles)
% hObject    handle to envelope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of envelope


% --- Executes on button press in bkgr.
function bkgr_Callback(hObject, eventdata, handles)
% hObject    handle to bkgr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bkgr


% --- Executes on button press in filter.
function filter_Callback(hObject, eventdata, handles)
% hObject    handle to filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filter



function flo_Callback(hObject, eventdata, handles)
% hObject    handle to flo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of flo as text
%        str2double(get(hObject,'String')) returns contents of flo as a double


% --- Executes during object creation, after setting all properties.
function flo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fhi_Callback(hObject, eventdata, handles)
% hObject    handle to fhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fhi as text
%        str2double(get(hObject,'String')) returns contents of fhi as a double


% --- Executes during object creation, after setting all properties.
function fhi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles) %begin multichannel
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global selected_channels;

%survey status of buttons
ch1 = get(handles.ch1,'value');
ch2 = get(handles.ch2,'value');
ch3 = get(handles.ch3,'value');
ch4 = get(handles.ch4,'value');

channels = [ch1 ch2 ch3 ch4];
names = ["CH1", "CH2", "CH3", "CH4"];

sctr=1;
for k = 1:numel(channels)
   if(channels(k))
       selected_channels{sctr} = names{k};
       sctr=sctr+1;
   end
end
N_chan = sctr-1;

set(handles.multisweep,'string','IN SWEEP');


obj1 = instrfind('Type', 'visa-gpib', 'RsrcName', 'GPIB0::7::INSTR', 'Tag', '');
disp('Entering sweep');
% Create the VISA-GPIB object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = visa('AGILENT', 'GPIB0::7::INSTR');
else
    fclose(obj1);
    obj1 = obj1(1);
end
set(obj1,'InputBufferSize',2^20);
fopen(obj1);

% Communicating with instrument object, obj1.
fprintf(obj1, 'DATa:ENCdg ASCII');
fprintf(obj1, 'DATa:WIDth 8');
fprintf(obj1, 'DATa:STARt 1');
fprintf(obj1, 'DATa:STOP 1000000');
wfmdataq = query(obj1, 'CURVe?');
ymult = query(obj1, 'WFMPre:YMUlt?'); 
xinc = query(obj1, 'WFMPre:XINcr?');

wfmdata = str2double(strsplit(wfmdataq,','));
ymult = str2double(ymult);
xinc = str2double(xinc);

%%%%%%%%%%% FOR TEST ONLY 
% dir = "C:\Users\wagner58\Desktop\2019-04-07_ERLVD_Array_X_Y_Sweeps\";
% S = csvread(dir+ "erlvd_array_y_4_plastic.csv");
% wfmdata = S(:,1);
% xinc = 2.44e-11;
%%%%%%%%%%%%%%%%%%%%%%%%


N_averages = handles.N_averages;
if(~isfloat(handles.N_averages))
    N_averages = str2double(get(handles.N_averages,'string'));
end
N_expected_traces=80;
multi_data_c = N_expected_traces;

t = 0:xinc:(length(wfmdata)-1)*xinc;


global multi_data;
global frametimes;
global time;
multi_data = zeros(length(wfmdata),N_expected_traces,N_chan);
frametimes = zeros(1,N_expected_traces);
time = t;

traceCtr=1;
%while(strcmp(get(handles.multisweep,'string'), 'IN SWEEP'))
while(strcmp(get(handles.multisweep,'string'), 'IN SWEEP') && traceCtr <= 82)
    tic;
    
    for k = 1:N_chan
        fprintf(obj1, 'DATa:SOUrce ' + string(selected_channels(k))); %uncomment
        %for proper operation
        %fprintf('DATa:SOUrce ' + string(selected_channels(k))+'\n');
    
        dataMat = zeros(length(wfmdata),1);
    
        for j = 1:N_averages
            dataMat = dataMat + sscanf(query(obj1, 'CURVe?'),"%f,");
            %uncomment
            %dataMat = dataMat + S(:,traceCtr);
        end
        meanMatrix = dataMat./N_averages.*ymult;
        %meanMatrix = dataMat./N_averages;
        meanMatrix = meanMatrix - mean(meanMatrix);
        multi_data(:,traceCtr,k) = meanMatrix; 
    end
    
    pause(1e-2)

    frametimes(traceCtr) = toc;
    traceCtr=traceCtr+1;
    
    if(traceCtr == multi_data_c)
       tmp = multi_data;
       multi_data = zeros(length(wfmdata), multi_data_c+20,N_chan);
       multi_data(:,1:multi_data_c,:)=tmp;
       
       tmpf = frametimes;
       frametimes = zeros(1,multi_data_c+20);
       frametimes(1:multi_data_c) = tmpf;
       
       multi_data_c=multi_data_c+20;
    end
end



% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles) %end multichannel
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.multisweep,'string','NOT IN SWEEP');

global multi_data;
global frametimes;
global time;
global scanx;
[n,c,N_chan] = size(multi_data);

frametimes = frametimes(frametimes~=0);


%uncomment
speed = str2double(get(handles.speed,'string')).*0.0254;
scanx = cumsum(0.25.*ones(size(frametimes))).*speed;


for j = 1:c
      if(sum(abs(multi_data(:,j,N_chan))) ~= 0)
            edited_multi_data(:,j,:) = multi_data(:,j,:);
      end
end



multi_data = edited_multi_data;
[n,c,N_chan] = size(multi_data);

cla(handles.axes2);
axes(handles.axes2); hold on; colormap(gray(2^12))
imagesc(1:c,time,flipud(multi_data(:,:,1)));
xlim([1 c]);
ylim([0 max(time)]);


cla(handles.axes4);
axes(handles.axes4); hold on;
xlabel('Time, s');
histogram(frametimes);




% --- Executes on button press in ch1.
function ch1_Callback(hObject, eventdata, handles)
% hObject    handle to ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ch1


% --- Executes on button press in ch2.
function ch2_Callback(hObject, eventdata, handles)
% hObject    handle to ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ch2


% --- Executes on button press in ch3.
function ch3_Callback(hObject, eventdata, handles)
% hObject    handle to ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ch3


% --- Executes on button press in ch4.
function ch4_Callback(hObject, eventdata, handles)
% hObject    handle to ch4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ch4


% --- Executes on button press in bscan1.
function bscan1_Callback(hObject, eventdata, handles)
% hObject    handle to bscan1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global large_data;
global multi_data;
global selected_channels;


k = min(find(selected_channels=="CH1"));
if(~isempty(k))
    large_data = multi_data(:,:,k);
    update_bscan(hObject, eventdata, handles);
else
    disp('Channel unused.');
end



% --- Executes on button press in bscan2.
function bscan2_Callback(hObject, eventdata, handles)
% hObject    handle to bscan2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global large_data;
global multi_data;
global selected_channels;


k = min(find(selected_channels=="CH2"));
if(~isempty(k))
    large_data = multi_data(:,:,k);
    update_bscan(hObject, eventdata, handles);
else
    disp('Channel unused.');
end

% --- Executes on button press in bscan3.
function bscan3_Callback(hObject, eventdata, handles)
% hObject    handle to bscan3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global large_data;
global multi_data;
global selected_channels;


k = min(find(selected_channels=="CH3"));
if(~isempty(k))
    large_data = multi_data(:,:,k);
    update_bscan(hObject, eventdata, handles);
else
    disp('Channel unused.');
end

% --- Executes on button press in bscan4.
function bscan4_Callback(hObject, eventdata, handles)
% hObject    handle to bscan4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global large_data;
global multi_data;
global selected_channels;


k = min(find(selected_channels=="CH4"));
if(~isempty(k))
    large_data = multi_data(:,:,k);
    update_bscan(hObject, eventdata, handles);
else
    disp('Channel unused.');
end

function update_bscan(hObject, eventdata, handles)
global large_data;
global time;
[n,c] = size(large_data);

cla(handles.axes2);
axes(handles.axes2); hold on; colormap(gray(2^12))
imagesc(1:c,time,flipud(large_data));
xlim([1 c]);
ylim([0 max(time)]);


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles) %save multichannel
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global multi_data;
global selected_channels;

N_chan = numel(unique(selected_channels));
basefn = get(handles.savemulti,'string');

for k = 1:N_chan
   data = multi_data(:,:,k);
   fn = basefn+string(selected_channels{k})+".csv";
   csvwrite(fn,data)
    
end





function savemulti_Callback(hObject, eventdata, handles)
% hObject    handle to savemulti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savemulti as text
%        str2double(get(hObject,'String')) returns contents of savemulti as a double


% --- Executes during object creation, after setting all properties.
function savemulti_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savemulti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in jitter.
function jitter_Callback(hObject, eventdata, handles)
% hObject    handle to jitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of jitter


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles) %pp only
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global large_data;
global time;
global scanx;
pp_large_data = real(preprocess(large_data, time, scanx, hObject, eventdata,handles));
[n,c] = size(pp_large_data);

cla(handles.axes3);
axes(handles.axes3); hold on; colormap(gca,gray(2^12));
imagesc((pp_large_data));
xlim([1 c]);
ylim([1 n]);
mI = min(min(pp_large_data));
mxI = max(max(pp_large_data));
caxis([mI mxI])


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles) %save metadata
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global scanx;
global time;
global frametimes;

fn = get(handles.metafile,'string');

dlmwrite(fn+"-meta_x_t_times.csv", scanx,'delimiter',',');
dlmwrite(fn+"-meta_x_t_times.csv", time,'delimiter',',','-append');
dlmwrite(fn+"-meta_x_t_times.csv", frametimes,'delimiter',',','-append');






function metafile_Callback(hObject, eventdata, handles)
% hObject    handle to metafile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of metafile as text
%        str2double(get(hObject,'String')) returns contents of metafile as a double


% --- Executes during object creation, after setting all properties.
function metafile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to metafile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton22. %GENERATE PULSE (AWG)
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pusghbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global obj2;
obj2 = visa('agilent', 'TCPIP0::localhost::hislip0::INSTR');

fopen(obj2); %open AWG

%import pulse into segment
%first pulse - 200 ps (sigma) differential approx DC-2.5GHz
disp("Reading pulse from " + handles.AWGFilename.String);
try
    fprintf(obj2,"*RST"); %reset settings

    fprintf(obj2,":SOUR:ROSC:RANG RANG1;:SOUR:ROSC:SOUR INT;:OUTP:ROSC:SOUR SCLK1"); %set "Ref Clk Out" switch
    fprintf(obj2,":OUTP:ROSC:SCD 400"); %Set divider by 400 (= 5 MHz)
    fprintf(obj2,":ARM:TRIG:SLOP POS"); %trigger on positive slope
    fprintf(obj2,":ARM:TRIG:LEV 0.25"); %trigger at 0.25V
    fprintf(obj2,":ARM:TRIG:OPER SYNC");%synchronous mode

    fprintf(obj2,":INIT:CONT ON");       %set cont mode
    fprintf(obj2,':TRAC1:DEF 1,1280,0'); %define segment with 1280 samples

    fprintf(obj2,':OUTP1 ON');           %set output 1 on
    fprintf(obj2,':INIT:IMM');           %turn on immediate

    fprintf(obj2,":TRAC1:IMP 1, "+handles.AWGFilename.String+", CSV, IONly, OFF, ALEN");
    fprintf(obj2,':SOUR:FREQ:RAST 64000000000'); %set 64 GS/s sampling rate
    fprintf(obj2,':OUTP1:FILT:FRAT:SCAL 1');   
    fprintf(obj2,':INIT:IMM');           %turn on immediate

    %IF ISSUE HERE, use instrreset (instrfind)
catch
    warning('Error importing waveform');
    fclose(obj2);
end
  







function AWGFilename_Callback(hObject, eventdata, handles)
% hObject    handle to AWGFilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AWGFilename as text
%        str2double(get(hObject,'String')) returns contents of AWGFilename as a double





% --- Executes during object creation, after setting all properties.
function AWGFilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AWGFilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton23. %KILL
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global obj2;
fprintf(obj2,'ABOR');  
fclose(obj2);


% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles) %start
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a;
writeDigitalPin(a,'D26',1);
writeDigitalPin(a,'D28',0);


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles) %stop
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a;
writeDigitalPin(a,'D26',0);
writeDigitalPin(a,'D28',0);



% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles) %reverse
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a;
writeDigitalPin(a,'D26',0);
writeDigitalPin(a,'D28',1);
