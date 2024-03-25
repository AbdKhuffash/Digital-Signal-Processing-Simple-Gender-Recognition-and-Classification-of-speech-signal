function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 18-Jan-2024 13:30:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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

end
% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hObject = audiorecorder(44100, 24, 1);
set(handles.rec, 'String' ,'Start speaking for audio');
set(handles.rec, 'foregroundcolor' ,[1 0 0]);
recordblocking(hObject, 2); % record 2 seconds
set(handles.rec, 'String' ,'Audio ended');
set(handles.rec, 'foregroundcolor' ,[0 0 0]);
y = getaudiodata(hObject);
y = y - mean(y);
file_name = sprintf('C:/Users/Abd/Desktop/DSPassingment/GUI/record.wav');
audiowrite(file_name, y, hObject.SampleRate);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.results, 'String' ,'Results will show Here');
training_files_zero_male = dir('C:\Users\Abd\Desktop\DSPassingment\Train\44100M\*.wav');
training_files_zero_female = dir('C:\Users\Abd\Desktop\DSPassingment\Train\44100F\*.wav');
record_file=dir('C:\Users\Abd\Desktop\DSPassingment\GUI\record.wav');
%MALE
data_zero_male_energy = [];
for i = 1:length(training_files_zero_male)
 file_path = strcat(training_files_zero_male(i).folder,'\',training_files_zero_male(i).name);
[y,fs] = audioread(file_path);
energy_zero_male=sum(y .^2);

data_zero_male_energy = [data_zero_male_energy energy_zero_male];
end
energy_zero_male=mean(data_zero_male_energy);

%FEMALE
data_zero_female_energy = [];
for i = 1:length(training_files_zero_female)
 file_path = strcat(training_files_zero_female(i).folder,'\',training_files_zero_female(i).name);
[y,fs] = audioread(file_path);
energy_zero_female=sum(y .^2);

data_zero_female_energy = [data_zero_female_energy energy_zero_female];
end
energy_zero_female=mean(data_zero_female_energy);

%evaluation:
[x,fs] = audioread(file_path);
y_energy  = sum(x.^2);

 %>|
    % Plot in Time Domain
    figure;
    subplot(2,1,1); % Subplot 1 for time domain
    t = linspace(0, length(x)/fs, length(x)); % Time vector
    plot(t, x);
    title('Record File');
    xlabel('Time (seconds)');
    ylabel('Amplitude');

    % Plot in Frequency Domain
    X = fft(x); % Fourier Transform
    P2 = abs(X/length(x));
    P1 = P2(1:floor(length(x)/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(length(x)/2))/length(x); % Frequency vector
    subplot(2,1,2); % Subplot 2 for frequency domain
    plot(f, P1);
    title('Record');
    xlabel('Frequency (Hz)');
    ylabel('|P1(f)|');
    
    ZCR_r1 = mean(abs(diff(sign(y(1:floor(end/3))))))./2;
    ZCR_r2 = mean(abs(diff(sign(y(floor(end/3):floor (end*2/3))))))./2;
    ZCR_r3 = mean(abs(diff(sign(y(floor(end*2/3):end)))))./2;
    
    y_ZCR = [ZCR_r1 ZCR_r2 ZCR_r3];
    
    disp(y_energy);
    disp(y_ZCR);
    

% test if the energy of this file is closer to YES or NO average energies
    if(abs(y_energy-energy_zero_male) < abs(y_energy-energy_zero_female)) 
        %fprintf('Test file [Zero Male] #%d classified as Male saying Zero \n');
        info=['Male:  Energy:' num2str(y_energy)  '                      Zero Crossing 3 parts ' num2str(y_ZCR)]; 
        set(handles.results, 'String' ,info);
    else
        %fprintf('Test file [Zero Male] #%d classified as Female Saying Zero\n');
        info=['Female:  Energy:' num2str(y_energy)  '                      Zero Crossing 3 parts ' num2str(y_ZCR)]; 
        set(handles.results, 'String' ,info);
        
    end

end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
end

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

end
