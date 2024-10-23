% compareROMvsFOM: MATLAB code for compareROMvsFOM.fig
%      compareROMvsFOM, by itself, creates a new compareROMvsFOM 
%      instance or raises the existing singleton.
%
%      H = compareROMvsFOM returns the handle to a new 
%          compareROMvsFOM or the handle to the existing
%          singleton.

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

function varargout = compareROMvsFOM(varargin)

  % Begin initialization code - DO NOT EDIT
  gui_Singleton = 1;
  gui_State = struct('gui_Name',       mfilename, ...
                     'gui_Singleton',  gui_Singleton, ...
                     'gui_OpeningFcn', @CIV_GUI_OpeningFcn, ...
                     'gui_OutputFcn',  @CIV_GUI_OutputFcn, ...
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

% --- Executes just before compareROMloadedvsFOMloaded is visible.
function CIV_GUI_OpeningFcn(hObject, ~, handles, varargin)
  % This function has no output args, see OutputFcn.
  % hObject    handle to figure
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  % varargin   command line arguments to compareROMloadedvsFOMloaded 

  % Default command line output for compareROMloadedvsFOMloaded
  handles.output = hObject;

  % Update handles structure
  guidata(hObject, handles);

  % RMS Error Text Box Initialization
  global FOMloaded ROMloaded
  FOMloaded = 0; ROMloaded = 0; 
  set(handles.RMS_Error_Text,'String','0.0000 mV');
end

%% --- Outputs from this function are returned to the command line.
function varargout = CIV_GUI_OutputFcn(~, ~, handles) 
  % varargout  cell array for returning output args (see VARARGOUT);
  % hObject    handle to figure
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)

  % Get default command line output from handles structure
  varargout{1} = handles.output;
end


%% Executes on button press: Load FOMloaded results.
function COM_Results_Callback(~, ~, handles)  %#ok<DEFNU>
  % handles    structure with handles and user data (see GUIDATA)
  global Vcell_FOM Vcell_ROM FOMloaded ROMloaded

  [filename, pathname] = uigetfile('*.mat','Pick a MATLAB data file');
  if isequal(filename,0) || isequal(pathname,0)
     disp('User pressed cancel'); return
  end
  disp(['Loading file ', fullfile(pathname, filename)])
  S = load(fullfile(pathname, filename));

  if ~isfield(S,'FOM'),
    disp('File does not contain COMSOL FOM data');
    return % oops -- loaded wrong kind of file?
  end
  
  FOMloaded = 1; Vcell_FOM = S.FOM.Vcell;

  t = 0:length(S.FOM.Vcell)-1;
  axes(handles.axes1); % voltage
  plot(t,S.FOM.Vcell,'r--','linewidth',2); grid on; hold on;
  xlabel('Time (sec)','fontsize',10);
  ylabel('Voltage (DC Volts)','fontsize',10);

  axes(handles.axes2); % cse neg
  plot(t,S.FOM.cse_neg(:,1),'r--', ...
       t,S.FOM.cse_neg(:,end),'r:','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Concentration (mol m^{-3})','fontsize',10);
  title('c_{s,e}','fontsize',12);

  axes(handles.axes3); % cse pos
  plot(t,S.FOM.cse_pos(:,1),'r:', ...
       t,S.FOM.cse_pos(:,end),'r--','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Concentration (mol m^{-3})','fontsize',10);
  title('c_{s,e}','fontsize',12);

  axes(handles.axes4); % ce neg
  [~,indCC] = min(S.FOM.locs.ce_neg_locs);
  [~,indSep] = max(S.FOM.locs.ce_neg_locs);
  plot(t,S.FOM.ce(:,indCC),'r--', ...
       t,S.FOM.ce(:,indSep),'r:','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Concentration (mol m^{-3})','fontsize',10);
  title('c_{e}','fontsize',12);

  axes(handles.axes5); % ce pos
  L = length(S.FOM.locs.ce_pos_locs);
  [~,indCC] = max(S.FOM.locs.ce_pos_locs);
  [~,indSep] = min(S.FOM.locs.ce_pos_locs);
  plot(t,S.FOM.ce(:,end-L+indCC),'r--', ...
       t,S.FOM.ce(:,end-L+indSep),'r:','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Concentration (mol m^{-3})','fontsize',10);
  title('c_{e}','fontsize',12);

  axes(handles.axes6); % phi_s neg
  plot(t,S.FOM.phis_neg(:,end),'r:','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Voltage (V)','fontsize',10);
  title('\phi_{s}','fontsize',12);

  axes(handles.axes7); % phi_s pos
  plot(t,S.FOM.phis_pos(:,1),'r:', ...
       t,S.FOM.phis_pos(:,end),'r--','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Voltage (V)','fontsize',10);
  title('\phi_{s}','fontsize',12);

  axes(handles.axes8); % phi_e neg
  [~,indCC] = min(S.FOM.locs.phie_neg_locs);
  [~,indSep] = max(S.FOM.locs.phie_neg_locs);
  plot(t,S.FOM.phie(:,indCC),'r--',...
       t,S.FOM.phie(:,indSep),'r:','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Voltage (V)','fontsize',10);
  title('\phi_{e}','fontsize',12);

  axes(handles.axes9); % phi_e pos
  L = length(S.FOM.locs.phie_pos_locs);
  [~,indCC] = max(S.FOM.locs.phie_pos_locs);
  [~,indSep] = min(S.FOM.locs.phie_pos_locs);
  plot(t,S.FOM.phie(:,end-L+indCC),'r--', ...
       t,S.FOM.phie(:,end-L+indSep),'r:','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Voltage (V)','fontsize',10);
  title('\phi_{e}','fontsize',12);          

  axes(handles.axes10); % j neg
  plot(t,S.FOM.j_neg(:,1),'r--', ...
       t,S.FOM.j_neg(:,end),'r:','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Flux (mol m^{-2} s^{-1})','fontsize',10);
  title('j','fontsize',12);

  axes(handles.axes11);
  plot(t,S.FOM.j_pos(:,1),'r:', ...
       t,S.FOM.j_pos(:,end),'r--','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Flux (mol m^{-2} s^{-1})','fontsize',10);
  title('j','fontsize',12);

  % RMS Error Display
  if (FOMloaded == 1) && (ROMloaded == 1)
    if length(Vcell_ROM) == length(Vcell_FOM),
      Vrms = sqrt(mean((Vcell_ROM(:)-Vcell_FOM(:)).^2));
      set(handles.RMS_Error_Text,'String',...
          sprintf('%2.4f mV',Vrms*1e3));
    end
  end
end

%% Executes on button press: Load ROMloaded results.
function ROM_Results_Callback(~, ~, handles) %#ok<DEFNU>
  % handles    structure with handles and user data (see GUIDATA)
  global Vcell_FOM Vcell_ROM FOMloaded ROMloaded

  [filename, pathname] = uigetfile('*.mat','Pick a MATLAB data file');
  if isequal(filename,0) || isequal(pathname,0)
    disp('User pressed cancel'); return
  end
  disp(['Loading file ', fullfile(pathname, filename)])

  S = load(fullfile(pathname, filename));
  if ~isfield(S,'ROM'),
    disp('File does not contain simMATLAB data');
    return % oops -- wrong kind of file loaded
  end
  
  ROMloaded = 1; Vcell_ROM = S.ROM.Vcell;
  
  t = 0:length(Vcell_ROM)-1; t=t(:);
  axes(handles.axes1); % cell voltage
  plot(S.ROM.time,S.ROM.Vcell,'b--','linewidth',2); 
  grid on; hold on; xlabel('Time (sec)','fontsize',10);
  ylabel('Voltage (DC Volts)','fontsize',10);

  axes(handles.axes2); % cse in negative electrode
  [minCC,indCC] = min(S.ROM.locs.cse_neg_locs);
  if minCC ~= 0, warning('No cse_neg at current collector'); end
  [maxSep,indSep] = max(S.ROM.locs.cse_neg_locs);
  if maxSep ~= 1, warning('No cse_neg at separator'); end
  plot(t,S.ROM.cse_neg(:,indCC),'b--', ...
       t,S.ROM.cse_neg(:,indSep),'b:','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Concentration (mol m^{-3})','fontsize',10);
  title('c_{s,e}','fontsize',12);

  axes(handles.axes3); % cse in positive electrode
  [minCC,indCC] = min(S.ROM.locs.cse_pos_locs);
  if minCC ~= 0, warning('No cse_pos at current collector'); end
  [maxSep,indSep] = max(S.ROM.locs.cse_pos_locs);
  if maxSep ~= 1, warning('No cse_pos at separator'); end
  plot(t,S.ROM.cse_pos(:,indCC),'b--', ...
       t,S.ROM.cse_pos(:,indSep),'b:','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Concentration (mol m^{-3})','fontsize',10);
  title('c_{s,e}','fontsize',12);

  axes(handles.axes4); % ce in neg
  [~,indCC] = min(S.ROM.locs.ce_neg_locs);
  [~,indSep] = max(S.ROM.locs.ce_neg_locs);
  plot(t,S.ROM.ce(:,indCC),'b--', ...
       t,S.ROM.ce(:,indSep),'b:','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Concentration (mol m^{-3})','fontsize',10);
  title('c_{e}','fontsize',12);

  axes(handles.axes5); % ce in pos
  L = length(S.ROM.locs.ce_pos_locs);
  [~,indCC] = max(S.ROM.locs.ce_pos_locs);
  [~,indSep] = min(S.ROM.locs.ce_pos_locs);
  plot(t,S.ROM.ce(:,end-L+indCC),'b--', ...
       t,S.ROM.ce(:,end-L+indSep),'b:','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Concentration (mol m^{-3})','fontsize',10);
  title('c_{e}','fontsize',12);

  axes(handles.axes6); % phi_s in neg
  [minCC,indCC] = min(S.ROM.locs.phis_neg_locs);
  if minCC == 0, % this is usually missing, but plot if it is there
    plot(t,S.ROM.phi_s_neg(:,indCC),'b--','linewidth',2); hold on
  end
  [maxSep,indSep] = max(S.ROM.locs.phis_neg_locs);
  if maxSep == 1, 
    plot(t,S.ROM.phis_neg(:,indSep),'b:','linewidth',2); hold on
  end
  grid on; xlabel('Time (s)','fontsize',10);
  ylabel('Voltage (V)','fontsize',10);
  title('\phi_{s}','fontsize',12);

  axes(handles.axes7); % phi_s in pos
  [minCC,indCC] = min(S.ROM.locs.phis_pos_locs);
  if minCC == 0, % this is usually missing, but plot if it is there
    plot(t,S.ROM.phis_pos(:,indCC),'b--','linewidth',2); hold on
  end
  [maxSep,indSep] = max(S.ROM.locs.phis_pos_locs);
  if maxSep ~= 1, 
    warning('No phi_s_pos at separator'); 
  else
    plot(t,S.ROM.phis_pos(:,indSep),'b:','linewidth',2); hold on
  end
  grid on;  xlabel('Time (s)','fontsize',10);
  ylabel('Voltage (V)','fontsize',10);
  title('\phi_{s}','fontsize',12);

  axes(handles.axes8); % phi_e neg
  [~,indCC] = min(S.ROM.locs.phie_neg_locs);
  [~,indSep] = max(S.ROM.locs.phie_neg_locs);
  plot(t,S.ROM.phie(:,indCC),'b--',...
       t,S.ROM.phie(:,indSep),'b:','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Voltage (V)','fontsize',10);
  title('\phi_{e}','fontsize',12);

  axes(handles.axes9); % phi_e pos
  L = length(S.ROM.locs.phie_pos_locs);
  [~,indCC] = max(S.ROM.locs.phie_pos_locs);
  [~,indSep] = min(S.ROM.locs.phie_pos_locs);
  plot(t,S.ROM.phie(:,end-L+indCC),'b--', ...
       t,S.ROM.phie(:,end-L+indSep),'b:','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Voltage (V)','fontsize',10);
  title('\phi_{e}','fontsize',12);          

  axes(handles.axes10); % j in neg
  [minCC,indCC] = min(S.ROM.locs.j_neg_locs);
  if minCC ~= 0, warning('No j_neg at current collector'); end
  [maxSep,indSep] = max(S.ROM.locs.j_neg_locs);
  if maxSep ~= 1, warning('No j_neg at separator'); end
  plot(t,S.ROM.j_neg(:,indCC),'b--', ...
       t,S.ROM.j_neg(:,indSep),'b:','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Flux (mol m^{-2} s^{-1})','fontsize',10);
  title('j','fontsize',12);

  axes(handles.axes11); % j in pos
  [minCC,indCC] = min(S.ROM.locs.j_pos_locs);
  if minCC ~= 0, warning('No j_pos at current collector'); end
  [maxSep,indSep] = max(S.ROM.locs.j_pos_locs);
  if maxSep ~= 1, warning('No j_pos at separator'); end
  plot(t,S.ROM.j_pos(:,indCC),'b--', ...
       t,S.ROM.j_pos(:,indSep),'b:','linewidth',2); 
  grid on; hold on; xlabel('Time (s)','fontsize',10);
  ylabel('Flux (mol m^{-2} s^{-1})','fontsize',10);
  title('j','fontsize',12);

  % RMS Error Display
  if (FOMloaded == 1) && (ROMloaded == 1)
    if length(Vcell_ROM) == length(Vcell_FOM),
      Vrms = sqrt(mean((Vcell_ROM(:)-Vcell_FOM(:)).^2));
      set(handles.RMS_Error_Text,'String',...
          sprintf('%2.4f mV',Vrms*1e3));
    end
  end
end

%% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(~, ~, handles) %#ok<DEFNU>
  % handles    structure with handles and user data (see GUIDATA)
  global FOMloaded ROMloaded

  cla(handles.axes1,'reset'); cla(handles.axes2,'reset');
  cla(handles.axes3,'reset'); cla(handles.axes4,'reset');
  cla(handles.axes5,'reset'); cla(handles.axes6,'reset');
  cla(handles.axes7,'reset'); cla(handles.axes8,'reset');
  cla(handles.axes9,'reset'); cla(handles.axes10,'reset');
  cla(handles.axes11,'reset');

  % Clear RMS Error Flags
  FOMloaded = 0; ROMloaded = 0; set(handles.RMS_Error_Text,'String','0.0000 mV');
end
