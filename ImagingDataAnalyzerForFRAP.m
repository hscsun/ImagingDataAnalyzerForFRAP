function varargout = ImagingDataAnalyzerForFRAP(varargin)
% ImagingDataAnalyzerForFRAP MATLAB code for ImagingDataAnalyzerForFRAP.fig
%      ImagingDataAnalyzerForFRAP, by itself, creates a new ImagingDataAnalyzerForFRAP or raises the existing
%      singleton*.
%
%      H = ImagingDataAnalyzerForFRAP returns the handle to a new ImagingDataAnalyzerForFRAP or the handle to
%      the existing singleton*.
%
%      ImagingDataAnalyzerForFRAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ImagingDataAnalyzerForFRAP.M with the given input arguments.
%
%      ImagingDataAnalyzerForFRAP('Property','Value',...) creates a new ImagingDataAnalyzerForFRAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImagingDataAnalyzerForFRAP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImagingDataAnalyzerForFRAP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ImagingDataAnalyzerForFRAP

% Last Modified by GUIDE v2.5 20-Aug-2015 16:40:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImagingDataAnalyzerForFRAP_OpeningFcn, ...
                   'gui_OutputFcn',  @ImagingDataAnalyzerForFRAP_OutputFcn, ...
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


% --- Executes just before ImagingDataAnalyzerForFRAP is made visible.
function ImagingDataAnalyzerForFRAP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImagingDataAnalyzerForFRAP (see VARARGIN)

% Choose default command line output for ImagingDataAnalyzerForFRAP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ImagingDataAnalyzerForFRAP wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ImagingDataAnalyzerForFRAP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function limitvalue_Callback(hObject, eventdata, handles)
% hObject    handle to limitvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of limitvalue as text
%        str2double(get(hObject,'String')) returns contents of limitvalue as a double


% --- Executes during object creation, after setting all properties.
function limitvalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to limitvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function insteadvalue_Callback(hObject, eventdata, handles)
% hObject    handle to insteadvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of insteadvalue as text
%        str2double(get(hObject,'String')) returns contents of insteadvalue as a double


% --- Executes during object creation, after setting all properties.
function insteadvalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to insteadvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global frame maxv fig m n
[cvname cvpath] = uigetfile({'*.tif';'*.lsm';'*.xlsx';'*.mat';'*.xls'},'File Selector');
 fig=strcat(cvpath,cvname);
 m=str2num(get(handles.channels,'string'));
 n=str2num(get(handles.channel,'string'));

frame=imread(fig,1);
maxv=unique(max(max(frame)));

 axes(handles.axes1);
 imshow(frame,'DisplayRange',[0 maxv]);
 colormap('hot');
 set(handles.result,'string',fig);


% --- Executes on button press in showimage.
function showimage_Callback(hObject, eventdata, handles)
% hObject    handle to showimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 global fig maxv frames imgcor m n
Imageno=str2num(get(handles.imageno,'string'));

frames=imread(fig,((Imageno-1)*m+n));
if isempty(imgcor)
    frames=frames;
else
    frames=frames*imgcor(ceil(Imageno/2));
end

 axes(handles.axes1);
 
 imshow(frames,'DisplayRange',[0 maxv]);
 colormap('hot');
 figure;imshow(frames,'DisplayRange',[0 maxv],'Border','tight');
 colormap('hot');
 set(handles.imageno,'string',Imageno);
 %h=imcontrast(gca);
 


function imageno_Callback(hObject, eventdata, handles)
% hObject    handle to imageno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imageno as text
%        str2double(get(hObject,'String')) returns contents of imageno as a double


% --- Executes during object creation, after setting all properties.
function imageno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function pixelno_Callback(hObject, eventdata, handles)
% hObject    handle to pixelno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixelno as text
%        str2double(get(hObject,'String')) returns contents of pixelno as a double


% --- Executes during object creation, after setting all properties.
function pixelno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixelno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function channels_Callback(hObject, eventdata, handles)
% hObject    handle to channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of channels as text
%        str2double(get(hObject,'String')) returns contents of channels as a double


% --- Executes during object creation, after setting all properties.
function channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function channel_Callback(hObject, eventdata, handles)
% hObject    handle to channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of channel as text
%        str2double(get(hObject,'String')) returns contents of channel as a double


% --- Executes during object creation, after setting all properties.
function channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function regionNo_Callback(hObject, eventdata, handles)
% hObject    handle to regionNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of regionNo as text
%        str2double(get(hObject,'String')) returns contents of regionNo as a double


% --- Executes during object creation, after setting all properties.
function regionNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to regionNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function roi_Callback(hObject, eventdata, handles)
% hObject    handle to roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roi as text
%        str2double(get(hObject,'String')) returns contents of roi as a double


% --- Executes during object creation, after setting all properties.
function roi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function background_Callback(hObject, eventdata, handles)
% hObject    handle to background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of background as text
%        str2double(get(hObject,'String')) returns contents of background as a double


% --- Executes during object creation, after setting all properties.
function background_CreateFcn(hObject, eventdata, handles)
% hObject    handle to background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function flip_Callback(hObject, eventdata, handles)
% hObject    handle to flip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of flip as text
%        str2double(get(hObject,'String')) returns contents of flip as a double


% --- Executes during object creation, after setting all properties.
function flip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function radius_Callback(hObject, eventdata, handles)
% hObject    handle to radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of radius as text
%        str2double(get(hObject,'String')) returns contents of radius as a double


% --- Executes during object creation, after setting all properties.
function radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Loaddata.
function Loaddata_Callback(hObject, eventdata, handles)

%get/read file
global cvpath cvname meanmas nnormaliseddata  measurementNo AftprebleachNo mode imgcor mm 

[cvname cvpath] = uigetfile({'*.txt';'*.xlsx';'*.xls';'*.csv'},'File Selector');
FileName=strcat(cvpath,cvname);
format=strfind(cvname,'txt');
format=size(format);

if format==0
    [data, text, alldata] = xlsread(FileName);
else
    data=textread(FileName,'','whitespace','\t','headerlines',1);
end

set(handles.result,'string',FileName);
imgcor=[]; mm=[];


%get parameters
region= str2num(get(handles.regionNo,'string'));
roi= str2num(get(handles.roi,'string'));
flip= str2num(get(handles.flip,'string'));
bg= str2num(get(handles.background,'string'));
region=region+1;
r= str2num(get(handles.radius,'string'));
AftprebleachNo=r+1;
roi=roi+1;
bg=bg+1;
flip=flip+1;



%load default parameters
dsize=size(data);
dcol=dsize(1,2);
drow=dsize(1,1);
ro1=data(1,:);



%fix the parameters
if bg==roi || bg==flip || bg==1
    bgdata=0;
elseif bg==region
    bg=0;
else
    bg=bg;
end

if flip==roi || bg==flip || flip==1
    flipdata=0;
elseif flip==region
    flip=0;
else
    flip=flip;
end

if roi==region
    roi=0;
else
    roi=roi;
end
ii=0;

%divide data into groups: ROI, reference and background 
for col=1:dcol
    if isnan(ro1(1,col))
        empty=col;
    else
        ii=ii+1;
        iii=ceil(ii/region);
        if rem(ii,region)==1        
            t(:,iii)=data(:,col);    
        elseif rem(ii,region)==roi  
            sig(:,iii)=data(:,col);
        elseif rem(ii,region)==bg  
            bgdata(:,iii)=data(:,col);
        elseif rem(ii,region)==flip  
           flipdata(:,iii)=data(:,col);
        end
    end
end

%pass the prebleach fluorescence intensity into saving function

Imax=mean(sig(1:4,:));
FLIPmax=mean(flipdata(1:4,:));


%prepare the default values for reference and background
if flipdata==0
    flipdata=ones(size(sig))*max(Imax);
end
if bgdata==0
    bgdata=zeros(size(sig));
end

%debackground and normalise the frap data
Imbg=mean(bgdata(1:r,:));
Imflip=mean(flipdata(1:r,:));
dif=Imax-Imbg;
diff=Imflip-Imbg;
ssize=size(dif);
dsize=ssize(1,2);
sigg=sig-bgdata;
siggbg=sigg(AftprebleachNo,:);
ddif=mean(sigg(1:r,:))-siggbg;
totalsig=flipdata-bgdata;

%bring to 0
for jj=1:drow
   sigbg(jj,:)=sigg(jj,:)-siggbg; 
end


%normalise your data
for hh=1:dsize
nnormalisedsig(:,hh)=sigbg(:,hh)/ddif(hh);
To(:,hh)=totalsig(:,hh)/diff(hh);
end
%%normaliseddata=normalisedsig./To;

%correct your data
nnormaliseddata=nnormalisedsig./To;



%calculate the averaged curve
meannorm=mean(nnormaliseddata,2);

%calculate the final recovery level
final=mean(nnormaliseddata((drow-20):drow,:));
meansfinal=mean(final);
hfinal=final/2;
meanfinal=mean(meannorm((drow-20):drow,:));


global  anoname val1 val2 positionsB  FIcell posflips 
tr=size(sig);
hmeanfinal=meanfinal/2;
tsize=size(t);
trow=tsize(1,2);
roend= t(drow,:);
jj=0;


%calculate the half recovery time and Diffusion
for row=1:trow
    if isnan(roend(1,row))
    empty=trow;
    else
        jj=jj+1;
        tt(:,jj)=t(:,row);
    end
    ma(row)=max(find(nnormaliseddata(:,row)<=hfinal(row)));
    minsincludestart=find(nnormaliseddata(:,row)>=hfinal(row));
    minsexcludestart=minsincludestart(minsincludestart>10);    
    mi(row)=min(minsexcludestart);    
end

meanmas=max(find(meannorm<=hmeanfinal));
minincludestart=find(meannorm>=hmeanfinal);
minexcludestart=minincludestart(minincludestart>r);
meanmi=min(minexcludestart);

mm=mean(tt,2);
ma=mm(ma);
mi=mm(mi);
thalf=(ma+mi)/2-mm(r);
thalf=thalf';
meansthalf=mean(thalf);
Ds=0.22*2.5*2.5/meansthalf;

measurementNo=tr(1,1);


%pass data to saving function
if tr(1,2)>1
gStrin = get(handles.controlnameANOVA,'string');
Strins(1:tr(1,2))=cellstr(gStrin);
anoname = [anoname, Strins];
val1 = [val1 thalf];
val2 = [val2 final];
FIcell=[FIcell (Imax+FLIPmax)/2];
else
    set(handles.result,'string',['T-half is ', num2str(thalf), '; Final recovery level is ', num2str(final)]);
end


meanma=mm(meanmas);
meanmi=mm(meanmi);
meanthalf=(meanma+meanmi)/2-mm(r);
D=0.22*2.5*2.5/meanthalf;
drl=D/meanfinal;

%rude fit curve
tcorrect=mm(AftprebleachNo);

yyyy=meanfinal-meanfinal*exp(-0.6931/(meanthalf-tcorrect)*(mm-tcorrect));
yyyy(1:r)=1;


axes(handles.axes3);
plot(mm,nnormaliseddata,'LineWidth',2)
%hold on
%plot(mm,yyyy,'color','k','LineWidth',2)
axis([0,1.02*max(mm),-0.1,1.2])
ylabel('Normalised Fluorescence Intensity (A.U.)','FontSize',10,'FontWeight','normal','Color','k');
xlabel('Time (s)','FontSize',12,'FontWeight','normal','Color','k');
hold off

scrsz = get(groot,'ScreenSize');
rudefitting=figure('Name','FRAP Curve Rude Fitting','NumberTitle','off','Position',[scrsz(3)/10 scrsz(4)/10 scrsz(3)/2 scrsz(4)/2]);
axes('FontSize',12);
plot(mm,nnormaliseddata,'LineWidth',2,'LineStyle','-')
hold on
plot(mm,yyyy,'color','k','LineWidth',2)
axis([0,1.02*max(mm),-0.1,1.2])
ylabel('Normalised Fluorescence Intensity (A.U.)','FontSize',14,'FontWeight','normal','Color','k');
xlabel('Time (s)','FontSize',14,'FontWeight','normal','Color','k');
hold off


%show photobleaching
flip=figure('Name','FLIP','NumberTitle','off','Position',[scrsz(3)/10 scrsz(4)/10 scrsz(3)/2 scrsz(4)/2]);
axes('FontSize',12);
plot(mm,flipdata,'LineWidth',2,'LineStyle','-')
%axis([0,1.02*max(mm),0,1.2])
ylabel('Averaged Fluorescence Intensity (A.U.)','FontSize',14,'FontWeight','normal','Color','k');
xlabel('Time (s)','FontSize',14,'FontWeight','normal','Color','k');
mkdir([cvpath,'anaFRAP/',cvname]);
savenamef=[cvpath, 'anaFRAP\',cvname,'\FLIP'];
print(flip,'-dpng','-r600',savenamef) 

if tr(1,2)>1
%range plot

    mean_nnormaliseddata = mean(nnormaliseddata,2);
    std_nnormaliseddata = std(nnormaliseddata,0,2);
    scrsz = get(groot,'ScreenSize');
    rangeplot=figure('Name','FRAP Curve + std','NumberTitle','off','Position',[scrsz(3)/10 scrsz(4)/10 scrsz(3)/2 scrsz(4)/2]);
    axes('FontSize',14);
    box on
    hold on;
    H1 = plot(mm', mean_nnormaliseddata', 'Color', 'k', 'LineWidth', 2);
    H3 = plot(mm', [mean_nnormaliseddata' - std_nnormaliseddata'; mean_nnormaliseddata' + std_nnormaliseddata'], 'Color', 'm');
    hold on;
    H(2) = shadedErrorBar(mm', nnormaliseddata', {@mean, @(mm) 1*std(mm)  }, '-m', 0);
    legend([H(2).mainLine, H.patch], '\mu', '\sigma', 'Location', 'Northwest','TextColor','g','FontSize',10,'FontWeight','bold');
    ylim([-0.1 1.1])
    ylabel('Normalised Fluorescence (A.U.)','FontSize',14,'FontWeight','normal','Color','k');
    xlabel('Time (s)','FontSize',14,'FontWeight','normal','Color','k');    

    
    savename=[cvpath, 'anaFRAP\',cvname,'\orgRP'];
    print(rangeplot,'-dpng','-r600',savename) 

end



if tr(1,2)==1
    flot=0.1:0.1:1;
    sizeflot=max(size(flot));
    for jjj=1:1:sizeflot
        times(jjj)=(log(1-flot(jjj)))/(-0.6931)*(meanthalf-tcorrect)+tcorrect;
        pos(jjj)=max(find(mm<=times(jjj)));
        posflip(jjj)=flipdata(1,1)/flipdata(pos(jjj),1);
    end
    
    positionsB=[1 AftprebleachNo pos];
    posflips=[1 flipdata(1,1)/flipdata(AftprebleachNo,1) posflip];

    for imgc=1:length(flipdata(:,1))
        imgcor(imgc)=flipdata(1,1)/flipdata(imgc,1);
    end
     
end





function result_Callback(hObject, eventdata, handles)
% hObject    handle to result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of result as text
%        str2double(get(hObject,'String')) returns contents of result as a double


% --- Executes during object creation, after setting all properties.
function result_CreateFcn(hObject, eventdata, handles)
% hObject    handle to result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in profile.
function profile_Callback(hObject, eventdata, handles)
% hObject    handle to profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fig poflips positions  positionsB xxlabe sFD xlabe sree re posflips rer sreer posflipsS AftprebleachNo cvpath cvname
m=str2num(get(handles.channels,'string')); n=str2num(get(handles.channel,'string')); 
cx=str2num(get(handles.cx,'string'));
cy=str2num(get(handles.cy,'string'));
rno=str2num(get(handles.rno,'string'));
psz=str2num(get(handles.pixelsize,'string'));
smr=str2num(get(handles.smoothrange,'string'));
fitting=1;

positionsC=AftprebleachNo:10:positionsB(12);
positionsC=[1,positionsC];

if fitting==0
    [sree,ri,re,sre]=MoRscan(fig,m,n,positions,cx,cy,rno,poflips,smr,fitting);
    
elseif fitting==1
    [sree,ri,re,sre]=MoRscan(fig,m,n,positionsB,cx,cy,rno,posflips,smr,fitting); 
end

    sd=re(:,2);
    se=sree(:,2);
    
    for i=1:12
        rer(:,i)=re(:,i)-sd;
        sreer(:,i)=sree(:,i)-se;
    end

    si=size(re);
    sii=size(ri);
    m1=(-1*rno+1):1:-1;
    m2=1:(rno-1);
    xxlabe=[m1 0 m2]*psz/1000;
    xlabe=(1:si(1,1))*psz/1000;

    axes(handles.axes2);
    plot(xxlabe,sree(:,1:7),'Marker','none','LineWidth',1.2);hold on; plot(xxlabe,sree(:,8:12),'Marker','.','LineWidth',1.2); hold off
    hleg1 = legend('PreBleach','AfterBleach','10% FR','20% FR','30% FR','40% FR','50% FR','60% FR','70% FR','80% FR','90% FR','100% FR');
    axis([1.015*min(xxlabe),1.015*max(xxlabe),-0.05,1.2])
    ylabel('Normalised Fluorescence Intensity (A.U.)','FontSize',12,'FontWeight','normal','Color','k');
    xlabel('Distance to bleaching center (\mum)','FontSize',12,'FontWeight','normal','Color','k');
    set(hleg1,'Location','SouthEast')
    set(hleg1,'Interpreter','none')
    ssi=size(sree);


%saving the radial profiles which I used in the FGF Binding and diffuse
%paper (Raw profile data were used), you can activie this if need.
    out(:,1)=ri(:,1);
    out(:,2)=ri(:,2);
    out(:,3)=ri(:,7);
    out(:,4)=ri(:,12);
    name='profile';
    mkdir([cvpath,'anaFRAP/',cvname]);
    FileName=[cvpath, 'anaFRAP\',cvname,'\rPlot.xlsx'];%[Rpath name Rname];
    xlswrite(FileName,out);




%save half recovery time & final recovery level for many different samples
function saveParameters_Callback(hObject, eventdata, handles)

global anoname val1 val2 FIcell

name='summary';
[Rname Rpath] = uiputfile({'*.csv','*.xlsx'},'save as');
FileName=[Rpath name Rname];
%1 thalf 2 final
%anoname=anoname';
%val1=val1';val2=val2';
xlswrite(FileName,anoname',1,'A1');
xlswrite(FileName,val1',1,'B1');
xlswrite(FileName,val2',1,'C1');
xlswrite(FileName,FIcell',1,'D1');


% --- Executes on button press in zero.
function zero_Callback(hObject, eventdata, handles)
% hObject    handle to zero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global anoname val1 val2 FIcell
anoname=[]; val1=[]; val2=[]; FIcell=[];

function controlnameANOVA_Callback(hObject, eventdata, handles)
% hObject    handle to controlnameANOVA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of controlnameANOVA as text
%        str2double(get(hObject,'String')) returns contents of controlnameANOVA as a double


% --- Executes during object creation, after setting all properties.
function controlnameANOVA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to controlnameANOVA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function cx_Callback(hObject, eventdata, handles)
% hObject    handle to cx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cx as text
%        str2double(get(hObject,'String')) returns contents of cx as a double


% --- Executes during object creation, after setting all properties.
function cx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cy_Callback(hObject, eventdata, handles)
% hObject    handle to cy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cy as text
%        str2double(get(hObject,'String')) returns contents of cy as a double


% --- Executes during object creation, after setting all properties.
function cy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rno_Callback(hObject, eventdata, handles)
% hObject    handle to rno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rno as text
%        str2double(get(hObject,'String')) returns contents of rno as a double


% --- Executes during object creation, after setting all properties.
function rno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pixelsize_Callback(hObject, eventdata, handles)
% hObject    handle to pixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixelsize as text
%        str2double(get(hObject,'String')) returns contents of pixelsize as a double


% --- Executes during object creation, after setting all properties.
function pixelsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% raw profile plot for the bleaching area
function rawplot_Callback(hObject, eventdata, handles)
% hObject    handle to rawplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global xlabe re
%raw half plot
scrsz = get(groot,'ScreenSize');
figure('Name','Raw radial profile for the bleaching area','NumberTitle','off','Position',[scrsz(3)/10 scrsz(4)/10 scrsz(3)/2 scrsz(4)/2]);
axes('FontSize',12);
plot(xlabe,re(:,1:7),'Marker','none','LineWidth',1.2);hold on; plot(xlabe,re(:,8:12),'Marker','.','LineWidth',1.2)
hleg1 = legend('Prebleach','AfterBleach','10% FR','20% FR','30% FR','40% FR','50% FR','60% FR','70% FR','80% FR','90% FR','100% FR');
axis([0,1.015*max(xlabe),0,1.2])
set(hleg1,'Location','SouthEast')
set(hleg1,'Interpreter','none')
ylabel('Normalised Fluorescence Intensity (A.U.)','FontSize',14,'FontWeight','normal','Color','k');
xlabel('Distance to bleaching center (\mum)','FontSize',14,'FontWeight','normal','Color','k');
hold off;



% smoothed profile plot
function smp_Callback(hObject, eventdata, handles)
% hObject    handle to smp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global xxlabe sree
scrsz = get(groot,'ScreenSize');
figure('Name','Smoothed radial profile for the bleaching area','NumberTitle','off','Position',[scrsz(3)/10 scrsz(4)/10 scrsz(3)/2 scrsz(4)/2]);
axes('FontSize',12);
plot(xxlabe,sree(:,1:7),'Marker','none','LineWidth',1.2);hold on; plot(xxlabe,sree(:,8:12),'Marker','.','LineWidth',1.2)
hleg1 = legend('PreBleach','AfterBleach','10% FR','20% FR','30% FR','40% FR','50% FR','60% FR','70% FR','80% FR','90% FR','100% FR');
axis([1.015*min(xxlabe),1.015*max(xxlabe),0,1.2])
set(hleg1,'Location','SouthEast')
set(hleg1,'Interpreter','none')
ylabel('Normalised Fluorescence Intensity (A.U.)','FontSize',14,'FontWeight','normal','Color','k');
xlabel('Distance to bleaching center (\mum)','FontSize',14,'FontWeight','normal','Color','k');
hold off;

function smoothrange_Callback(hObject, eventdata, handles)
% hObject    handle to smoothrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smoothrange as text
%        str2double(get(hObject,'String')) returns contents of smoothrange as a double


% --- Executes during object creation, after setting all properties.
function smoothrange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smoothrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function frange_Callback(hObject, eventdata, handles)
% hObject    handle to frange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frange as text
%        str2double(get(hObject,'String')) returns contents of frange as a double


% --- Executes during object creation, after setting all properties.
function frange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


