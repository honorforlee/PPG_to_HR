% --- Executes on slider movement. 
function slider_t_int_Callback(hObject, eventdata, handles) 
% hObject handle to slider_t_int (see GCBO) 
% eventdata reserved - to be defined in a future version of MATLAB 
% handles structure with handles and user data (see GUIDATA)

%obtains the slider value from the slider component

sliderValue = get(handles.slider_t_int,'Value');

%puts the slider value into the edit text component 
set(handles.sliderValue_editText,'String', num2str(sliderValue)); 

% Update handles structure 
guidata(hObject, handles);

% Hints: get(hObject,'Value') returns position of slider 
% get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties. 
function slider_t_int_CreateFcn(hObject, eventdata, handles) 
% hObject handle to slider1 (see GCBO) 
% eventdata reserved - to be defined in a future version of MATLAB 
% handles empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background. 
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor')) 
set(hObject,'BackgroundColor',[.9 .9 .9]); 
end

function sliderValue_editText_Callback(hObject, eventdata, handles) 
% hObject handle to sliderValue_editText (see GCBO) 
% eventdata reserved - to be defined in a future version of MATLAB 
% handles structure with handles and user data (see GUIDATA) 
%get the string for the editText component 
sliderValue = get(handles.sliderValue_editText,'String'); 

%convert from string to number if possible, otherwise returns empty 
sliderValue = str2num(sliderValue); 

%if user inputs something is not a number, or if the input is less than 0 
%or greater than 100, then the slider value defaults to 0 
if (isempty(sliderValue) || sliderValue < 920 || sliderValue > 1000) 
set(handles.slider_t_int,'Value',940); 
set(handles.sliderValue_editText,'String','940'); 
else 
set(handles.slider_t_int,'Value',sliderValue); 
end

% Hints: get(hObject,'String') returns contents of sliderValue_editText as text 
% str2double(get(hObject,'String')) returns contents of sliderValue_editText as a double

% --- Executes during object creation, after setting all properties. 
function sliderValue_editText_CreateFcn(hObject, eventdata, handles) 
% hObject handle to sliderValue_editText (see GCBO) 
% eventdata reserved - to be defined in a future version of MATLAB 
% handles empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows. 
% See ISPC and COMPUTER. 
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor')) 
set(hObject,'BackgroundColor','white'); 
end 