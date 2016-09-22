function data_list = data_list( folder , filter )
    
    if nargin == 0
        fl = what;
    else
        fl = what(folder);
    end
    if nargin < 2
        ppg_optional = 0; ecg_optional = 1;
    else
        filter = upper(filter);
        ppg_optional = isempty(strfind(upper(filter),'PPG'));
        ecg_optional = isempty(strfind(upper(filter),'ECG'));
    end
    
    data_list = {};
    for f = fl.mat'
        f = f{1}; %#ok<FXSET>
        if exist([f(1:end-3) 'info'],'file')
            d = HBData(f(1:end-4));
            if (d.has_ppg || ppg_optional) && (d.has_ecg || ecg_optional)
                data_list = [data_list d]; %#ok<AGROW>
            end
        end
    end

end
