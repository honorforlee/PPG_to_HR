classdef HBData < handle
    
    properties
        
    end
    properties( SetAccess = private )
        dt
        file
        ppg_k
        ecg_k
        ppg
        ecg
    end
    properties( GetAccess = private )
    end
    properties( Dependent )
        has_ppg
        has_ecg
    end
    
    methods
        
        function obj = HBData( file )
            obj.file = file;
            fid = fopen([file '.info'], 'rt');
            fgetl(fid); fgetl(fid); fgetl(fid);
            obj.dt = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
            obj.dt = obj.dt(2);
            fgetl(fid); s = fgetl(fid);
            obj.ppg_k = 0; obj.ecg_k = 0;
            while s
                [k,signal,~,~,~] = strread(s,'%d%s%f%f%s','delimiter','\t'); %#ok<DSTRRD>
                if strcmp(signal,'II');    obj.ecg_k = k; end
                if strcmp(signal,'PLETH'); obj.ppg_k = k; end
                s = fgetl(fid);
            end
            fclose(fid);
        end
        
        function has_ppg = get.has_ppg( obj )
            has_ppg = ( obj.ppg_k > 0 );
        end
        function has_ecg = get.has_ecg( obj )
            has_ecg = ( obj.ecg_k > 0 );
        end
        
        function load( obj , to_load )
            load([obj.file '.mat']);
            if nargin == 1
                load_ppg = 1; load_ecg = 1;
            else
                load_ppg = ~isempty(strfind(upper(to_load),'PPG'));
                load_ecg = ~isempty(strfind(upper(to_load),'ECG'));
            end
            if ( load_ppg && obj.ppg_k > 0 && isempty(obj.ppg))
                obj.ppg = val(obj.ppg_k, :); %#ok<NODEF>
                obj.ppg = obj.ppg( 1:(find(obj.ppg ~= obj.ppg(end),1,'last')) );
                obj.ppg  = (obj.ppg  - mean(obj.ppg))/sqrt(var(obj.ppg));
            end
            if ( load_ecg && obj.ecg_k > 0  && isempty(obj.ecg))
                obj.ecg = val(obj.ecg_k, :);
                obj.ecg = obj.ecg( 1:(find(obj.ecg ~= obj.ecg(end),1,'last')) );
                obj.ecg  = (obj.ecg  - mean(obj.ecg))/sqrt(var(obj.ecg));
            end
        end
        
        function plot( obj )
            obj.load;
            if obj.has_ppg
                if obj.has_ecg
                    plot( obj.dt*(1:length(obj.ppg)), obj.ppg, 'x-' , obj.dt*(1:length(obj.ecg)), obj.ecg, 'x-' )
                else
                    plot( obj.dt*(1:length(obj.ppg)), obj.ppg, 'x-' )
                end
            else
                if obj.has_ecg
                    plot( obj.dt*(1:length(obj.ecg)), obj.ecg, 'x-' )
                end
            end
        end
        
        function t = ecg_events( obj )
            if ~obj.has_ecg;      return;           end
            if isempty(obj.ecg);  obj.load('ECG');  end
            k = 1 + find( (obj.ecg(2:end-1) > 1) & (obj.ecg(2:end-1) > obj.ecg(1:end-2)) & (obj.ecg(2:end-1) > obj.ecg(3:end)) );
            t = obj.dt*(.5 + k + (obj.ecg(k+1)-obj.ecg(k))./(2*obj.ecg(k)-obj.ecg(k-1)-obj.ecg(k+1)) );
        end
        
    end
    
end
