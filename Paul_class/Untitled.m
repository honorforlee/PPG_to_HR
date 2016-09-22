%%

addpath('ResultsCheck');
data_folder = '.';

dl = data_list('.','PPG ECG');

%%

disp('------------------------------------------------');
disp('------------------------------------------------');
for filter = {'ppg' 'ecg' 'ecg_ppg' 'ppgecg' 'xxx'}
    filter = filter{1};
    filter = upper(filter);
    ppg_optional = isempty(strfind(filter,'PPG'));
    ecg_optional = isempty(strfind(filter,'ECG'));
    s = [filter '            -> '] ;
    if (~ppg_optional) s = [s 'PPG ']; end
    if (~ecg_optional) s = [s 'ECG ']; end
    disp('------------------------------------------------');
    disp(s);
    for has_ppg = [0 1]
        for has_ecg = [0 1]
            s = ' >  ';
            if ( has_ppg )  s = [s 'PPG ']; end
            if ( has_ecg )  s = [s 'ECG ']; end
            if (has_ppg || ppg_optional) && (has_ecg || ecg_optional)
                disp(s)
            end
        end
    end
end

%%

for filter = {'foo' 'bar'}
    filter
end

%%
if strfind(filter,'ppg') x = 1; end
if strfind(filter,'ppg') x = 1; end
[x y]




%% 
d = dl(1);
d.load;
y = d.ecg;
dt = d.dt;
t = dt*(1:length(y));

k = 1 + find( (y(2:end-1) > 1) & (y(2:end-1) > y(1:end-2)) & (y(2:end-1) > y(3:end)) );
T = dt*(.5 + k + (y(k+1)-y(k))./(2*y(k)-y(k-1)-y(k+1)) );

plot( t,y,'x-',   dt*k , y(k) ,'o' , t(1:end-1)+.5*dt,conv(y,[1;-1],'valid'),'x-',   kron(T,[1 1 1]),repmat([min(y) max(y) nan],size(T)),'-k' )

%%
ttrans = [];
ytrans = [];
dtrans = [];
for i = 1:(length(T)-1)
    if i == 1;  tmin = 0;  else  tmin = (T(i-1)+2*T(i)) / 3;  end
    tmax = (T(i+1)+2*T(i)) / 3;
    k = (t >= tmin) & (t <= tmax); size(y(k))
    ttrans = [ttrans (t(k)+(T(i+1)-T(i))) nan];
    ytrans = [ytrans y(k) nan];
end
plot( t,y,'x-' ...
    , ttrans,ytrans,'-' ...
    , conv(t,     [.5 .5],'valid') , -6 + conv(y,     [1 -1],'valid') , 'x-' ...
    , conv(ttrans,[.5 .5],'valid') , -6 + conv(ytrans,[1 -1],'valid') , '-' );

%%

figure
hold on
for d = dl
    t = d.ecg_events
    plot( conv(t,[.5 .5],'valid') , conv(t,[1 -1],'valid') )
end

%%
d = dl(4)
t = d.ecg_events;
sp(1) = subplot(2,1,1);
d.plot
hold on
plot( conv( d.dt*(1:length(d.ecg)) , [.5 .5] ,'valid') , conv( d.ecg , [1 -1] ,'valid') ,'x-g') 
plot( kron(t,[1 1 1]) , repmat( [ylim nan] , length(t) ) , '-k' )
sp(2) = subplot(2,1,2);
plot( conv(t,[.5 .5],'valid') , conv(t,[1 -1],'valid') )
linkaxes(sp, 'x');

%%
plot(T(2:end)-T(1:end-1),'x-')


%%
d = dl(4);
e = Events(d.ecg);

plot(d.ecg); hold on; e.plot

%%
subplot(2,1,1); histogram(e.top)
subplot(2,1,2); histogram(e.dhi-e.dlo)
%%
d = dl(1); d.load('ecg');
e = Events(d.ecg);
plot(e.top,e.dhi-e.dlo,'.')

%%
for d = dl
    d.load('ecg'); figure
    ax1 = subplot(2,1,1); plot(d.ecg,'x-')
    ax2 = subplot(2,1,2); plot(d.ecg(2:end)-d.ecg(1:end-1),'x-');
    linkaxes([ax1,ax2],'x')
    %e = Events(d.ecg); e_ = Events(-d.ecg); plot(e.top,e.dhi-e.dlo,'.',e_.top,e_.dhi-e_.dlo,'.')
end


