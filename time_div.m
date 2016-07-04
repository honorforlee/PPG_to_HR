function [t_ s_] = time_div(t,s,dt, frame_length,frame)
range = (1 : (frame_length/dt)) * dt;

%   - Divide timeline -
for k = 0 : length(s) / length(range) - 1
    
    t_div (k+1,:) =  t(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
    s_div(k+1,:)= s(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
    
end

t_ = t_div(frame,:);    s_ = s_div(frame,:);    