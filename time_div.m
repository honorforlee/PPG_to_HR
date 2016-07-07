function [t0_ s0_ t_ s_ ] = time_div(t0,s0,dt0, t,s,dt, frame_length,frame)
range0 = (1 : (frame_length/dt0)) * dt0;
range = (1 : (frame_length/dt)) * dt;

%   - Divide timeline -
for k = 0 : (length(s0) / length(range0)) - 1
    
    t0_div (k+1,:) =  t0(  k*(length(range0)) + 1 : (k+1)*(length(range0)) ) ;
    s0_div(k+1,:)= s0(  k*(length(range0)) + 1 : (k+1)*(length(range0)) ) ;
    
end
t0_ = t0_div(frame,:);    s0_ = s0_div(frame,:);   

for k = 0 : (length(s) / length(range)) - 1
    
    t_div (k+1,:) =  t(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
    s_div(k+1,:)= s(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
    
end
t_ = t_div(frame,:);    s_ = s_div(frame,:);    