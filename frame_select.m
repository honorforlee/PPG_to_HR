% Ivan NY HANITRA - Master thesis
%       -- Select frame of timeline and event  --

function [t0_frame,s0_frame, t_frame,s_frame, kx_frame,tx_frame,sx_frame,sx_N_frame, dhi_frame,dlo_frame,note_x_frame] = frame_select(t0,s0, t,s, kx,tx,sx,sx_N, dhi,dlo,note_x, frame_init,frame_end)
index0 = find(t0 >= frame_init & t0 <= frame_end);
t0_frame = t0(index0);
s0_frame = s0(index0);

index = find(t >= frame_init & t <= frame_end);
t_frame = t(index);
s_frame = s(index);

index_x = find(tx >= frame_init & tx <= frame_end);
kx_frame = kx(index_x);
tx_frame = tx(index_x);
sx_frame = sx(index_x);
sx_N_frame = sx_N(index_x);

dhi_frame = dhi(index_x);
dlo_frame = dlo(index_x);
note_x_frame = note_x(index_x);

