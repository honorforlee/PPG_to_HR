% Ivan NY HANITRA - Master thesis
%       -- Select events in the frame  --

function [kx_frame,tx_frame,sx_frame, note_x_frame] = frame_select(kx,tx,sx,note_x, frame_init,frame_end)
index_x = find(tx >= frame_init & tx <= frame_end);
kx_frame = kx(index_x);
tx_frame = tx(index_x);
sx_frame = sx(index_x);
note_x_frame = note_x(index_x);



