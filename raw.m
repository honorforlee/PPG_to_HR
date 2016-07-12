


% i = 1;
% for k = i:length(tx_neg)
%     if tx_neg(k) < T - T*0.5
%         
%         if note_x(kx==kx_major(k)) > note_x(kx==kx_major(k+1))      % compare which peak is more relevant
%             kx_major(k+1) = [];
%             tx_major(k+1) = [];
%             sx_major(k+1) = [];
%         else
%             kx_major(k) = [];
%             tx_major(k) = [];
%             sx_major(k) = [];
%         end
%         
%         tx_neg = delta_tx(tx_major);        % recompute tx_neg and T
%         T = mean(delta_tx(tx_major));
%         i=k;                                % start after peak removal
%         break
%     end
%     
% end

kx_major = [1:3:18];
note_x = [1 0.5 1.5 1.6 0.5 2];
tx_major = [0 0.5 1.5 3 3.5 4.5];
tx_neg = delta_tx(tx_major);        % recompute tx_neg and T
T = mean(delta_tx(tx_major));

i = 1;
loop = 0;
loop_ = length(tx_neg);

while loop < loop_ 
for k = i:length(tx_neg)
    if tx_neg(k) < T 
        
        if note_x(k) > note_x(k+1)      % compare which peak is more relevant
            kx_major(k+1) = [];
            tx_major(k+1) = [];
            note_x(k+1) = [];
        else
            kx_major(k) = [];
            tx_major(k) = [];
            note_x(k) = [];
        end
        
        tx_neg = delta_tx(tx_major);        % recompute tx_neg and T
        T = mean(delta_tx(tx_major));
        i=k;                                % start after peak removal
        break
    end
    loop = loop + 1;
end
loop = loop + 1;
end