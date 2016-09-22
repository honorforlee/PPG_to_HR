classdef Events < handle
    
    properties( SetAccess = private )
        t
        top
        dhi
        dlo
    end
    
    methods
        
        function obj = Events( s , dt )
            if nargin < 2; dt = 1; end
            d = s(2:end) - s(1:end-1);
            kx = d > 0;
            kx = find(kx(1:end-1) & ~kx(2:end));
            obj.t   = dt * (kx + d(kx)./(d(kx)-d(kx+1)));
            obj.top = s(kx+1);
            obj.dhi = d(kx);
            obj.dlo = d(kx+1);
            for k = 1:length(kx)
                i = kx(k)-1;   while i > 0         && d(i) >= obj.dhi(k); obj.dhi(k) = d(i); i = i-1; end
                i = kx(k)+2;   while i < length(d) && d(i) <= obj.dlo(k); obj.dlo(k) = d(i); i = i+1; end
            end
        end
% 
%         function obj = Events( s , dt )
%             if nargin < 2; dt = 1; end
%             d = s(2:end) - s(1:end-1);
%             obj.top = s([false (d(1:end-1) > 0 & 0 > d(2:end))]);
%             k  = d(1:end-1) >= d(2:end);
%             kt = find( [ k(1)  , ~k(1:end-1) &  k(2:end) , false  ] );
%             kb = find( [ false ,  k(1:end-1) & ~k(2:end) , k(end) ] );
%             k = (d(kt) > 0 & d(kb) < 0);
%             kt = kt(k);
%             kb = kb(k);
%  
%             
% t = 1:80;
% d = sin(.002*pi*t.^2);
% k0 = find( d(1:end-1) > 0 & 0 > d(2:end) );
% k  = d(1:end-1) >= d(2:end);
% kt = find( [ k(1)  , ~k(1:end-1) &  k(2:end) , false  ] );
% kb = find( [ false ,  k(1:end-1) & ~k(2:end) , k(end) ] );
% k = (d(kt) > 0 & d(kb) < 0);
% kt = kt(k);
% kb = kb(k);
% plot([0 80],[0 0],':k',t,d,'x-',t(k0),d(k0+1),'x-',t(kt),d(kt),'o',t(kb),d(kb),'o')
% 
%             
%             
%             
%             
%             kx = d > 0;
%             kx = find(kx(1:end-1) & ~kx(2:end));
%             obj.t   = dt * (kx + d(kx)./(d(kx)-d(kx+1)));
%             obj.top = s(kx+1);
%             obj.dhi = d(kx);
%             obj.dlo = d(kx+1);
%             for k = 1:length(kx)
%                 i = kx(k)-1;   while i > 0         && d(i) >= obj.dhi(k); obj.dhi(k) = d(i); i = i-1; end
%                 i = kx(k)+2;   while i < length(d) && d(i) <= obj.dlo(k); obj.dlo(k) = d(i); i = i+1; end
%             end
%         end

        function plot( obj , col )
            if nargin < 2; col = 'r'; end
            l = nan(1,3*length(obj.dhi)); l(1:3:end) = obj.dlo; l(2:3:end) = obj.dhi;
            if isrow(obj.t); tl = kron(obj.t,[1 1 nan]); else tl = kron(obj.t',[1 1 nan]); end
            plot( obj.t,obj.top,'+',obj.t,obj.dhi,'^',obj.t,obj.dlo,'v' ...
                , tl, l, ':' ...
                ,'MarkerSize',4,'MarkerEdgeColor',col,'Color',col)
        end
        
        
        
    end
end