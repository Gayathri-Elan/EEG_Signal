classdef butterTF
    properties
        n;
        wc;
        wl;
        wh;
    end
    methods
        function obj = butterTF(n,wc,wl,wh)
            obj.n = n;
            obj.wc = wc;
            obj.wh = wh;
            obj.wl = wl;
        end
        function [a,b] = butter_low(obj)
            a = 1;
            if(mod(obj.n,2) == 0)
                b = 1;
            else
                b = [1/obj.wc 1];
            end
            for k = 1:1:obj.n/2
                poly = [1 (2/obj.wc)*sin((pi*(2*k-1))/(2*obj.n)) 1/obj.wc^2];
                b = conv(b, poly);
                a = conv(a,[0 0 1]);
            end
        end

        function [a,b] = butter_high(obj)
            a = 1;
            if(mod(obj.n,2) == 0)
                b = 1;
            else
                b = [obj.wc 1];
                a = conv(a,[1 0]);
            end
            for k = 1:1:obj.n/2
                poly = [1 (2*obj.wc)*sin((pi*(2*k-1))/(2*obj.n)) obj.wc^2];
                b = conv(b, poly);
                a = conv(a,[1 0 0]);
            end
        end

        function [a,b] = butter_band(obj)
            a = 1;
            if(mod(obj.n,2) == 0)
                b = 1;
            else
                b = [1 (obj.wh-obj.wl) obj.wh*obj.wl];
                a = conv(a,[(obj.wh-obj.wl) 0]);
            end
            for k = 1:1:obj.n/2
                poly = [1 (2*(obj.wh-obj.wl))*sin((pi*(2*k-1))/(2*obj.n)) (obj.wh^2+obj.wl^2) 2*obj.wl*obj.wh*(obj.wh-obj.wl)*sin((pi*(2*k-1))/(2*obj.n)) (obj.wl*obj.wh)^2];
                b = conv(b, poly);
                a = conv(a,[(obj.wh-obj.wl)^2 0 0]);
            end
        end
    end
end