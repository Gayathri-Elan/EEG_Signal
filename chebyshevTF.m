classdef chebyshevTF
    properties
        n;
        wc;
        wh;
        wl;
        ripple;
        type;
    end
    methods
        function obj = chebyshevTF(n,wc,wh,wl,ripple,type)
            obj.n = n;
            obj.wc = wc;
            obj.wh = wh;
            obj.wl = wl;
            obj.ripple = ripple;
            obj.type = type;
        end
        function [a,b] = cheby_TF(obj)
            e = sqrt(10^(0.1*obj.ripple)-1);
            Tn2 = 1;
            Tn1 = [1 0];
            for k = 1:obj.n-1
                temp = 2*conv([1 0],Tn1);
                size1 = size(temp,2);
                size2 = size(Tn2,2);
                Tn2 = [zeros(1,abs(size1-size2)) Tn2]; %#ok<AGROW>
                T = temp - Tn2;
                Tn2 = Tn1;
                Tn1 = T;
            end
            if(obj.n==1)
                T = Tn1;
            end
            size3 = size(T,2);
            b = e^2*T + [zeros(1,size3-1) 1];
            if(obj.type == "low")
                [a,b] = createLow(obj,b, obj.wc);
            elseif(obj.type == "high")
                [a,b] = createHigh(obj,b, obj.wc);
            else
                [a,b] = createBand(obj,b, obj.wh, obj.wl);
            end
        end

        function [num,den] = createLow(~,b, wc)
            syms w;
            den = poly2sym(b, w);
            den = subs(den, w, w/wc);
            den = sym2poly(den);
            num = 1;
        end

        function [num,den] = createHigh(~,b, wc)
            syms w;
            den = poly2sym(b, w);
            den = subs(den, w, -wc/w);
            [den,num] = numden(den);
            den = sym2poly(den);
            num = sym2poly(num);
        end

        function [num,den] = createBand(~,b, wh, wl)
            syms w;
            den = poly2sym(b, w);
            den = subs(den, w, (w^2-(wl*wh))/(w*(wh-wl)));
            [den,num] = numden(den);
            den = sym2poly(den);
            num = sym2poly(num);
        end
    end
end