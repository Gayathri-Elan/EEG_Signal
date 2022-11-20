classdef chebyTF
    methods
        function [a,b] = chebyTF(obj,n, wc, wh, wl, ripple, type)
            e = sqrt(10^(0.1*ripple)-1);
            Tn2 = 1;
            Tn1 = [1 0];
            for k = 1:n-1
                temp = 2*conv([1 0],Tn1);
                size1 = size(temp,2);
                size2 = size(Tn2,2);
                Tn2 = [zeros(1,abs(size1-size2)) Tn2]; %#ok<AGROW>
                T = temp - Tn2;
                Tn2 = Tn1;
                Tn1 = T;
            end
            if(n==1)
                T = Tn1;
            end
            size3 = size(T,2);
            b = e^2*T + [zeros(1,size3-1) 1];
            if(type == "low")
                [a,b] = createLow(obj,b, wc);
            elseif(type == "high")
                [a,b] = createHigh(obj,b, wc);
            else
                [a,b] = createBand(obj,b, wh, wl);
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