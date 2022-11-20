classdef Bi_Linear_Transform
    properties
        T;
    end
    methods
        function obj = Bi_Linear_Transform(T)
            obj.T = T;
        end
        function [num, den] = calcBLT(obj,a,b)
            syms s z;
            num1 = poly2sym(a,s);
            den1 = poly2sym(b,s);
            num1 = subs(num1,s,(2/obj.T)*((1-z^-1)/(1+z^-1)));
            den1 = subs(den1,s,(2/obj.T)*((1-z^-1)/(1+z^-1)));
            tf = num1/den1;
            [num,den] = numden(tf);
            num = sym2poly(num);
            den = sym2poly(den);
        end
    end
end

