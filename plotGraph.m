classdef plotGraph
    properties
        T;
    end
    methods
        function obj = plotGraph(T)
            obj.T = T;
        end
        function [] = plotTime(obj,X,name)
            figure('Name',name);
            t = calcTime(obj,X);
            plot(t,X);
            xlabel('time');
            ylabel('amplitude');
            hold on
        end
        function [] = plotFreq(obj,X)
            [a,b] = calcZ(obj,X);
            freqz(a,b,'whole');
        end
        function [a,b] = calcZ(~,X)
            TF = X(1);
            syms z;
            for n = 2:length(X)
                TF = TF + X(n)*z^-(n-1);
            end
            [a,b]=numden(TF);
            a =  sym2poly(a);
            b = sym2poly(b);
        end
        function t = calcTime(obj,X)
            t = 0:obj.T:obj.T*(size(X,2)-1);
        end
    end
end