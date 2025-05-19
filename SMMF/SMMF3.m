classdef SMMF3 < PROBLEM
% <2025> <multi> <real> <multimodal>
% SMMF3
% Multi-modal Multi-objective test Function with shifted local Pareto front

%------------------------------- Reference --------------------------------
% Block optimization and switchable hybrid clustering for multimodal
% multiobjective evolutionary optimization with shifted local Pareto front
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = public)
        POS;    % Pareto optimal set for IGDX calculation
    end
    methods
        %% Initialization
        function Setting(obj)
            obj.M = 2;
            obj.D = 2;
            [obj.lower, obj.upper] = obj.ParameterSet([-20, -20], [20, 20]);
            obj.encoding = ones(1,obj.D);            
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            x = PopDec;

            a = 1;
            b = 10;
            c = 8;
            for i=1:size(PopDec,1)
                temp_t1 = sign(x(i,1))*ceil((abs(x(i,1))-(a+c/2))/(2*a+c));
                temp_t2 = sign(x(i,2))*ceil((abs(x(i,2))-b/2)/b);
                t1 = sign(temp_t1)*min(abs(temp_t1),1);
                t2 = sign(temp_t2)*min(abs(temp_t2),1);
                x1 = x(i,1)-t1*(c+2*a);
                x2 = x(i,2)-t2*b;
                y = fun([x1,x2],a);
                PopObj(i,:) = y;
            end
            
            ref_point = [-10,10;-10,0;-10,-10;0,10;0,0;0,-10;10,10;10,0;10,-10;];
            dist1 = pdist2(PopDec,ref_point);
            [~,idx1] = min(dist1,[],2);
            PopObj(idx1 == 1 | idx1 == 2 | idx1 == 3,:) = PopObj(idx1 == 1 | idx1 == 2 | idx1 == 3,:) + 1;
            PopObj(idx1 == 4 | idx1 == 5 | idx1 == 6,:) = PopObj(idx1 == 4 | idx1 == 5 | idx1 == 6,:) + 0;
            PopObj(idx1 == 7 | idx1 == 8 | idx1 == 9,:) = PopObj(idx1 == 7 | idx1 == 8 | idx1 == 9,:) + 1;
        end
        %% Sample reference points on Pareto front
        function P = GetOptimum(obj,N)
            temp = load(strrep(mfilename, '_modified', ''));
            obj.POS = temp.PS;
            P = temp.PF;
            P = [P;P+1];
        end
    end
end

function y = fun(x,a)
    y = zeros(2,1);
    y(1) = (x(1)+a)^2+x(2)^2;
    y(2) = (x(1)-a)^2+x(2)^2;
end