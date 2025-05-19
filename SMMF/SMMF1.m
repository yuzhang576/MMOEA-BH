classdef SMMF1 < PROBLEM
% <2025> <multi> <real> <multimodal>
% SMMF1
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
            obj.D = 3;
            [obj.lower, obj.upper] = obj.ParameterSet(zeros(1, obj.D), [4 4 4]);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            for j=1:size(PopDec,1)
                x = PopDec(j,:);
                f = zeros(2,1);
                n = length(x);
                for i=1:n
                    f(1) = f(1)+sin(pi*x(i));
                    f(2) = f(2)+cos(pi*x(i));
                end
                PopObj(j,:) = f;
            end
            ref_point = [       1.2,1.2,1.2;    1.2,3.2,1.2;    1.2,5.2,1.2;    3.2,1.2,1.2;
                3.2,3.2,1.2;    3.2,5.2,1.2;    5.2,1.2,1.2;    5.2,3.2,1.2;    5.2,5.2,1.2;
                1.2,1.2,3.2;    1.2,3.2,3.2;    1.2,5.2,3.2;    3.2,1.2,3.2;    3.2,3.2,3.2;
                3.2,5.2,3.2;    5.2,1.2,3.2;    5.2,3.2,3.2;    5.2,5.2,3.2;];
            dist1 = pdist2(PopDec,ref_point);
            [~,idx1] = min(dist1,[],2);

            PopObj(idx1 <= 9 ,:) = PopObj(idx1 <= 9 ,:) + 2;
            PopObj(idx1 > 9 ,:) = PopObj(idx1 > 9 ,:) + 0;            
        end
        
        %% Sample reference points on Pareto front
        function P = GetOptimum(obj,N)
            temp = load(strrep(mfilename, '_modified', ''));
            obj.POS = temp.PS;
            P = temp.PF;
            P = [P;P+2];
        end
    end
end