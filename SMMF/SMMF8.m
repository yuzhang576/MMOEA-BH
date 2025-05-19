classdef SMMF8 < PROBLEM
% <2025> <multi> <real> <multimodal>
% SMMF8
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
            [obj.lower, obj.upper] = obj.ParameterSet([1, -1], [3, 3]);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopDec_old  = PopDec;
            PopDec(:,2) = PopDec(:,2)-2.*(PopDec(:,2)>1);
            PopObj(:,1) = abs(PopDec(:,1)-2); 
            PopObj(:,2) = 1 - sqrt(PopObj(:,1)) + 2.*( PopDec(:,2)-sin(6*pi* abs(PopObj(:,1))+pi)).^2;

            idx1 = PopDec_old(:,2) >= 1.00001;
            idx2 = PopDec_old(:,2) < 1.00001;
            PopObj(idx1,:) = PopObj(idx1,:) + 1.5;
            PopObj(idx2,:) = PopObj(idx2,:) + 0;
        end
        %% Sample reference points on Pareto front
        function P = GetOptimum(obj,N)
            temp = load(strrep(mfilename, '_modified', ''));
            obj.POS = temp.PS;
            P = temp.PF;
            P = [P;P+1.5];
        end
    end
end