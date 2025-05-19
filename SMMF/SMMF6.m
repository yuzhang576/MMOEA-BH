classdef SMMF6 < PROBLEM
% <2025> <multi> <real> <multimodal>
% SMMF6
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
            [obj.lower,obj.upper] = obj.ParameterSet([0, 0], [1, 2]);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopDec_old = PopDec;
            PopObj(:,1) = PopDec(:,1);
            PopDec(:,2) = PopDec(:,2)-(PopDec(:,2)>1);
            y2 = PopDec(:,2)-PopDec(:,1).^0.5;
            PopObj(:,2) = 1-sqrt(PopDec(:,1))+2.*((4.*y2.^2)-2*cos(20.*y2.*pi/sqrt(2))+2);            
            
            ref_point = [0.1,0.3; 0.43,0.66; 0.9,0.9; 0.1,1.3; 0.43,1.66; 0.9,1.9; ];
            dist1 = pdist2(PopDec_old,ref_point);
            [~,idx1] = min(dist1,[],2);
            PopObj(idx1 == 1 | idx1 == 2 | idx1 == 3,:) = PopObj(idx1 == 1 | idx1 == 2 | idx1 == 3,:) + 0;
            PopObj(idx1 == 4 | idx1 == 5 | idx1 == 6,:) = PopObj(idx1 == 4 | idx1 == 5 | idx1 == 6,:) + 2;
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