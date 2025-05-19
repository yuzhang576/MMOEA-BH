classdef UECLP < PROBLEM
% <2025> <multi> <real>
% Case study on urban emergency center location planning
% lower --- -10 --- Lower bound of decision variables
% upper ---  10 --- Upper bound of decision variables

%------------------------------- Reference --------------------------------
% Block optimization and switchable hybrid clustering for multimodal
% multiobjective evolutionary optimization with shifted local Pareto front
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = public)
        POS;    % Pareto optimal set for IGDX calculation
        Points1 = [-2,-1; -1,-1]; % Target points in lower-left quadrant
        Points2 = [2,1; 1,1];    % Target points in upper-right quadrant
        Points3 = [-2,1; -1,1];  % Target points in upper-left quadrant
        Points4 = [2,-1; 1,-1];  % Target points in lower-right quadrant
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            [lower,upper] = obj.ParameterSet(-10,10);
            obj.M        = 2;
            obj.D        = 2;
            obj.lower    = zeros(1,2) + lower;
            obj.upper    = zeros(1,2) + upper;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            % Initialize result matrix
            PopObj = zeros(size(PopDec,1),2);
            
            % Determine quadrants and create corresponding target point matrices
            idx1 = PopDec(:,1) < 0 & PopDec(:,2) < 0;  % Lower-left quadrant
            idx2 = PopDec(:,1) > 0 & PopDec(:,2) > 0;  % Upper-right quadrant
            idx3 = PopDec(:,1) < 0 & PopDec(:,2) > 0;  % Upper-left quadrant
            idx4 = ~(idx1 | idx2 | idx3);              % Lower-right quadrant
            
            % Batch calculate distances
            if any(idx1)
                PopObj(idx1,:) = pdist2(PopDec(idx1,:),obj.Points1);
            end
            if any(idx2)
                PopObj(idx2,:) = pdist2(PopDec(idx2,:),obj.Points2);
            end
            if any(idx3)
                PopObj(idx3,:) = pdist2(PopDec(idx3,:),obj.Points3);
            end
            if any(idx4)
                PopObj(idx4,:) = pdist2(PopDec(idx4,:),obj.Points4);
            end

            idx1 = PopDec(:,1) >= 0;
            idx2 = PopDec(:,1) < 0;
            PopObj(idx1,:) = PopObj(idx1,:) + 1;
            PopObj(idx2,:) = PopObj(idx2,:) + 0;
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            N = 200;
            % Generate points for each quadrant
            N1 = ceil(N/4); N2 = ceil(N/4); N3 = ceil(N/4); N4 = N - 3*ceil(N/4);
            
            % Lower-left quadrant
            [X1,Y1] = ndgrid(linspace(-2,-1,N1),linspace(-1,-1,N1));
            % Upper-right quadrant
            [X2,Y2] = ndgrid(linspace(1,2,N2),linspace(1,1,N2));
            % Upper-left quadrant
            [X3,Y3] = ndgrid(linspace(-2,-1,N3),linspace(1,1,N3));
            % Lower-right quadrant
            [X4,Y4] = ndgrid(linspace(1,2,N4),linspace(-1,-1,N4));
            
            % Merge all points
            X = [X1(:); X2(:); X3(:); X4(:)];
            Y = [Y1(:); Y2(:); Y3(:); Y4(:)];
            
            % Calculate distance to nearest points
            R = zeros(length(X),2);
            idx1 = X < 0 & Y < 0;
            idx2 = X > 0 & Y > 0;
            idx3 = X < 0 & Y > 0;
            idx4 = ~(idx1 | idx2 | idx3);
            
            if any(idx1)
                R(idx1,:) = pdist2([X(idx1),Y(idx1)],obj.Points1);
            end
            if any(idx2)
                R(idx2,:) = pdist2([X(idx2),Y(idx2)],obj.Points2);
            end
            if any(idx3)
                R(idx3,:) = pdist2([X(idx3),Y(idx3)],obj.Points3);
            end
            if any(idx4)
                R(idx4,:) = pdist2([X(idx4),Y(idx4)],obj.Points4);
            end
            R = [R;R+1];
            
            % Save Pareto optimal set for IGDX calculation
            obj.POS = [X,Y];
        end
        %% Display a population in the decision space
        function DrawDec(obj,Population)
            % Draw target points for four quadrants
            Draw(obj.Points1([1:end,1],:),'-k','LineWidth',1.5);
            Draw(obj.Points2([1:end,1],:),'-k','LineWidth',1.5);
            Draw(obj.Points3([1:end,1],:),'-k','LineWidth',1.5);
            Draw(obj.Points4([1:end,1],:),'-k','LineWidth',1.5);
            
            % Draw target point markers
            Draw(obj.Points1,'o','MarkerSize',6,'Markerfacecolor',[1 1 1]);
            Draw(obj.Points2,'o','MarkerSize',6,'Markerfacecolor',[1 1 1]);
            Draw(obj.Points3,'o','MarkerSize',6,'Markerfacecolor',[1 1 1]);
            Draw(obj.Points4,'o','MarkerSize',6,'Markerfacecolor',[1 1 1]);
            
            % Draw population decision variables
            Draw(Population.decs);
        end
    end
end