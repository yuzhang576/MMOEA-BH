classdef MMOEA_BH < ALGORITHM
% <2025> <multi> <real> <multimodal>

%------------------------------- Reference --------------------------------
% Block optimization and switchable hybrid clustering for multimodal
% multiobjective evolutionary optimization with shifted local Pareto front
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK GrouProblem. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [n_PBA] = Algorithm.ParameterSet(5);

            %% Generate random population
            Population = Problem.Initialization();
            PBA = cell(1, numel(Population));
            LBA = cell(1, numel(Population));
            for i = 1:numel(Population)
                [~, temp_PBA,temp_PBA_SCD] = nd_pccs_sort_2(Population(i));
                PBA{i} = temp_PBA(1:min(numel(temp_PBA),n_PBA));
                PBA_SCD{i} = temp_PBA_SCD(1:min(numel(temp_PBA),n_PBA));
            end
            [idx_APC, maxCluster] = APC(Population.decs);
            for i = 1:maxCluster
                [~, LBA{i}, LBA_SCD{i}]= nd_pccs_sort_2(Population(i==idx_APC));
            end
            Archive = [LBA{1:maxCluster}];

            [maxCluster_all, Archive_all] = deal([]);
            gen = 1;
            count_cso  = 1;
            %% Optimization
            while Algorithm.NotTerminated(Archive)
                if Problem.FE < Problem.maxFE*0.5 % Stage 1:PSO
                    [Pbest, Nbest] = REPSelection_PSO1(Problem, PBA, PBA_SCD, LBA, LBA_SCD, idx_APC, maxCluster);
                    Population = PSO_SM(Problem,Population,Pbest,Nbest);
                    if mod(gen,10) == 0
                        [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx_APC, maxCluster] = Environmental_Selection_APC_small(Problem, Archive, Population, LBA, PBA, n_PBA);
                    else
                        [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx_APC, maxCluster] = Environmental_Selection_kmeans(Problem, Archive, Population, LBA, PBA, n_PBA, maxCluster);
                    end
                else                              % Stage 2:CSO
                    if mod(count_cso,2) == 1
                        count_cso                 = count_cso + 1;
                        pop_cso                   = [];
                        [orig_vol, center]        = get_orig_vol(Problem, maxCluster, LBA);
                        for i = 1:maxCluster
                            % get_upper_lower
                            pop = LBA{i};
                            [new_up_s, new_low_s] = get_upper_lower(Problem, pop, orig_vol);
                            [new_up_s, new_low_s] = clip_upper_lower(new_up_s, new_low_s, center);
                            % niche_CSO
                            N                     = ceil(Problem.N/maxCluster);
                            FEs                   = N*25;
                            Problem.FE            = Problem.FE + FEs;
                            result                = platemo('algorithm',{@niche_LMOCSO_without_initial,pop},'problem',{eval(['@',class(Problem)]),new_low_s,new_up_s},'N',N,'maxFE',FEs);
                            pop_cso               = [pop_cso,result{end}];
                            temp{i}               = result{end};
                            % t1
                        end
                        if numel(pop_cso) < Problem.N
                            repmat_num = Problem.N - numel(pop_cso);
                            pop_cso    = [pop_cso,pop_cso(randi(numel(pop_cso),repmat_num,1))];
                        end
                        [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx_APC, maxCluster] = Environmental_Selection_DBSCAN(Problem, Archive, pop_cso(1:Problem.N), LBA, PBA, n_PBA, Problem.N, maxCluster);
                    else
                        count_cso = count_cso + 1;
                        [Pbest, Nbest] = REPSelection_PSO1(Problem, PBA, PBA_SCD, LBA, LBA_SCD, idx_APC, maxCluster);
                        Population = PSO_SM(Problem,Population,Pbest,Nbest);                        
                        [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx_APC, maxCluster] = Environmental_Selection_APC_small(Problem, Archive, Population, LBA, PBA, n_PBA);
                    end
                end
                Archive_all = [Archive_all;{Archive}];
                maxCluster_all = [maxCluster_all;maxCluster];
                gen = gen + 1;
            end
        end
    end
end