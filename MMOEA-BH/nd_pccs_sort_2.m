 %%
function [selected_population, sorted_population, sorted_SCD]= nd_pccs_sort_2(population,aSize)
    % Modified the criteria for combining x and y
    % Special pccs on x and y, same as SCD

    try
        aSize;
    catch
        aSize = 1000;
    end
    PopDec = population.decs;
    PopObj = population.objs;
    [N,~] = size(PopDec);
    if N == 1
        selected_population = population;
        sorted_population = population;
        sorted_SCD = 1;
    else        
        CD_x = map0to1( -calc_pccs(PopDec) );
        CD_f = map0to1( -calc_pccs(PopObj) );
        % For cases where x is large y is small, x is small y is large, x is large y is large - sum x and y
        % For case where x is small y is small - take the minimum of x and y
        idx_max = bitor(CD_x>mean(CD_x), CD_f>mean(CD_f));
        idx_min = ~idx_max;
        CD(idx_max) = CD_x(idx_max) + CD_f(idx_max);
        CD(idx_min) = min(CD_x(idx_min), CD_f(idx_min));
        
        [sorted_SCD, idx] = sort(CD,'descend');
        selected_population = population(idx(1));
        if numel(idx)>aSize
            sorted_population = population(idx(1:aSize));
            sorted_SCD = sorted_SCD(1:aSize);
        else
            sorted_population = population(idx);
        end
    end
end