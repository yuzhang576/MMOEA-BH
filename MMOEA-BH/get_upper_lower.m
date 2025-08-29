function [new_up_s, new_low_s] = get_upper_lower(Problem, pop, orig_vol)
    pop_dec     = pop.decs;
    upper_small = max(pop_dec,[],1);
    lower_small = min(pop_dec,[],1);
    
    % Calculate the side lengths (or edge lengths) of the original small rectangle
    dimensions_small = upper_small - lower_small;
    % Calculate the ratio of the original small rectangle (ratio between dimensions)
    ratio = dimensions_small ./ max(dimensions_small);
    % Adjust the ratio to ensure each dimension's ratio is not below the minimum threshold
    adjusted_ratio = max(ratio, 0.2);
    % Normalize the adjusted ratio
    ratio = adjusted_ratio ./ sum(adjusted_ratio);
    % The new rectangle's volume (or area) is the same as the original rectangle
    new_volume = orig_vol;
    % Calculate new side lengths (or edge lengths) based on new volume and ratio, maintaining the same proportions
    new_dimensions = (new_volume / prod(ratio))^(1/length(ratio)) * ratio;
    
    % Calculate the center point coordinates of the new rectangle
    center = (upper_small + lower_small) / 2;
    % Calculate new upper-right and lower-left corner coordinates
    new_low_s = center - new_dimensions / 2;
    new_up_s = new_low_s + new_dimensions;
    % Boundary check
    new_up_s = max(min(new_up_s,Problem.upper),Problem.lower);
    new_low_s = max(min(new_low_s,Problem.upper),Problem.lower);
end