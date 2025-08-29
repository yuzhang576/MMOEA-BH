
function [new_new_up_s, new_new_low_s] = clip_upper_lower(new_up_s, new_low_s, center)
% Create logical array indicating if each point's dimensions are within bounds
inside_each_dim = (center >= new_low_s) & (center <= new_up_s);

% Use all function to check if all dimensions of each point are within bounds
inside_rectangle = all(inside_each_dim, 2);

% Get points inside the rectangle
points_inside = center(inside_rectangle, :);
num_inside = sum(inside_rectangle); % Number of points inside the rectangle

if num_inside < 2
    new_new_up_s = new_up_s;
    new_new_low_s = new_low_s;
else
    % Calculate original hypervolume (area)
    original_volume = prod(new_up_s - new_low_s);
    
    % Target hypervolume is 1/a of original hypervolume
    target_volume = original_volume / num_inside;
    
    % Perform linear regression fitting
    X = points_inside(:, 1:end-1);
    y = points_inside(:, end);
    mdl = fitlm(X, y);
    
    % Get linear regression coefficients (excluding constant term)
    coefficients = round(1./abs([mdl.Coefficients.Estimate(2:end)', 1]), 2);
    coefficients(isinf(coefficients)) = 10^6;

    %%    
    % Original side lengths (different dimensions)
    original_dimensions = new_up_s - new_low_s;
   
    % Define symbolic variable a
    syms a
    
    % Construct volume equation
    equation_terms = 1;
    for i = 1:length(original_dimensions)
        equation_terms = equation_terms * (original_dimensions(i) - coefficients(i) * a);
    end
    
    % Set volume equation equal to target volume
    equation = equation_terms == target_volume;
    
    % Solve the equation
    solutions = double(solve(equation, a));
    
    % Select valid solutions (greater than 0 and less than original length divided by coefficient)
    valid_solution = [];
    for i = 1:length(solutions)
        if all(solutions(i) > 0) && all(original_dimensions ./ coefficients > solutions(i))
            valid_solution = solutions(i);
            break;
        end
    end
    
    if isempty(valid_solution)
        error('No valid solution found');
    end
    
    % Calculate new side lengths
    adjusted_target_side_length = original_dimensions - coefficients * valid_solution;

    %%
    % Calculate centroid
    centroid = mean([new_up_s + new_low_s] / 2, 1);
    
    % Calculate new upper and lower bounds
    new_new_low_s = centroid - adjusted_target_side_length / 2;
    new_new_up_s = centroid + adjusted_target_side_length / 2;
    
    % Ensure new bounds remain valid
    new_new_up_s = max(new_new_up_s, new_low_s);
    new_new_low_s = min(new_new_low_s, new_up_s);
end
end