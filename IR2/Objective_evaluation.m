function [score] = Objective_evaluation(x, Candidate_location, Target_location)
% This function takes in two stl files and provides minimum distance between the points of the two stl files as a score.  
% Rounding the integer variables before sending to geometry creator

x(4) = round(x(4)); x(7) = round(x(7));

try
    Stent_generator(x, Candidate_location);
    full_path = fullfile(Candidate_location,"input.txt");
    writematrix(x, full_path);
    Candidate_stl_location = fullfile(Candidate_location,"Stent.STL");
    Target_stl_location = fullfile(Target_location,"Stent.STL"); 

    % Reading STLs in Matlab, generating pointsdata and evaluating score
    Candidate_data = stlread(Candidate_stl_location);
    Candidate_points = Candidate_data.Points;
    %Candidate_points = Candidate_points - min(Candidate_points);
    Target_data = stlread(Target_stl_location);
    Target_points = Target_data.Points;
    %Target_points = Target_points - min(Target_points);
    dist_matrix = pdist2(Candidate_points, Target_points);
    score = sum(min(dist_matrix,[],1));
catch
    score = 1E6;
end
return