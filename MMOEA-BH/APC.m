%% affinity propagation
function varargout = APC(data)
%% parameter settings
N = size(data,1);
DM = pdist2(data,data);
S = -DM;
 
meds = median(S,'all');
n = 5;

S(logical(eye(N))) = n*meds;

%% update A(i,k) and R(i,k)
A = zeros(N,N);
R = zeros(N,N); % Initialize messages
lam = 0.6; % Set damping factor
count = 0;
for iter = 1:500
    Eold = A + R;
    % Compute responsibilities
    Rold = R;
    AS = A+S;
    [Y,I]=max(AS,[],2);% Get the maximum value of each row as a column vector
    for i=1:N
        AS(i,I(i))=-inf;% Set the maximum value of each row to -inf
    end
    Y2 = max(AS,[],2);
    R = S-Y;
    for i=1:N
        R(i,I(i))=S(i,I(i))-Y2(i);
    end
    R = (1-lam)*R+lam*Rold; % Dampen responsibilities
    % Compute availabilities
    Aold = A;
    Rp = max(R,0);    
    Rp(logical(eye(N))) = R(logical(eye(N)));

    A = sum(Rp,1)-Rp;
    dA = diag(A);
    A = min(A,0);    
    A(logical(eye(N))) = dA;
    A = (1-lam)*A+lam*Aold; % Dampen availabilities

    %terminate the algorithm when these decisions did not change for several iterations.
    E = A + R;
    if diag(Eold) == diag(E)
        count = count + 1;
        if count == 10   
            break;
        end
    else
        count=0;
    end
end
%% outputs
E = R+A; % Pseudomarginals
I = find(diag(E)>0);
K = length(I);% Indices of exemplars
[~, c]=max(E(:,I),[],2);
c(I) = 1:K;                 %% Replace its own position
idx = I(c); % Assignments   %% Find the center points corresponding to N samples
maxCluster = 0;
for i = unique(idx)'
    maxCluster = maxCluster + 1;
    idx(idx==i) = maxCluster;
end
varargout = {idx, maxCluster};
end