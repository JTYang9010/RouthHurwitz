% Optimize the RouthHurwithz algorithm
function [numRHP, isStable] = routhOptimized(coeff)

    tol = 1e-10;

    hasZeroRoot = (coeff(end) == 0);
    % 初始化前先分離原點根
    while numel(coeff) > 1 && coeff(end) == 0
        coeff(end) = [];
    end

    n = numel(coeff) - 1;

    %======================
    if n == 0
        numRHP = 0;
        isStable = ~hasZeroRoot;
        return;
    end
    %======================

    numCols = ceil((n+1)/2);

    prev2 = zeros(1, numCols);
    prev1 = zeros(1, numCols);

    oddCoeff = coeff(1:2:end);
    evenCoeff = coeff(2:2:end);

    prev2(1:numel(oddCoeff)) = oddCoeff;
    prev1(1:numel(evenCoeff)) = evenCoeff;

    current = zeros(1, numCols);

    hasZeroRow = false;

    % Handle initial all-zero row
    if all(prev1 == 0)
    
        hasZeroRow = true;
    
        powers = (n-1) - 2*(0:numCols-1);
        powers = max(powers, 0);
    
        prev1 = powers .* prev2;
    
    end
    
    % Handle initial zero pivot
    
    if prev1(1) == 0 && any(prev1 ~= 0)
        prev1(1) = tol;
    end
    
    % Count initial sign change
    numRHP = double(sign(prev2(1)) ~= sign(prev1(1)));
    
    for i = 3:n+1
        
        pivot = prev1(1);

        %  Zero pivot handling
        if pivot == 0
            pivot = tol;
            prev1(1) = pivot;
        end
            
        previous = prev2(1);
    
        current(1:numCols-1) = (pivot * prev2(2:numCols) ...
            - previous * prev1(2:numCols)) / pivot;

        % 先修正 current row 的 zero pivot，再進行 sign-change counting
        % Case 1: All-zero row
        if all(current == 0)
        
            hasZeroRow = true;
        
            % prev1 corresponds to s^(n-i+2)
            order = n - i + 2;
        
            powers = order - 2*(0:numCols-1);
            powers = max(powers, 0);
        
            current = powers .* prev1;
        
        % Case 2: Zero pivot
        elseif current(1) == 0
        
            current(1) = tol;
        
        end
        
        % check whether sign change 
        if sign(prev1(1)) ~= sign(current(1))
            numRHP = numRHP + 1;
        end

        % Rolling update/Rotate buffers
        temp = prev2;
        prev2 = prev1;
        prev1 = current;
        current = temp;

        
    end

    isStable = (numRHP == 0) && ~hasZeroRow && ~hasZeroRoot;
    
end