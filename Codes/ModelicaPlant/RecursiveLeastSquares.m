classdef RecursiveLeastSquares < handle
    % This Recursive Least Squares is an object for solving general
    %   computations required by an RLS algorithm. 
    properties
        e   % Residual
        r   % Regressor
        t   % Parameter
        w   % Weights
        s   % Schur complement
        L   % Learning rate
        P   % Uncertainty
        K   % Gain
    end
    methods
        function regression(rls,r,y)
            rls.e = y - r'*rls.t;
            rls.s = rls.L + r'*rls.P*r;
            rls.K = rls.P*r/rls.s;
            rls.t = rls.t + rls.K*rls.e*rls.w;
            rls.P = (rls.P - rls.K*rls.s*rls.K')/rls.L;
        end
        function initialize(rls,ForgettingFactor,ResidualWeight,Initial)
            rls.t = Initial.t;
            rls.P = Initial.P;
            rls.L = ForgettingFactor;
            rls.w = ResidualWeight;
        end
    end
end