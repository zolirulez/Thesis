classdef RecursiveLeastSquares < matlab.mixin.Copyable
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
        std % Variance of perturbation
    end
    methods
        function t = regression(rls,r,y)
            rls.e = y - r'*rls.t;
            rls.s = rls.L + r'*rls.P*r;
            rls.K = rls.P*r/rls.s;
            rls.t = rls.t + rls.K*rls.e*rls.w;
            rls.P = (rls.P - rls.K*rls.s*rls.K')/rls.L;
            t = rls.t;
        end
        function e_t = regression_simulink(rls,r,y,batch)
            rls.e = y - r'*rls.t;
            rls.s = rls.L + r'*rls.P*r;
            rls.K = rls.P*r/rls.s;
            rls.t = rls.t + rls.K*rls.e*rls.w;
            if ~batch
                rls.P = (rls.P - rls.K*rls.s*rls.K')/rls.L;
            end
            e_t = [rls.e'; reshape(rls.t,size(rls.t,1)*size(rls.t,2),1)];
%             rls.t = rls.t + rls.std*randn(size(rls.t,1),size(rls.t,2)).*rls.t;
%             rls.std = rls.std*0.99;
%             rls.std = max(rls.std,eps);
        end
        function initialize(rls,ForgettingFactor,ResidualWeight,Initial)
            rls.t = Initial.t;
            rls.P = Initial.P;
            rls.L = ForgettingFactor;
            rls.w = ResidualWeight;
            rls.std = 0.1; % For simulink
            rls.e = zeros(1,length(rls.w));
        end
    end
end