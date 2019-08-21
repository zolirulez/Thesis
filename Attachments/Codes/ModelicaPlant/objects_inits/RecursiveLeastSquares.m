classdef RecursiveLeastSquares < matlab.mixin.Copyable
    % This Recursive Least Squares is an object for solving general
    %   computations required by an RLS algorithm. 
    % Zoltan Mark Pinter, Master Thesis, 2019
    properties
        e   % Residual
        r   % Regressor
        t   % Parameter
        s   % Schur complement
        L   % Learning rate
        P   % Uncertainty
        K   % Gain
    end
    methods
        function t = regression(rls,r,y)
            % Function to be used in matlab script
            rls.e = y - r'*rls.t;
            rls.s = rls.L + r'*rls.P*r;
            rls.K = rls.P*r/rls.s;
            rls.t = rls.t + rls.K*rls.e;
            rls.P = (rls.P - rls.K*rls.s*rls.K')/rls.L;
            t = rls.t;
        end
        function e_t = regression_simulink(rls,r,y,batch)
            % Function to be used in simulink
            rls.e = y - r'*rls.t;
            if ~batch
                rls.s = rls.L + r'*rls.P*r;
            else
                rls.s = 1 + r'*rls.P*r;
            end
            rls.K = rls.P*r/rls.s;
            rls.t = rls.t + rls.K*rls.e;
            % Batch is a binary variable. If 1, the training is more
            % careful
            if ~batch
                P = (rls.P - rls.K*rls.s*rls.K')/rls.L;
            else
                P = (rls.P - rls.K*rls.s*rls.K');
            end
            if ~any(diag(P)<0)
                rls.P = P;
            end
            e_t = [rls.e'; reshape(rls.t,size(rls.t,1)*size(rls.t,2),1)];
        end
        function initialize(rls,ForgettingFactor,Initial)
            % Initialization of object
            rls.t = Initial.t;
            rls.P = Initial.P;
            rls.L = ForgettingFactor;
            rls.e = zeros(1,size(rls.t),2);
        end
    end
end