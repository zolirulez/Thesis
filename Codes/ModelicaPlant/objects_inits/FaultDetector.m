classdef FaultDetector < matlab.mixin.Copyable
    % This Fault Detector is an object for solving general
    %   computations required by a CUSUM or GLR 
    properties
        Method % String of the type
        g   % g function
        m0  % Normal operation mean
        m1  % Fault operation mean
        M   % Window length
        v   % Variance
        h   % Threshold
        % Expectation Maximization
        vf  % Faulty variance (for EM)
        tau % Responsibility time constant (for EM)
        Ts  % Sampling time (for EM)
        % Generalized Likelihood Ratio
        rv  % Residual vector 
        G   % Gain of sum of squares
        G1
        mds % Sum of mean differences
        Pfalse      % Probability of false alarm
        Pmissed     % Probability of missed detection
        % CUSUM
        md  % Mean difference
        ma  % Mean average
        Tdetect     % Average time of detection
        Tfalse      % Average time of false alarm
    end
    methods
        function [g,fault] = detect(fd,z,h)
            switch fd.Method
                case 'CUSUM'
                    % CUSUM algorithm, for both scalar and vector case
                    if nargin > 2
                        fd.h = h;
                    end
                    fd.g = max(0,fd.g + fd.md'/fd.v*(z - fd.ma));
                    g = fd.g;
                    fault = g > fd.h;
                case 'GLR'
                    % The function is given by
                    %           g(k) = 1/(2*sigma^2*M)* (sum(i=k-M+1,k)(z(i)-mu_0))^2
                    %
                    
                    if nargin > 2
                        fd.h = h;
                    end
                    fd.mds = fd.mds + z - fd.rv(1);
                    fd.rv = [fd.rv(1,2:end) z];
                    fd.g = fd.G*fd.mds^2;
                    g = fd.g;
                    fault = g > fd.h;
                case 'EM'
                    if nargin > 2
                        fd.h = h;
                    end
                    fd.rv = [fd.rv(:,2:end) z];
                    g = min(max(fd.g,0.01),0.99);
                    gain0 = (1-g)*fd.G;
                    gain1 = fd.g*fd.G1;
                    % Since the matrix is diagonal:
                    r0 = max(gain0*exp(-0.5*sum((fd.rv-fd.m0)'.^2/fd.v,2)),1e-20);
                    r1 = max(gain1*exp(-0.5*sum((fd.rv-fd.m1)'.^2/fd.vf,2)),1e-20);
                    rn1 = r1./(r0+r1);
                    Nk1 = sum(rn1);
                    
                    fd.m1 = 1/Nk1*fd.rv*rn1;
                    fd.g = Nk1/fd.M;%r1/(r0 + r1);
                    g = fd.g;
                    fault = g > fd.h;
                    %fd.g = (1-fd.Ts/fd.tau)*fd.g + fd.Ts/fd.tau*r1/(r0 + r1); % pi1
%   for k = 1:K
%       r(:,k) = pi_k(k)*mvnpdf(data, mu{k}, Sigma{k});
%   end
%   rn = r./sum(r,2);
%   %% Update parameters
%   Nk = sum(rn,1);
%   % XXX: FILL ME IN!
%   mu_new = cell(K,1);
%   Sigma_new = cell(K,1);
%   mu_new(:) = {zeros(1,2)};
%   Sigma_new(:) = {zeros(2,2)};
%   for k = 1:K
%       for n = 1:N
%           mu_new{k} = mu_new{k} + 1/Nk(k)*rn(n,k)*data(n,:);
%       end
%       for n = 1:N
%           Sigma_new{k} = Sigma_new{k} + 1/Nk(k)*rn(n,k)*((data(n,:)-mu_new{k})'*(data(n,:)-mu_new{k}));
%       end
%   end
%   pi_k_new = Nk'/N;
%   mu = mu_new;
%   Sigma = Sigma_new;
%   pi_k = pi_k_new;
                    
            end
        end
        function initialize(fd,Mean,Variance,Method,FurtherParameters)
            fd.Method = Method;
            fd.m0 = Mean.m0;
            fd.v = Variance;
            switch fd.Method
                case 'CUSUM'
                    fd.m1 = Mean.m1;
                    fd.g = 0;
                    fd.ma = (fd.m0 + fd.m1)/2;
                    fd.md = fd.m1 - fd.m0;
                case 'GLR'
                    fd.M = FurtherParameters.WindowLength;
                    fd.G = 1/(2*fd.v*fd.M);
                    fd.rv = fd.m0*ones(1,fd.M);
                    fd.mds = sum(fd.rv - fd.m0);
                    if isfield(Mean,'m1')
                        fd.m1 = Mean.m1;
                    end
                case 'EM'
                    D = length(fd.m0);
                    fd.g = 0.05;
                    fd.vf = power(1/FurtherParameters.MeanDensityRatio-1,2/D)*fd.v;
                    fd.tau = FurtherParameters.ResponsibilityTimeConstant;
                    fd.Ts = FurtherParameters.SamplingTime;
                    fd.m1 = fd.m0;
                    fd.h = 0.5;
                    fd.M = FurtherParameters.WindowLength;
                    fd.G = 1/(2*pi*det(fd.v));
                    fd.G1 = 1/(2*pi*det(fd.vf));
                    fd.rv = fd.m0.*ones(length(fd.v),fd.M);
                otherwise
                    warning('Method required is unknown')
            end
            if isfield(FurtherParameters,'Threshold')
                fd.h = FurtherParameters.Threshold;
            else
                fd.thresholdInitialization(FurtherParameters);
            end
            switch fd.Method
                case 'CUSUM'
                    fd.ARL;
                case 'GLR'
                    fd.GLRdesign;
            end
        end
        function thresholdInitialization(fd,FurtherParameters)
            switch fd.Method
                case 'GLR'
                    fd.Pfalse = FurtherParameters.FalseAlarmProbability;
                    InitialGuess = FurtherParameters.InitialGuess;
                    m1_h = fsolve(@fd.GLRdesign4fsolve,InitialGuess);
                    fd.m1 = m1_h(1);
                    fd.h = m1_h(2);
                case 'CUSUM'
                    fd.Tfalse = FurtherParameters.FalseAlarmTime;
                    InitialGuess = FurtherParameters.InitialGuess;
                    m1_h = fsolve(@fd.ARL4fsolve,InitialGuess);
                    fd.m1 = m1_h(1);
                    fd.h = m1_h(2);
            end
        end
        function resid = ARL4fsolve(fd,m1_h)
            % Input demux
            m1 = m1_h(1);
            h = m1_h(2);
            FalseAlarmTime = fd.Tfalse;
            % ARL function
            sigs = sqrt((m1-fd.m0)^2/fd.v);
            mus = sqrt((m1-fd.m0)^2/(2*fd.v));
            L = @(mus,sigs,h)(sigs^2/2/mus^2*...
                (exp(-2*(mus*h/sigs^2+1.166*mus/sigs))-1+2*(mus*h/sigs^2+1.166*mus/sigs)));
            logTdetect = log(L(mus^2,sigs^2,h));
            logTfalse = L(-mus^2,sigs^2,h);
            % Checking how close we are to the solution
            detectresid = logTdetect;
            falseresid = logTfalse - log(FalseAlarmTime);
            resid = [detectresid; falseresid];
        end
        function [Tdetect,Tfalse] = ARL(fd)
            sigs = sqrt((fd.m1-fd.m0)^2/fd.v);
            mus = sqrt((fd.m1-fd.m0)^2/2/fd.v);
            L = @(mus,sigs,h)(sigs^2/2/mus^2*...
                (exp(-2*(mus*h/sigs^2+1.166*mus/sigs))-1+2*(mus*h/sigs^2+1.166*mus/sigs)));
            fd.Tdetect = L(mus,sigs,fd.h);
            fd.Tfalse = L(-mus,sigs,fd.h);
            Tdetect = fd.Tdetect;
            Tfalse = fd.Tfalse;
        end
        function resid = GLRdesign4fsolve(fd,m1_h)
            % Input demux
            m1 = m1_h(1);
            h = m1_h(2);
            FalseAlarmProbability = fd.Pfalse;
            %Under H0:
            dof = 1;
            %Under H1:
            dof = 1;
            lambda = fd.M *(m1-fd.m0)^2/fd.v;
            % choose threshold h (threshold for 2*g is 2*h !!)
            % compute PF and PD:
            PF = 1 - chi2cdf(2*h, dof);
            PD = 1 - ncx2cdf(2*h, dof,lambda);
            % Returns
            logPmissed = log(1 - PD);
            logPfalse = log(PF);
            residMissed = logPmissed;
            residFalse = logPfalse - log(FalseAlarmProbability);
            resid = [residMissed; residFalse];
        end
        function [Pmissed, Pfalse] = GLRdesign(fd)
            %Under H0:
            dof = 1;
            %Under H1:
            dof = 1;
            lambda = fd.M *(fd.m1-fd.m0)^2/fd.v;
            % choose threshold h (threshold for 2*g is 2*h !!)
            % compute PF and PD:
            PF = 1 - chi2cdf(2*fd.h, dof);
            PD = 1 - ncx2cdf(2*fd.h, dof,lambda);
            % Returns
            fd.Pmissed = 1 - PD;
            fd.Pfalse = PF;
            Pmissed = fd.Pmissed;
            Pfalse = fd.Pfalse;
        end
    end
end