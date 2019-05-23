classdef FaultDetector < matlab.mixin.Copyable
    % This Fault Detector is an object for solving general
    %   computations required by a CUSUM or GLR 
    properties
        g   % g function
        m0  % Normal operation mean
        m1  % Fault operation mean
        M   % Window length
        v   % Variance
        vf  % Faulty variance (for EM)
        h   % Threshold
        rv  % Residual vector (for GLR)
        G   % Gain of sum of squares (for GLR)
        md  % Mean difference
        ma  % Mean average
        mds % Sum of mean differences (for GLR)
        tau % Responsibility time constant (for EM)
        Ts  % Sampling time (for EM)
    end
    methods
        function [g,fault] = GLR(fd,z,h)
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
        end
        function [g,fault] = CUSUM(fd,z,h)
            % CUSUM algorithm, for both scalar and vector case
            if nargin > 2
                fd.h = h;
            end
            fd.g = max(0,fd.g + fd.md'/fd.v*(z - fd.ma));
            g = fd.g;
            fault = g > fd.h;
        end
        function [g,fault] = EM(fd,z,h)
            % CUSUM algorithm, for both scalar and vector case
            if nargin > 2
                fd.h = h;
            end
            gain0 = (1-fd.g+0.01)/sqrt(2*pi*det(fd.v));
            gain1 = (fd.g+0.01)/sqrt(2*pi*det(fd.vf));
            r0 = gain0*exp(-0.5*(z-fd.m0)'/fd.v*(z-fd.m0));
            r1 = gain1*exp(-0.5*(z-fd.m1)'/fd.vf*(z-fd.m1));
            fd.g = (1-fd.Ts/fd.tau)*fd.g + fd.Ts/fd.tau*r1/(r0 + r1); % pi1
            g = fd.g;
            fault = g > fd.h;
        end
        function initialize(fd,Mean,Variance,Method,FurtherParameters)
            fd.m0 = Mean.m0;
            fd.v = Variance;
            switch Method
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
                case 'EM'
                    D = length(fd.m0);
                    fd.g = 0.05;
                    fd.vf = power(1/FurtherParameters.MeanDensityRatio-1,2/D)*fd.v;
                    fd.tau = FurtherParameters.ResponsibilityTimeConstant;
                    fd.Ts = FurtherParameters.SamplingTime;
                    fd.m1 = fd.m0;
                    fd.h = 0.5;
                otherwise
                    warning('Method required is unknown')
            end
            if isfield(FurtherParameters,'Threshold')
                fd.h = FurtherParameters.Threshold;
            else
                fd.thresholdInitialization(Method,FurtherParameters);
            end
        end
        function thresholdInitialization(fd,Method,FurtherParameters)
            % Let M be window size
            % lambda is the non-centrality parameter in the non central
            % Chi square distribution, see (7.39).
            % The GLRT has the test statistic given by (7.35)
            % where z(i) is Normal with mean mu_1 and variance sigma^2.
            %
            % g(k) = 1/(2*sigma^2*M)* (sum(i=k-M+1,k)(z(i)-mu_0))^2
            % the noncentrality parameter for the distribution of 2*g_M is (7.39):
            switch Method
                case 'GLR'
                case 'CUSUM'
            end

%             % range for g(k)
% %             g = 0:0.1:15;
%             %Under H0:
%             dof = 1;
%             pdfH0 = chi2pdf(2*g,dof);
%             cdfH0 = chi2cdf(2*g,dof);
%             
%             %Under H1:
%             dof = 1;
%             lambda = fd.M *(mu_1-mu_0)^2/sigma^2;
%             pdfH1 = ncx2pdf(2*g,dof,lambda);
%             cdfH1 = ncx2cdf(2*g,dof,lambda);
%             
%             % choose threshold h (threshold for 2*g is 2*h !!)
%             
%             h = threshold;
%             % compute PF and PD:
%             PF = 1 - chi2cdf(2*h, dof);
%             PD = 1 - ncx2cdf(2*h, dof,lambda);
%             
%             nullh  = [0:0.1:2*h];
%             PFline = (1-PF)*ones(size(nullh));
%             PDline = (1-PD)*ones(size(nullh));
%             
%             hrangep = [0:0.1:0.3];
%             hlinep  = 2*h*ones(size(hrangep));
%             hrangec = [0:0.1:1];
%             hlinec  = 2*h*ones(size(hrangec));
%             
%             % Returns
%             Pmissed = 1 - PD;
%             Pfalse = PF;
        end
    end
end