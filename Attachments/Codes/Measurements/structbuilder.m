function [nstruc,losse,losst]=structbuilder(nmax,Ze,Zt)
% This function is to build all the structures up to model order nmax. In
% the input, e denotes estimation set, t denotes test set

% -- Determine models --
nstruc=[]; losse=[]; losst=[];
for ia=1:nmax,
    for ib=1:nmax+1,
        for ic=1:nmax,
            ns=[ia ib ic 0 0 1];
            disp(ns)
            nstruc=[nstruc; ns];
            th=pem(Ze,ns);
            res=pe(Ze,th);
            losse=[losse; th.EstimationInfo.LossFcn];% res'*res/length(res)]; % Loss fcn for estimation data
            res=pe(Zt,th);
            losst=[losst; res'*res/length(res)]; % Loss fcn for test data
        end % ic
    end % ib
end % ia

