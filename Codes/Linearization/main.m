% Comment or uncomment the lines
clearvars
close all

%% Model creation, nonlinear model, linear model
% modelcreation
modelcreation_simplified
%% Substitution and normalization
% substitution
substitution_simplified
% normalization
%% Modal, controllability, observability analysis with plots
analysis
%% Creation of Kalman filter, which saves the filter into kf.mat and produces a transient plot
kfinit
%% Check the surface plot of partial differentials along a pressure-enthalpy range
% differentials
