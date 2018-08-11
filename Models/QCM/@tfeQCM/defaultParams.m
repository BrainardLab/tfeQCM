function [params,paramsLb,paramsUb] = defaultParams(obj,varargin)
% [params,paramsLb,paramsUb] = defaultParams(obj,varargin)
%
% Set objects params to default values as well as provide reasonable lower
% and upper bournds.
%
% All three returns are in struct form, use paramsToVec on the structs to
% get vector form.

%% Default parameters
% Quadratic parameters
%
% We only store 5 because we handle the amplitude of the response
% in the amplitude parameter below.  The first axis of the ellipse
% has an implicit value of 1.
params.Qvec = [1 1 1 0 0];

% Let's have a Naka-Rushton sigmoidal contrast response function
params.crfAmp = 1;
params.crfSemi = 1;
params.crfExponent = 2;
params.offset = 0; 
% Exponential falloff
params.expFalloff = 0.3;

% Noise level (used for simulations
params.noiseSd = 0.2;

%% Lower bounds
paramsLb.Qvec = [1e-3 1e-3 0 0 0];
paramsLb.crfAmp = 1e-3;
paramsLb.crfSemi = 1e-3;
paramsLb.crfExponent = 1e-2;
paramsLb.expFalloff = 1e-1;
paramsLb.noiseLevel = 0;
paramsLb.offset = -2;
%% Upper bounds
paramsUb.Qvec = [1e3 1e3 2*pi 2*pi 2*pi];
paramsUb.crfAmp = 1e3;
paramsUb.crfSemi = 1e3 ;
paramsUb.crfExponent = 1e2;
paramsUb.expFalloff = 1e1;
paramsUb.noiseLevel = 100;
paramsUb.offset = 2; 
end