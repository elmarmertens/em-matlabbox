function [iwishDraws, SigmaT, dof, cholSigmaT] = bayesVCVgibbsDraw1(Sigma0T, dof0, resid, rndStream, diagFlag)
% BAYESVARRESIDVCV iwishDraws variance-covariance matrix
% [iwishDraws, SigmaT, dof, cholSigmaT] = bayesVCVgibbsDraw1(Sigma0T, dof0, resid, rndStream, diagFlag)
% Prior is (Sigma0T, dof0) invWishart and posterior is (SigmaT, dof) invWishart
% Recall: Mean of inverse Wishart is SigmaT / (dof - N - 1)
% dof = dof0 + T
%
% resid is T x N
% 
% Optional Argument:
%     diagFlag: if true it will assume diagonal VCV and treat input/output of Sigma0T, 
%               SigmaT and iwishDraws as column vectors (default is false)
%
% Optional Output: wishDraws and icholSigmaT, if inverse draws are desired
%
% See also iwishdraw

%   Coded by  Elmar Mertens, em@elmarmertens.com


%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 04-Mar-2009 15:56:44 $
% $Revision : 1.00 $
% DEVELOPED : 7.7.0.471 (R2008b)
% FILENAME  : bayesVARresidVCV.m

%% parse inputs
if nargin < 4 || isempty(rndStream)
    rndStream = getDefaultStream;
end
if nargin < 5
    diagFlag = false;
end

[T, Ny] = size(resid);

if ~isscalar(dof0)
   error('dof must be a scalar')
end

%% draw random normals
dof         = dof0 + T;
z           = randn(rndStream, Ny, dof);

if ~diagFlag
   
   %% Posterior Update
   SigmaT      = Sigma0T + resid' * resid;
   cholSigmaT  = chol(SigmaT)';
   
   %    if doInvDraws
   %       icholSigmaT   = eye(Ny) / cholSigmaT;
   %       wishDraws     = zeros(Ny, Ny, Ndraw);
   %    end

   
   %% compute invWishart iwishDraws from iW(Sigma, dof)
	iwishDraws	= cholSigmaT / (z * z') * cholSigmaT';
   
else % diagonal case -- assumes vector input/output instead of VCV-matrix
   
   %% Posterior Update
   SigmaT = Sigma0T(:) + sum(resid.^2)'; % or: diag(resid' * resid)
   
   %% construct inverse-gammas
   zz           = sum(z.^2,2);
   % iwishDraws   = bsxfun(@rdivide, SigmaT, zz);
   iwishDraws   = SigmaT ./ zz;
end
