function [Kvstat,Kpval,CVMvstat,CvMpval,KS95cv,T,rvec] = PITtest(PITvec)
%% PIT test of Rossi and Sekhposyan (2019)

% get rid of NaNs - input can contain NaNs but the valid observations must be contiguous
PITvec = PITvec(~isnan(PITvec));
T      = size(PITvec, 1); % PITvec is (T x 1)
el     = floor(T^(1/3)); % block length, must be o(T^0.5), RS used T^1/3 and T^1/4 in MC
bootMC = 1000;
rvec   = 0 : 0.001 : 1; % rvec is (1 x numr)
% numr   = size(rvec,2);

cumcumz = (PITvec < rvec) - rvec;
v       = sum(cumcumz,1)/sqrt(T);

Kvstat   = max(abs(v));
CVMvstat = mean((v.^2));
pgrid    = 0.01 : 0.01 :0.99;
tablecv  = bootstrapInoue(el,bootMC,PITvec,rvec,pgrid);
% 95th percentile of the bootstrap distribution of the Kolmogorov-Smirnov statistic, for plotting
KS95cv=tablecv(1,95);
% now look up the p-value
% KS stat
Kpval=1-interp1(tablecv(1,:)',pgrid',Kvstat,'linear','extrap');
if Kpval < 0
    Kpval = 0;
elseif Kpval > 1
    Kpval = 1;
end
% CvM stat
CvMpval=1-interp1(tablecv(2,:)',pgrid',CVMvstat,'linear','extrap');
if CvMpval < 0
    CvMpval = 0;
elseif CvMpval > 1
    CvMpval = 1;
end

end % function PITtest

function result = bootstrapInoue(el, bootMC, pit, rvec, pgrid)
%% bootstrapping PIT stats
KSv  = zeros(bootMC,1);
CVMv = zeros(bootMC,1);

P         = size(pit,1);
invsqrtP  = 1 ./ sqrt(P);

emp_cdf        = pit <= rvec;
emp_cdf        = emp_cdf - mean(emp_cdf,1);
emp_cdf_sumel  = NaN(size(emp_cdf));
for j = 1 : P-el+1
    emp_cdf_sumel(j,:) = sum(emp_cdf(j:j+el-1,:), 1);
end
bigNormMat = 1/sqrt(el) * randn(P-el+1,bootMC);

parfor bootrep = 1 : bootMC

    z      = bigNormMat(:, bootrep);
    K_star = zeros(1,length(rvec));
    
    for j = 1 : P-el+1
         K_star = K_star + invsqrtP .* z(j,1) .* emp_cdf_sumel(j,:); %#ok<PFBNS>
    end

    KSv(bootrep,1) =  max(abs(K_star'));
    CVMv(bootrep,1) = mean(K_star'.^2);
end

ndx  = round(bootMC*pgrid);
KSv  = sort(KSv,'ascend');
cvKv = KSv(ndx);
CVMv = sort(CVMv,'ascend');
cvMv = CVMv(ndx);

result = [cvKv'; cvMv'];
end % bootstrapInoue