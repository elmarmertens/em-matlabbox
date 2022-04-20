function stars = Zstarchi2(Z, dof)

%   Coded by  Elmar Mertens, em@elmarmertens.com
%   Extended by Gergely Ganics gergely.ganics@gmail.com


    if  Z > chi2inv(0.99, dof)
        stars = '^{\ast\ast\ast}';
    elseif Z > chi2inv(0.95, dof)
        stars = '^{\ast\ast}';
    elseif Z > chi2inv(0.9, dof)
        stars = '^{\ast}';
    else
        stars = '';
    end

end