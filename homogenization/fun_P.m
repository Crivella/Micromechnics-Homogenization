function [Pout] = fun_P(phaseobj,C0,checkiso,tol)
% define string for Hill tensor function
if ismember(phaseobj.shape.def,{'matrix','sphere'}) %SPHERE
    if strcmp(checkiso,'iso')
        Pout=fun_P_sphere_iso(C0);
    elseif strcmp(checkiso,'transiso')
        Pout=fun_P_spheroid_transiso(C0,1,1,tol);
        %the function  fun_P_sphere_transiso contains an error and is
        %therefore not used
    elseif strcmp(checkiso,'aniso')
        Pout=fun_P_ellipsoid_aniso(C0,1,1,tol);
    else
        error(['Hill Tensor for spherical phase ', phaseobj.name,' not calculated'])
    end
elseif ismember(phaseobj.shape.def,{'spheroid'}) % SPHEROID
    if strcmp(checkiso,'iso')
        Pout=fun_P_spheroid_iso(C0,phaseobj.shape.ar,phaseobj.shape.sr);
    elseif strcmp(checkiso,'transiso')
        if phaseobj.shape.sr<1e4 %needle=cyl
            Pout=fun_P_cyl_transiso(C0);
        elseif phaseobj.shape.sr>1e4 %disc
            Pout=fun_P_disc_transiso(C0);
        else
            Pout=fun_P_spheroid_transiso(C0,phaseobj.shape.ar,phaseobj.shape.sr,tol);
            %not programmed yet, needs to be included
        end
    elseif strcmp(checkiso,'aniso')
        Pout=fun_P_ellipsoid_aniso(C0,phaseobj.shape.ar,phaseobj.shape.sr,tol);
    else
        error(['Hill Tensor for spheroidal phase ', phaseobj.name,' not calculated'])
    end
elseif ismember(phaseobj.shape.def,{'ellipsoid'}) % ELLIPSOID
    if ismember(checkiso,{'iso','transiso'})
        Pout=fun_P_ellipsoid_transiso(C0,phaseobj.shape.ar,phaseobj.shape.sr,tol);
    elseif strcmp(checkiso,'aniso')
        Pout=fun_P_ellipsoid_aniso(C0,phaseobj.shape.ar,phaseobj.shape.sr,tol);
    else
        error(['Hill Tensor for ellipsoidal phase ', phaseobj.name,' not calculated']);
    end
end
end

