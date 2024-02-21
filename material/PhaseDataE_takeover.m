classdef PhaseDataE_takeover < PhaseData
    properties
        eigenstress % eigenstresss (global base)
    end
    methods
        % CONSTRUCTOR
        function newphaseE = PhaseDataE(phasename,matobj,phaseshape,phasevol,phaseigenstress)
            % Call PhaseData constructor
            newphaseE@PhaseData(phasename,matobj,phaseshape,phasevol); 
            if nargin>0
                % Add eigenstress
                newphaseE.eigenstress=phaseigenstress;
            end
        end
    end
end