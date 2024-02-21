function [obj] = fun_eigstraininfl(obj,chooseall)
%un_eigstraininfl computes eigenstress influence tensors Q for all Phases in the RVE object
% NOT TESTED, written MArch2019

if ~exist('chooseall','var')
    chooseall=0;
end
fld=fieldnames(obj.phases);
if chooseall
    fldpi=fld;
else
    %check whether phase eigenstress is nonzero
    for i=1:numel(fld)
        fldpi={};
        if ~isempty(obj.phases.(fld{i}).pi)
            fldpi{end+1}=fld{i}; %#ok<AGROW>
        end
    end
end
I4 = fun_I4();

% check whether all Hill tensors are identical
testhill=0;
for it=1:numel(fld)
    if ~ismember(obj.phases.(fld{it}).shape.def,{'sphere','matrix'})
        testhill=testhill+1;
    end
end

% compute Q tensors
if testhill==0;
    % if all Hill tensors are equal then simplified definition, Eq.27 of Pichler&Hellmich 2010
    for it=1:numel(fld)
        for jt=1:numel(fldpi)
            Q.(fld{it}).(fldpi{jt})=((it==ij)*I4.I-obj.A.(fld{it}))*obj.f.(fldpi{jt})*obj.Ainf.(fldpi{jt})*obj.P.(fldpi{jt});
        end
    end
else
    %calcuate new Hill tensors for all isotropic phases for which no stroud points were used before
    for it=1:numel(fld)
        if strcmp(obj.phases.(fld{it}).shape.dis,'iso') && isempty(obj.phases.(fld{it}).APstroud.P)
            
            [obj]=fun_Ainf_4iso_Q(obj.phases.(fld{i}),obj);
        end
    end
    
    % calculate sumfAinfP
    checkisomat=fun_check_isotropy(obj.C0);
    sumfAinfP=zeros(6);
    for it=1:numel(fld)
        if strcmp(obj.phases.(fld{it}).shape.dis,'single')
            sumfAinfP=sumfAinfP+obj.f.(fld{it})*obj.Ainf.(fld{it})*obj.P.(fld{it});
        elseif strcmp(obj.phases.(fld{it}).shape.dis,'iso')
            if strcmp(checkisomat,'iso')
                sumfAinfP=sumfAinfP+obj.f.(fld{it})*fun_isoav(obj.APstroud.(fld{it}).Ainfe3*obj.APstroud.(fld{it}).Pe3);
            else
                error('not programmed yet')
            end
        else
        end
    end
    
    if strcmp(obj.scheme,'SCS')
        % if SCS scheme, simplified formula, Eq.51 of Pichler&Hellmich 2010
        for it=1:numel(fld)
            for jt=1:numel(fldpi)
                Q.(fld{it}).(fldpi{jt})=((it==ij)*I4.I-obj.A.(fld{it}))*obj.f.(fldpi{jt})*obj.Ainf.(fldpi{jt})*obj.P.(fldpi{jt})...
                    +(obj.A.(fld{it})*sumfAinfP-obj.Ainf.(fld{it})*obj.P.(fld{it}))*obj.f.(fldpi{jt})*transpose(obj.A.(fldpi{jt}));
            end
        end
    else
        % if other scheme, generalized formula, Eq.25 of Pichler&Hellmich 2010
        % calculate sumfAinfP
        sumfCCAinfP=zeros(6);
        for it=1:numel(fld)
            if strcmp(obj.phases.(fld{it}).shape.dis,'single')
                sumfCCAinfP=sumfCCAinfP+obj.f.(fld{it})*(obj.Chom-obj.phases.(fld{it}).mat.C)*obj.Ainf.(fld{it})*obj.P.(fld{it});
            elseif strcmp(obj.phases.(fld{it}).shape.dis,'iso')
                if strcmp(checkisomat,'iso')
                    sumfCCAinfP=sumfCCAinfP+obj.f.(fld{it})*fun_isoav((obj.Chom-obj.phases.(fld{it}).mat.C)*obj.APstroud.(fld{it}).Ainfe3*obj.APstroud.(fld{it}).Pe3);
                else
                    error('not programmed yet')
                end
            else
            end
        end
        
        for it=1:numel(fld)
            for jt=1:numel(fldpi)
                Q.(fld{it}).(fldpi{jt})=((it==ij)*I4.I-obj.A.(fld{it}))*obj.f.(fldpi{jt})*obj.Ainf.(fldpi{jt})*obj.P.(fldpi{jt})...
                    +(obj.A.(fld{it})*sumfAinfP-obj.Ainf.(fld{it})*obj.P.(fld{it}))*sumfCCAinfP*obj.f.(fldpi{jt})...
                    *(transpose(I4.I-obj.A.(fldpi{jt}))+(obj.Chom-obj.phases.(fldpi{jt}).mat.C)*obj.Ainf.(fldpi{jt})*obj.P.(fldpi{jt}));
            end
        end
    end
end
obj.Q=Q;
end

