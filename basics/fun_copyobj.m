function b = fun_copyobj(a)
   b = eval(class(a));  %create default object of the same class as a. one valid use of eval
   for pr =  properties(a).'  %copy all public properties
      try   %may fail if property is read-only
         b.(pr) = a.(pr);
      catch
         %warning('failed to copy property: %s', pr);
      end
   end
end
