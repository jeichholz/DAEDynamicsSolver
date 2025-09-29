function soln=timeStepODESystem(ODESystem,tspan,InitialConditions,options)
%function soln=timeStepODESystem(ODESystem,vars,InitialConditions,tspan,PostProcessFunctions)
%
%Take the system of ODE defined in ODESystem and plug it in to ODE23s.
%Optionally, evaluate the quantities of interest defined in
%PostProcessFunctions. 
%
%INPUT
%
%ODESystem -- a system of pure ODEs, no DAEs.  Should be given as a vector
%of symbolic functions/expressions.  Not a cell. 
%
%InitialConditions -- A cell array of initial conditions.  Initial conditions can be
%specified in a number of different ways.  Examples:
%   x==1  <-- set symbolic function x to a value
%   x(t) == 1 <-- set expression x(t) to a value
%   diff(x) == 1 <-- set derivative of symbolic function x to 1. 
%   "Dxt==1" <-- set derivative of x to 1.  The naming convention is
%   "D<function name><derivative order>.  So "Dthetatt==2" means the second
%   derivative of theta with respect to t. 
%
%tspan -- a vector containing the time span to integrate over. 
%
%PostProcessingFunctions -- optional.  If none, omit or set to [].  If
%included, a list of explicit equations detailing how to calculate
%quantities of interest from the state functions.  For instance, in a
%single pendulum, if the state variables are theta and diff(theta), one
%post-processing function might be 
%x==sin(theta). 
%
%doPlot -- optional.  If omitted or set to [] then defaults to 1.  If the
%value is 1 then a plot of all calculated quantities is calculated. 
%
%%doDiagnositics -- Optional, with default value 1.  If true, then print
%diagnostic information as you proceed. 
%
%OUTPUT
%
%soln -- a struct with one field for each quantity of interest. If
%postprocessing equations are included, those fields are included in soln. 

arguments
ODESystem
tspan
InitialConditions
options.PostProcessFunctions=[];
options.doPlot=1;
options.doDiagnostics=0;
end

doPlot=options.doPlot;
doDiagnostics=options.doDiagnostics;
PostProcessFunctions=options.PostProcessFunctions;


    if ~exist("doPlot","var") || isempty(doPlot)
        doPlot=1;
    end

    if ~exist("doDiagnostics","var") || isempty(doDiagnostics)
        doDiagnostics=1;
    end

    ODESystem=evalSymfun2Expr(ODESystem);
    InitialConditions=evalSymfun2Expr(InitialConditions);

    ICstruct=parseConditions(InitialConditions);

    [firstOrderSystem,S]=odeToVectorField(ODESystem);

    [~,~,v2f,~,baseFunctions]=symFunsToSymVars(ODESystem);

    vars=unique(v2f(baseFunctions.values));

    dfprintf("The unknown functions (just the base functions, not derivatives) are:\n%s",vars(1));
    for i=2:numel(vars)
        dfprintf(",\n%s",vars(i));
    end
    dfprintf("\n\n")

    %Create a map of state variables to indices. 
    varToIdx=createVarToIndexMap(vars,S);
    
    %Create initial condition vector. 
    Y0=structToVec(ICstruct,varToIdx);

    %Now solve. 
    [t,ys]=ode23s(matlabFunction(firstOrderSystem,'vars',{'t','Y'}),tspan,Y0);

    %Put it in to a structure. 
    soln=struct();
    soln.t=t;
    for field=reshape(fieldnames(ICstruct),1,[])
        field=field{1};
        soln=setfield(soln,field,ys(:,varToIdx(field)));
    end

    %Massivly ugly way to evaluate all the post-processing functions 
    %using the data in the solved ODEs.  
    ICstructnames=fieldnames(ICstruct);


    %Build a string to evaluate the function RH with the fields in soln
    %plugged in. 
    evalRHString=['RH(soln.t,soln.',ICstructnames{1}];
    for j=2:numel(ICstructnames)
        evalRHString=[evalRHString,',soln.',ICstructnames{j}];
    end
    evalRHString=[evalRHString,');'];

    %Now, unlike all other second derivatives, all second derivatives of
    %the state variables won't be included in the PP equations. So if our
    %state variables are theta, theta', phi, phi', there are no PP
    %equations for theta'' and phi'', because those are in the governing
    %DEs.  Nevertheless, Simon says we should add that information to the
    %solution structure, because acceleration is frequently of interest to
    %engineers. So regardless of whether there are PP equations or not,
    %fill in these higher derivatives using the DEs. 
    
    DEs=symFunsToSymVars(ODESystem);
    for i=1:numel(DEs)
        LH=lhs(DEs(i));
        RH=matlabFunction(rhs(DEs(i)),'vars',{'t',ICstructnames{:}});
        val=eval(evalRHString);
        soln=setfield(soln,string(LH),val);
    end



    %Create a new RH for each different PP equation, and evaluate. 
    if exist("PostProcessFunctions","var") && ~isempty(PostProcessFunctions)
        PP=symFunsToSymVars(PostProcessFunctions);
        for i=1:numel(PP)
            LH=lhs(PP(i));
            RH=matlabFunction(rhs(PP(i)),'vars',{'t',ICstructnames{:}});
            val=eval(evalRHString);
            soln=setfield(soln,string(LH),val);
        end
    end

    %Plotting!
    if doPlot
        %Plot yo-self!
        figure();
        hold on
        for i=1:numel(ICstructnames)
            plot(soln.t,getfield(soln,ICstructnames{i}),'LineWidth',2)
        end
        legend(ICstructnames,'Location','Best')

        grid on
        hold off
    end


    function dfprintf(varargin)
        if doDiagnostics
            fprintf(varargin{:})
        end
    end


    function symToIdx=createVarToIndexMap(vars,S)
        [~,~,varsToFuns,~,baseVars]=symFunsToSymVars(vars);
        baseVars=unique(baseVars.values);
        symToIdx=dictionary();
        %Convert between naming conventions :(
        for ii=1:numel(baseVars)
            v=baseVars(ii);
            symToIdx(string(v))=find(S==v);
            differentiationOrder=1;
            while (1)
                vv=sym(['D',char(v),repmat('t',1,differentiationOrder)]);
                if differentiationOrder==1
                    otherName=sym(['D',char(v)]);
                elseif differentiationOrder>1
                    otherName=sym(['D',num2str(differentiationOrder),char(v)]);
                end
                idx=find(S==otherName);
                if isempty(idx)
                    break;
                else
                    symToIdx(string(vv))=idx;
                    differentiationOrder=differentiationOrder+1;
                end
            end
        end
    end






%Inserted by dependency analyzer
function expr=evalSymfun2Expr(symfun)
    if class(symfun)=="symfun"
        %expr=formula(symfun);
        expr=symfun("t");
    else
        expr=symfun;
    end
end

%Inserted by dependency analyzer
function ICstruct=parseConditions(ics)



    %Sometimes the initial conditions might be interpreted as a
    %vector-valued symfun.  Stop that. 
    if class(ics)=="symfun"
        ics=ics("t");
    end

   
    %It is reasonable to pass in initial conditions as a vector of symbols,
    %but not what we want.  Convert to cell array. 
    if class(ics)=="sym"
        ics=sym2cell(ics);
    end


    %Store the initial conditions in a structure for later
    %substitution/manipulation. 
    ICstruct=struct();

    %Initial conditions may be in the form of a string, "y==4", "Dyt==4"
    %Character array,                                   'y==4', 'Dyt==4'
    %symbolic function value,                           y==4,
    %diff(y,t)==4
    %or symvol value,                                   y(t)==4, diff(y(t),t)==4 
    for i=1:numel(ics)
        if class(ics{i})=="string"
            L=ics{i}.extractBefore("==");
            R=sym(ics{i}.extractAfter("=="));
        elseif class(ics{i})=="char"
            equalSigns=find(ics{i}=="=");
            L=ics{i}(1:equalSigns(1)-1);
            R=sym(ics{i}(equalSigns(2)+1:end));
        else
            [~,asSymbol]=symFunsToSymVars(ics{i});
            tmp=lhs(ics{i});
            L=string(asSymbol(formula(tmp)));
            R=rhs(ics{i});
            R=formula(R);
        end
        ICstruct=setfield(ICstruct,L,R);
    end


end

%Inserted by dependency analyzer
function [vec,isSet]=structToVec(structs,varToIdx)
    if isstruct(structs)
        structs={structs};
    end
    vec=zeros(numel(varToIdx.keys),1);
    isSet=vec;
    for i=1:numel(structs)
        struc=structs{i};
        fields=fieldnames(struc);
        for j=1:numel(fields)
            f=fields{j};
            vec(varToIdx(f))=getfield(struc,f);
            isSet(varToIdx(f))=1;
        end
    end
end



%Inserted by dependency analyzer
function [expr,funs2vars,vars2funs,diffOrder,diffBase]=symFunsToSymVars(expr)

    if isa(expr,"symfun")
        expr=formula(expr);
    end
    
    k=sym("z");
    F=symfun("f(z)",k);
    

    funs2vars=dictionary();
    vars2funs=dictionary();
    diffOrder=dictionary();
    diffBase=dictionary();
    
    expr=mapSymType(expr,"diff",@renamingDiffs);

    expr=mapSymType(expr,"symfun",@renamingSymFuns);

    if funs2vars.isConfigured
        baseFuns=diffBase.values;
        for ff = reshape(baseFuns,1,[])
            if ~ismember(ff,vars2funs.keys)
                g=symfun(str2sym([char(ff),'(t)']),str2sym('t'));
                renamingSymFuns(formula(g));
            else
                g=vars2funs(ff);
            end
            for i=1:diffOrder(ff)
                if ~ismember(diff(g,i),funs2vars.keys)
                    renamingDiffs(formula(diff(g,i)));
                end
            end
        end
    end

    function out=generalRenamer(A)

    end

    function out=renamingDiffs(A)
            %Find the name of the function which we are differentiating. 
            f=findSymType(A,'symfun');
            if numel(f) ~= 1
                error('Error in renaming %s, %d symbolic functions found when one expected.\n',string(A),numel(f));
            end
            %Get the name of the function in a string. 
            fname=symFunType(f);
            %Find the order of the derivative.  Interstingly, you can just
            %count the commas in the string version of A.  Maybe.  
            order=sum(char(string(A))==',');
            
           
            %Now make the string representation. 
            out=['D',char(fname),repmat('t',1,order)];
            sout=sym(out);
            funs2vars(A)=sout;
            vars2funs(sout)=A;
            basefun=sym(fname);
            if ~diffOrder.isConfigured || ~ismember(basefun,diffOrder.keys)
                diffOrder(basefun)=order;
            else
                diffOrder(basefun)=max(diffOrder(basefun),order);
            end
            diffBase(sout)=basefun;
    end

    function out=renamingSymFuns(A)
        
        %Ok, here it gets hard.  Undefined functions are of type symtype,
        %those are easy.  However, if we do diff(f(x(t)),t) we get back a
        %thing like D(f)(x(t)) whose type is also symfun.  :(

        %We might be able to tell the two apart using the class of A, but
        %I'm not sure.  For now, first get the string representaiton of A. 

        str=char(symFunType(A));
        
        %If str contains "D(" then we know we have a derivative.

        if str(1)=='D' && str(2) =='('
            %Then we need to figure out the order of the derivative. Count
            %the number of open parens.
            order=sum(str=='(');

            %Get the name of the thing we are differntiating. Find the
            %first close paren and the last open paren. 
            lo=find(str=='(');
            lo=lo(end);

            fc=find(str==')');
            fc=fc(1);

            %The string between lo and fc is the name of the function being
            %differentiated. 

            fname=str(lo+1:fc-1);

            %Append the appropriate number of p's for prime.
            out=[fname,repmat('p',1,order)];


        else
            %This case is easy.  We let symFunType get the name for us and
            %we are done. 
            out=str;
            order=0;
            fname=str;
        end
        sout=sym(out);
        
        funs2vars(A)=sout;
        vars2funs(sout)=A;
        basefun=sym(fname);
        if ~diffOrder.isConfigured || ~ismember(basefun,diffOrder.keys)
            diffOrder(basefun)=order;
        else
            diffOrder(basefun)=max(diffOrder(basefun),order);
        end
        diffBase(sout)=basefun;



    end

end

    



end
