function state=completeState(DEs,PositionConstraints,VelocityConstraints,partialState,options)

    arguments
        DEs
        PositionConstraints
        VelocityConstraints
        partialState
        options.hints=[];
        options.doDiagnostics=0;
    end


    state=timeStepDAESystem(DEs,PositionConstraints,VelocityConstraints,[0,0],partialState,hints=options.hints,doDiagnostics=options.doDiagnostics);






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
function expr=structsubs(expr,struct)
    fields=fieldnames(struct);
    for i=1:numel(fields)
        expr=subs(expr,fields{i},getfield(struct,fields{i}));
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

    




%Inserted by dependency analyzer
function [soln]=timeStepDAESystem(DEs,PositionConstraints,VelocityConstraints,tspan,ics,options)
%function [soln]=timeStepDAESystem(DEs,PositionConstraints,VelocityConstraints,vars,tspan,ics,hints,doPlot)
%
%Solve a DAE describing some kind of rigid-body motion.  This is purely
%numerical. 
%
%INPUTS: 
%
%DEs -- A system of DEs, and *other equations that don't need to be
%differentiated*.  This would generally include mx''=F and Ig theta''=M
%type equations, as well as other constituatitive relationns like F1=-F3.
%This should be a vector and not, e.g., a cell. 
%
%PositionConstraints -- optional.  If there are none then either omit or
%set to [].  A list of constraints that should be differentiated twice.
%These would tend to include only position-like functions. 
%
%VelocityConstraints -- optional.  If there are none then either omit or
%set to [].  A list of constraints that need to be differentiated once.
%These would tend to include only velocity-like functions. 
%
%tspan -- the time interval to integrate over, in the form [tmin,tmax]. 
%
%ics -- A cell array of initial conditions.  Initial conditions can be
%specified in a number of different ways.  Examples:
%   x==1  <-- set symbolic function x to a value
%   x(t) == 1 <-- set expression x(t) to a value
%   diff(x) == 1 <-- set derivative of symbolic function x to 1. 
%   "Dxt==1" <-- set derivative of x to 1.  The naming convention is
%   "D<function name><derivative order>.  So "Dthetatt==2" means the second
%   derivative of theta with respect to t. 
%
%hints -- Optional, if there are no hints then omit or set to [].  This tool 
%tries to warn you if there are finetely many initial
%states that satisfy your initial conditions.  For example, for a single
%pendulum, if your initial conditions are x==0.5 and diff(x)==0, there are
%two different suitable initial conditions for y.  Use the hints to select
%one of thost initial conditions. The tool will choose the initial
%condition that satisfies the constraints that is closest to your hints.
%Format is identical to ics. 
%
%doPlot -- Optional, if omitted or empty then defaults to 1.  If non-zero then create a plot at the end of the run.  
%
%doDiagnositics -- Options, with default value 1.  If true, then print
%diagnostic information as you proceed. 
%
%solverChoice -- Optional.  Set to one of "ode15i", "ode23t", or "ode15s".
%Default value is ode15s. 


%OUTPUTS
%
%soln -- a structure with a field for every calculated function.
%
%CAVEATS: 
%
% * At the moment the time variable MUST be named t. 
%
% * I intend for it to be ok to pass in either symbolic functions, like x,
% or symbolic expressions like x(t).  If trouble arises, it is safest to
% pass in expressions, like x(t) rather than x, or diff(x(t),t) rather than
% diff(x). 

arguments
    DEs
    PositionConstraints
    VelocityConstraints
    tspan
    ics
    options.hints = []
    options.doPlot = 1
    options.doDiagnostics =0
    options.solverChoice {mustBeMember(options.solverChoice,["ode15s","ode23t","ode15i"])} = "ode15s"
end

    hints=options.hints;
    doPlot=options.doPlot;
    doDiagnostics=options.doDiagnostics;
    solverChoice=options.solverChoice;

    t0=tspan(1);
    tmax=tspan(2);

    %Make sure that nothing is a symfun, but rather a vector of symbolic
    %expressions. 
    DEs=evalSymfun2Expr(DEs);
    PositionConstraints=evalSymfun2Expr(PositionConstraints);
    VelocityConstraints=evalSymfun2Expr(VelocityConstraints);

    ICstruct=parseConditions(ics);

    if exist("hints","var") && ~isempty(hints)
        hintsStruct=parseConditions(hints);
    else
        hintsStruct=[];
    end
    
    if ~exist("solverChoice","var")||isempty(solverChoice)
        solverChoice="ode15s";
    else
        validSolverChoices={'ode15s','ode23t','ode15i'};
        if ~ismember(solverChoice,validSolverChoices)
            error('Solver choice, %s, must be one of: %s.\n',solverChoice,sprintf("%s ",validSolverChoices{:}));
        end
    end


    %Turn the position constraints into two constraints each, and the
    %velocity constraints into one each; 
    constraints=sym([]);
    for i=1:numel(PositionConstraints)
        constraints=[constraints;[PositionConstraints(i); diff(PositionConstraints(i),"t")]];
        DEs(end+1)=diff(PositionConstraints(i),"t","t");
    end

    for i=1:numel(VelocityConstraints)
        constraints(end+1,1)=VelocityConstraints(i);
        DEs(end+1)=diff(VelocityConstraints(i),"t");
    end

    [~,~,varsToSyms,diffOrdersOriginal,baseVars]=symFunsToSymVars([DEs;constraints]);

    constraints=symFunsToSymVars(constraints);
    
    vars=unique(varsToSyms(baseVars.values));

    dfprintf("The unknown functions (just the base functions, not derivatives) are:\n%s",vars(1));
    for i=2:numel(vars)
        dfprintf("\n%s",vars(i));
    end
    dfprintf("\n\n")

    %Reduce the system of DEs to 1st order. 
    [firstOrderSystem,firstOrderFuns] = reduceDifferentialOrder(DEs,vars);

    %The variables in firstOrder are just plain syms, but they have (t) in
    %thier name, which is still an issue for solve.  So rename them to get
    %rid of that. 
    [firstOrderVars,firstOrderFuns2FirstOrderVars]=symFunsToSymVars(firstOrderFuns);
    [~,~,~,diffOrders]=symFunsToSymVars(firstOrderSystem);

    %We'll want to connect those variables to indices in an initial
    %condition soon. 
    firstOrderVarToIdx=dictionary();
    for i = 1:numel(firstOrderVars)
        firstOrderVarToIdx(firstOrderVars(i))=i;
    end

    numConstraints=numel(constraints);
    numSymbols=numel(firstOrderVars);
    numICsGiven=numel(fieldnames(ICstruct));
    numExpectedStates=sum(diffOrders.values)-numConstraints;

    if numICsGiven ~= numExpectedStates
        warning('This system is expected to require %d state variables, but you provided %d initial conditions.\n',numExpectedStates,numICsGiven);
    end

    %You need to find initial conditions that satisfy the constraints that
    %we have.  If there are no constraints, well then just use the initial
    %conditions. 

    if ~isempty(constraints)
        specificConstraints=constraints;

        %ASSUMES THE TIME VARIABLE IS NAMED t
        specificConstraints=subs(specificConstraints,"t",t0);

        %Substitute the given initial conditions into the constraints.
        specificConstraints=structsubs(specificConstraints,ICstruct);

        %Solve for what you can symbolically. It is nice to do it symbolically
        %in order to warn the user if there are multiple feasible choices.
        ICsoln=solve(specificConstraints,'Real',1);

        %It is possible that ICsoln has no solutions in some fields, which is
        %bad, or multiple solutions in some fields, in which case we need to
        %look at the hints to decide or decide arbitrarily.  Deal with that.
        ICsoln=selectIC(ICsoln,hintsStruct);

        %Now we have satisfied the initial constraints of the system, however,
        %that might not have pinned down all of the necessary initial
        %conditions, like unknown forces which will only show up in the DEs
        %themselves. First, build a partial initial condition vector out of the
        %ICs that we were given and that we determined in order to be feasible.
        %Then, use decic to find anything else that we need.
    else
        %If there were no constraints, then the ICsoln should just be an
        %empty struct.
        ICsoln=struct();
    end

    [Y0,Y0Fixed]=structToVec({ICsoln,ICstruct},firstOrderVarToIdx);

    %Use mass-matrix form for solving the ODE. 
    [Msym,Fsym] = massMatrixForm(firstOrderSystem,firstOrderFuns);

    %For some reason, when the inputs are given as functions like x, rather
    %than x(t), it becomes impossible to use the standard odeFunction to
    %turn M,F into matlab functions.  Convert to symVars rather than
    %anything involving names that even look like functions.  
    M=symFunsToSymVars(Msym);
    
    %Check for zero columns in M.  This will indicate free variables in YP0 in the
    %initial condition that we eventually find.  
    isfreeYP0=zeros(1,size(M,2));
    for i=1:size(M,2)
        if all(M(:,i)==0)
            isfreeYP0(i)=1;
        end
    end

    M = matlabFunction(M, "Vars",{sym("t"),firstOrderVars});

    F=symFunsToSymVars(Fsym);
    F =matlabFunction(F, "Vars",{sym("t"),firstOrderVars});

    implicitForm=@(t,y,yp) M(t,y)*yp-F(t,y);


    %Use decic to finish out any as-yet determined initial conditions. 
    [Y0,YP0]=decic(implicitForm,t0,Y0,Y0Fixed,zeros(size(Y0)),zeros(size(Y0Fixed)));


    dfprintf("Initial Conditions:\n")
    dfprintf("Y0/YP0:\n")
    for i=1:numel(Y0)
        if isfreeYP0(i)
            dfprintf("%s: %f\n",string(firstOrderVars(i)),Y0(i));
        else
            dfprintf("%s: %f %s:%f\n",string(firstOrderVars(i)),Y0(i),[char(firstOrderVars(i)),'t'],YP0(i));
        end
    end
    
    %If the timespan contains no time, assume the user is actually just
    %asking for you to complete the initial state. 
    if tspan(2)==tspan(1)
        soln=struct();
        for i=1:numel(firstOrderVars)
            soln=setfield(soln,string(firstOrderVars(i)),Y0(i));
        end
        return;
    end

    %Ok, solve!
    opt = odeset('Mass', M,'RelTol', 10.0^(-7), 'AbsTol' , 10.0^(-7),'InitialSlope',YP0);

    if string(solverChoice)=='ode15i'
        dfprintf('Using ode15i solver.\n')
        solnStruct=ode15i(implicitForm,[t0,tmax],Y0,YP0,opt);
    elseif string(solverChoice)=='ode23t'
        dfprintf('using ode23t solver.\n')
        solnStruct = ode23t(F, [t0, tmax], Y0, opt);
    elseif string(solverChoice)=='ode15s'
        solnStruct= ode15s(F, [t0, tmax], Y0, opt);
    end
    tSol=solnStruct.x';
    ySol=solnStruct.y';

    [~,ySolP]=deval(solnStruct,tSol);

    soln=struct();
    soln.t=tSol;
    for i=1:numel(firstOrderVars)
        soln=setfield(soln,string(firstOrderVars(i)),ySol(:,i));
    end

    %Now, compute which of the base variables is not specified in the
    %solution yet. 
    

    vars=diffOrdersOriginal.keys;
    for i=1:numel(vars)
        x=vars(i);
        order=diffOrdersOriginal(x);
        if (order>0)
            if order==1
                firstOrderVarName=string(x);
            else
                firstOrderVarName=strcat("D",string(x),repmat('t',1,order-1));
            end
            
            soln.(strcat("D",string(x),repmat('t',1,order))) = ySolP(firstOrderVarToIdx(firstOrderVarName),:)';
            
        end
    end

    


    if doPlot
        origVars=numel(vars);
        %Plot yo-self!
        plot(tSol,ySol(:,1:origVars),'LineWidth',2)
    
    
        for k = 1:numel(vars)
            S{k} = char(vars(k));
        end
    
        legend(S, 'Location', 'Best')
        grid on
    end



    function dfprintf(varargin)
        if doDiagnostics
            fprintf(varargin{:});
        end
    end


    function val=selectIC(soln,hintsStruct)

        %First, check to see if there are any empty fields and stop. 
        fields=reshape(fieldnames(soln),1,[]);

        for f=fields
            if isempty(getfield(soln,f{1}))
                error("No suitable initial value for %s found.\n",f{1});
            end
        end

        %See if there is a single unique solution. 
        numentries=zeros(numel(fields),1);
        for ii=1:numel(fields)
            numentries(ii)=numel(getfield(soln,fields{ii}));
        end
        
        %If everyone has just one unique solution, great, go with that.
        if all(numentries==1)
            val=soln;
        else
            if ~(all(numentries==numentries(1)))
                error('Error in selecting initial conditions. Possible initial conditions are %s\n',string(soln));
            end
            dfprintf('Multiple feasible initial conditions found:\n')
            for kk=1:numentries(1)
                for jj=1:numel(fields)
                    tmp=getfield(soln,fields{jj});
                    dfprintf('%s: %f    ',fields{jj},tmp(kk));
                end
                dfprintf('\n')
            end

            if ~isempty(hintsStruct)
                %Make an error vector for each possible initial condition. 
                err=zeros(numentries(1),1);
    
                %Loop over each hint entry.  Calculate the distance from each
                %possible IC to the hint. 
                hintFields=fieldnames(hintsStruct);
                for ll=1:numel(hintFields)
                    hf=hintFields{ll};
                    hfv=getfield(hintsStruct,hf);
                    for jj=1:numentries(1)
                        tmp=getfield(soln,hf);
                        err(jj)=err(jj)+(tmp(jj)-hfv)^2;
                    end
                end
    
                %Which possible IC has the lowest error?
                [~,midx]=min(err);

            else
                midx=1;
            end

            dfprintf("\nSelecting IC:\n")
            %Ok, go with midx. 
            for f=fields
                tmp=getfield(soln,f{1});
                soln=setfield(soln,f{1},tmp(midx));
                dfprintf("%s:%f    ",f{1},tmp(midx));
            end
            if isempty(hintsStruct)
                dfprintf('\nBecause it is first.  Use a hint if you want to select a different initial condition.\n\n')
            else
                dfprintf("\nBecause it is closest to satisfying the hints.\n\n")
            end
            val=soln;

                


        end
        
    end


end
end
