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
    options.AbsTol=[];
    options.RelTol=[];
end

soln=helpers.timeStepDAESystem(DEs,PositionConstraints,VelocityConstraints,tspan,ics,doPlot=options.doPlot,...
    doDiagnostics=options.doDiagnostics,solverChoice=options.solverChoice,hints=options.hints,AbsTol=options.AbsTol, RelTol=options.RelTol);




end