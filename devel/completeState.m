function state=completeState(DEs,PositionConstraints,VelocityConstraints,partialState,options)

    arguments
        DEs
        PositionConstraints
        VelocityConstraints
        partialState
        options.hints=[];
        options.doDiagnostics=0;
    end


    state=helpers.timeStepDAESystem(DEs,PositionConstraints,VelocityConstraints,[0,0],partialState,hints=options.hints,doDiagnostics=options.doDiagnostics);

end

