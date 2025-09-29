function makeDistrib()

    targets=["completeState", "generateGoverningDEs","timeStepDAESystem","timeStepODESystem"];

    supports=["evalSymfun2Expr","parseConditions","structsubs","structToVec","symFunsToSymVars"];


    targetDir="./distrib";

    dependencies=containers.Map;
    for targname=[targets,supports]
        %fprintf("%s:",targname)
        dependencies(targname)=[];
        for supportname=[supports(supports ~= targname),targets(targets ~= targname)]
            if doesDependOn(targname,supportname)
                %fprintf("%s ",supportname)
                dependencies(targname)=[dependencies(targname),supportname]
            end
        end
    end

    for targname=[targets,supports]
        for i=1:(numel(targets)+numel(supports))
            for dep=dependencies(targname)
                dependencies(targname)=unique([dependencies(targname), dependencies(dep)]);
            end
        end
    end

    for targname=targets
        target=readlines(strcat(targname,".m"));
        for j=numel(target):-1:1
            if strcmp(strip(target(j)),"end")
                lastline=j;
                break;
            end
        end
        target(lastline)='';
        target=[target;""];
        for depname=dependencies(targname)
            dep=readlines(strcat(depname,".m"));
            target=[target;"";"%Inserted by dependency analyzer";dep];
        end

        target=[target;"end"];
        writelines(target,strcat(targetDir,"/",targname,".m"))
    end 



            

    function yesno=doesDependOn(source,support)
        %Open the source and read it. 
        sourcelines=readlines(strcat(source,".m"));
        %Strip whitespace
        sourcelines=strip(sourcelines);
        %Begin with a percent sign?
        comments=arrayfun(@(x) x(1)=="%",sourcelines);
        %delete the comment lines
        sourcelines=sourcelines(~comments);
        %Now see if the source lines contain support(
        yesno=any(contains(sourcelines,strcat(support,"(")));
    end
end


