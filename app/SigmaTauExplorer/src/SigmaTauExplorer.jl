# The desktop explorer as a package, so it can be compiled into a standalone
# app. The code itself lives in app/StabilityExplorer.jl (which includes
# ExplorerCore.jl); this module only wraps it with the entry points the
# packaging tools expect:
#
#   julia_main()  PackageCompiler.create_app (see build/build.jl)
#   @main         `julia --project=app/SigmaTauExplorer -m SigmaTauExplorer`,
#                 which is also what AppBundler.jl launches
module SigmaTauExplorer

include(joinpath(@__DIR__, "..", "..", "StabilityExplorer.jl"))

# Backend activation is runtime state, so redo it when the module loads from
# a package image or sysimage.
__init__() = GLMakie.activate!()

function julia_main()::Cint
    try
        return explorer_main(ARGS)
    catch
        Base.invokelatest(Base.display_error, current_exceptions())
        return 1
    end
end

function (@main)(args)
    return explorer_main(args)
end

end # module
