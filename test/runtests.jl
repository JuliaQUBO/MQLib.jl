import Test
import TOML
import MQLib
import MQLib: MOI, QUBODrivers

function first_available_tool(names::Vector{String})
    for name in names
        path = Sys.which(name)
        if !isnothing(path)
            return path
        end
    end

    return nothing
end

Test.@testset "Compatibility metadata" begin
    root = dirname(dirname(@__FILE__))
    project = TOML.parsefile(joinpath(root, "Project.toml"))
    compat = project["compat"]

    Test.@test compat["julia"] == "1.10"
    Test.@test compat["QUBODrivers"] == "0.4, 0.5"
    Test.@test compat["QUBOTools"] == "0.12"

    ci = read(joinpath(root, ".github", "workflows", "ci.yml"), String)
    Test.@test occursin(r"version:\s*'1\.10'", ci)
    Test.@test occursin(r"version:\s*'1'", ci)
end

Test.@testset "QUBODrivers" begin
    QUBODrivers.test(MQLib.Optimizer) do model
        MOI.set(model, MOI.Silent(), true)
        MOI.set(model, MQLib.Heuristic(), first(MQLib.heuristics()))
    end
end

Test.@testset "C ABI contract" begin
    root = dirname(dirname(@__FILE__))
    include_dir = joinpath(root, "c_api", "include")
    header = joinpath(include_dir, "mqlib_c_api.h")
    source = joinpath(root, "c_api", "src", "mqlib_c_api.cpp")

    Test.@test isfile(header)
    Test.@test isfile(source)

    recipe = joinpath(root, "jll", "build_tarballs.jl")
    Test.@test isfile(recipe)

    header_text = read(header, String)
    for symbol in (
        "MQLIB_C_ABI_VERSION",
        "MQLibCQUBOInput",
        "MQLibCQUBOResult",
        "mqlib_solve_qubo",
        "mqlib_c_status_message",
        "MQLIB_STATUS_HYPERHEURISTIC_DATA_NOT_FOUND",
    )
        Test.@test occursin(symbol, header_text)
    end

    cc = first_available_tool(["cc", "gcc", "clang"])
    if isnothing(cc)
        @info "Skipping C ABI header syntax check because no C compiler was found"
    else
        mktempdir() do dir
            check = joinpath(dir, "check_mqlib_c_api.c")
            write(
                check,
                """
                #include "mqlib_c_api.h"

                int main(void) {
                    MQLibCQUBOInput input;
                    MQLibCQUBOResult result;
                    (void)input;
                    (void)result;
                    return mqlib_c_abi_version() == MQLIB_C_ABI_VERSION ? 0 : 0;
                }
                """,
            )
            Test.@test success(`$cc -std=c99 -I$include_dir -fsyntax-only $check`)
        end
    end

    recipe_text = read(recipe, String)
    for snippet in (
        "ExecutableProduct(\"MQLib\", :MQLib)",
        "LibraryProduct(\"libmqlib_c_api\", :libmqlib_c_api)",
        "c_api/include/mqlib_c_api.h",
        "c_api/src/mqlib_c_api.cpp",
        "! -name main.cpp",
        "-DMQLIB_C_BUILD_SHARED",
    )
        Test.@test occursin(snippet, recipe_text)
    end

    upstream_dir = get(ENV, "MQLIB_UPSTREAM_DIR", "")
    cxx = first_available_tool(["c++", "g++", "clang++"])
    if isempty(upstream_dir) || !isdir(joinpath(upstream_dir, "include"))
        @info "Skipping C ABI source syntax check because MQLIB_UPSTREAM_DIR is not set"
    elseif isnothing(cxx)
        @info "Skipping C ABI source syntax check because no C++ compiler was found"
    else
        upstream_include = joinpath(upstream_dir, "include")
        Test.@test success(
            `$cxx -std=c++11 -I$include_dir -I$upstream_include -fsyntax-only $source`,
        )
    end
end
