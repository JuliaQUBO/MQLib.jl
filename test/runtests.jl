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

Test.@testset "Julia C ABI bridge" begin
    linear_terms = Dict(1 => 5.0, 2 => 3.0, 3 => 1.0)
    quadratic_terms = Dict((1, 2) => -6.0, (2, 3) => -1.0)
    linear, quadratic_i, quadratic_j, quadratic_value =
        MQLib._mqlib_problem_data(3, linear_terms, quadratic_terms)

    Test.@test linear == [5.0, 3.0, 1.0]
    Test.@test quadratic_i == Int32[1, 2]
    Test.@test quadratic_j == Int32[2, 3]
    Test.@test quadratic_value == [-6.0, -1.0]
    Test.@test MQLib._mqlib_run_seed(65_536) == 0
    Test.@test endswith(
        MQLib._mqlib_hyperheuristic_data_dir(),
        joinpath("share", "mqlib", "hhdata"),
    )

    if MQLib._mqlib_has_c_api()
        result = MQLib._mqlib_solve_qubo(
            3,
            linear,
            quadratic_i,
            quadratic_j,
            quadratic_value;
            heuristic = "ALKHAMIS1998",
            random_seed = 1234,
            run_time_limit = 0.01,
        )

        Test.@test result.objective_value == 6.0
        Test.@test result.solution == Int32[1, 0, 1]
        Test.@test result.selected_heuristic == "ALKHAMIS1998"
        Test.@test result.runtime_seconds > 0.0
    else
        @info "Skipping direct Julia C ABI solve because MQLib_jll has no libmqlib_c_api product"
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
        "version = v\"0.1.2\"",
        "ExecutableProduct(\"MQLib\", :MQLib)",
        "LibraryProduct(\"libmqlib_c_api\", :libmqlib_c_api)",
        "c_api/include/mqlib_c_api.h",
        "c_api/src/mqlib_c_api.cpp",
        "mkdir -p \"\${bindir}\" \"\${libdir}\" \"\${includedir}\" \"\${datadir}/mqlib/hhdata\"",
        "install -vm 0644 hhdata/*.rf \"\${datadir}/mqlib/hhdata/\"",
        "! -name main.cpp",
        "MQLIB_LIBRARY_SOURCES=\"\$(find src -name '*.cpp' ! -name main.cpp | sort)\"",
        "-DMQLIB_C_BUILD_SHARED",
    )
        Test.@test occursin(snippet, recipe_text)
    end
    Test.@test !occursin("mapfile", recipe_text)

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
