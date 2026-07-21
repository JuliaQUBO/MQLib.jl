import Test
import TOML
import MQLib
import MQLib: MOI, QUBODrivers

const QUBOTools = MQLib.QUBOTools

_compat_entries(value::AbstractString) = strip.(split(value, ','))

function first_available_tool(names::Vector{String})
    for name in names
        path = Sys.which(name)
        if !isnothing(path)
            return path
        end
    end

    return nothing
end

function configure_public_c_api_smoke!(model)
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MQLib.Heuristic(), "ALKHAMIS1998")
    MOI.set(model, MQLib.RandomSeed(), 1234)
    MOI.set(model, MQLib.NumberOfReads(), 2)
    MOI.set(model, MOI.TimeLimitSec(), 0.02)

    return model
end

function configure_public_default_hyperheuristic_smoke!(model)
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MQLib.RandomSeed(), 1234)
    MOI.set(model, MQLib.NumberOfReads(), 1)
    MOI.set(model, MOI.TimeLimitSec(), 0.02)

    return model
end

function solution_metadata(model)
    raw = MOI.get(model, MOI.RawSolver())

    return QUBOTools.metadata(QUBOTools.solution(raw))
end

function test_public_default_hyperheuristic_succeeds()
    T = Float64
    n = 3
    model = MOI.instantiate(MQLib.Optimizer; with_bridge_type = T)
    x, _ = MOI.add_constrained_variables(model, fill(MOI.ZeroOne(), n))

    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{T}}(),
        MOI.ScalarQuadraticFunction{T}(
            [
                MOI.ScalarQuadraticTerm{T}(-6, x[1], x[2]),
                MOI.ScalarQuadraticTerm{T}(-1, x[2], x[3]),
            ],
            [
                MOI.ScalarAffineTerm{T}(5, x[1]),
                MOI.ScalarAffineTerm{T}(3, x[2]),
                MOI.ScalarAffineTerm{T}(1, x[3]),
            ],
            zero(T),
        ),
    )
    configure_public_default_hyperheuristic_smoke!(model)

    MOI.optimize!(model)

    Test.@test MOI.get(model, MOI.ResultCount()) > 0
    Test.@test length(MOI.get.(model, MOI.VariablePrimal(), x)) == n

    read_count_attr = QUBODrivers.QUBOTools_MOI.NumberOfReads
    total_reads = sum(
        MOI.get(model, read_count_attr(result_index))
        for result_index = 1:MOI.get(model, MOI.ResultCount())
    )
    Test.@test total_reads == 1

    return nothing
end

function test_public_c_api_bool_max_objectives()
    T = Float64
    n = 3
    Q = T[-1 2 2; 2 -1 2; 2 2 -1]
    model = MOI.instantiate(MQLib.Optimizer; with_bridge_type = T)
    x, _ = MOI.add_constrained_variables(model, fill(MOI.ZeroOne(), n))

    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{T}}(),
        MOI.ScalarQuadraticFunction{T}(
            [
                MOI.ScalarQuadraticTerm{T}(Q[i, j], x[i], x[j])
                for i = 1:n for j = 1:n if i != j
            ],
            [MOI.ScalarAffineTerm{T}(Q[i, i], x[i]) for i = 1:n],
            T(3),
        ),
    )
    configure_public_c_api_smoke!(model)

    MOI.optimize!(model)

    Test.@test MOI.get(model, MOI.ResultCount()) > 0
    for result_index = 1:MOI.get(model, MOI.ResultCount())
        xi = MOI.get.(model, MOI.VariablePrimal(result_index), x)
        expected = sum(Q[i, j] * xi[i] * xi[j] for i = 1:n for j = 1:n) + T(3)
        Test.@test MOI.get(model, MOI.ObjectiveValue(result_index)) ≈ expected
    end

    read_count_attr = QUBODrivers.QUBOTools_MOI.NumberOfReads
    total_reads = sum(
        MOI.get(model, read_count_attr(result_index))
        for result_index = 1:MOI.get(model, MOI.ResultCount())
    )
    Test.@test total_reads == 2

    metadata = solution_metadata(model)
    Test.@test metadata["seeds"]["sampler"] == 1234

    Test.@test isempty(QUBODrivers.validate_metadata(metadata))
    Test.@test metadata["origin"] == "MQLib.jl"
    Test.@test metadata["algorithm"]["name"] == "ALKHAMIS1998"
    Test.@test metadata["backend"]["name"] == "MQLib"
    Test.@test metadata["backend"]["version"] == MQLib.__VERSION__
    Test.@test metadata["reads"]["number_of_reads"] == 2
    Test.@test metadata["reads"]["final_number_of_reads"] == 2
    Test.@test metadata["time"]["effective"] > 0.0

    return nothing
end

function test_public_c_api_spin_max_objectives()
    T = Float64
    n = 3
    h = T[-1; -1; -1]
    J = T[0 4 4; 0 0 4; 0 0 0]
    model = MOI.instantiate(MQLib.Optimizer; with_bridge_type = T)
    s, _ = MOI.add_constrained_variables(model, fill(QUBODrivers.Spin(), n))

    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{T}}(),
        MOI.ScalarQuadraticFunction{T}(
            [
                MOI.ScalarQuadraticTerm{T}(J[i, j], s[i], s[j])
                for i = 1:n for j = 1:n
            ],
            [MOI.ScalarAffineTerm{T}(h[i], s[i]) for i = 1:n],
            T(5),
        ),
    )
    configure_public_c_api_smoke!(model)

    MOI.optimize!(model)

    Test.@test MOI.get(model, MOI.ResultCount()) > 0
    for result_index = 1:MOI.get(model, MOI.ResultCount())
        si = MOI.get.(model, MOI.VariablePrimal(result_index), s)
        expected =
            sum(J[i, j] * si[i] * si[j] for i = 1:n for j = 1:n) +
            sum(h[i] * si[i] for i = 1:n) +
            T(5)
        Test.@test MOI.get(model, MOI.ObjectiveValue(result_index)) ≈ expected
    end

    return nothing
end

Test.@testset "Compatibility metadata" begin
    root = dirname(dirname(@__FILE__))
    project = TOML.parsefile(joinpath(root, "Project.toml"))
    compat = project["compat"]

    Test.@test compat["julia"] == "1.10"
    Test.@test "0.1.2" in _compat_entries(compat["MQLib_jll"])
    Test.@test "0.6.1 - 0.6" in _compat_entries(compat["QUBODrivers"])
    Test.@test all(
        version -> version in _compat_entries(compat["QUBOTools"]),
        ("0.13", "0.14", "0.15", "0.16"),
    )

    ci = read(joinpath(root, ".github", "workflows", "ci.yml"), String)
    Test.@test occursin(r"version:\s*'1\.10'", ci)
    Test.@test occursin(r"version:\s*'1'", ci)
end

Test.@testset "QUBODrivers" begin
    Test.@test QUBODrivers.supports_seed(MQLib.Optimizer)
    Test.@test QUBODrivers.honors_final_reads(MQLib.Optimizer)
    Test.@test QUBODrivers.enforces_time_limit(MQLib.Optimizer)

    QUBODrivers.test(MQLib.Optimizer; benchmark_conformance = true) do model
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

    Test.@test MQLib._mqlib_has_c_api()
    Test.@test MQLib._mqlib_can_use_c_api("ALKHAMIS1998")
    if !MQLib._mqlib_has_hyperheuristic_data()
        Test.@test !MQLib._mqlib_can_use_c_api(nothing)
    end

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

    test_public_c_api_bool_max_objectives()
    test_public_c_api_spin_max_objectives()
    test_public_default_hyperheuristic_succeeds()
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
