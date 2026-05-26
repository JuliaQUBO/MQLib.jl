import MQLib
import MQLib: MOI, QUBODrivers
import Test

const QUBOTools = MQLib.QUBOTools

function test_model()
    return QUBOTools.Model{Int,Float64,Int}(
        Set(1:2),
        Dict(1 => 1.0, 2 => -1.0),
        Dict((1, 2) => 2.0);
        scale = 1.0,
        offset = 0.0,
        sense = :max,
        domain = :bool,
    )
end

Test.@testset "MQLib input" begin
    model = test_model()
    io = IOBuffer()

    QUBOTools.write_model(io, model, QUBOTools.Format{:qubo}(; style = :mqlib))

    expected = String(take!(io))

    mktempdir() do dir
        file_path = joinpath(dir, "model.qubo")
        input = MQLib._mqlib_input(model, file_path)

        if MQLib._supports_mqlib_stdin()
            Test.@test input isa MQLib._MQLibStdinInput
            Test.@test MQLib._mqlib_file_path(input) == MQLib._MQLIB_STDIN
            Test.@test String(copy(input.data)) == expected
            Test.@test !isfile(file_path)
        else
            Test.@test input isa MQLib._MQLibFileInput
            Test.@test MQLib._mqlib_file_path(input) == file_path
            Test.@test read(file_path, String) == expected
        end
    end
end

QUBODrivers.test(MQLib.Optimizer) do model
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MQLib.Heuristic(), first(MQLib.heuristics()))
end
