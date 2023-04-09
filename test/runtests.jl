import MQLib: MOI, MQLib, QUBODrivers

QUBODrivers.test(MQLib.Optimizer) do model
    MOI.set(model, MOI.Silent(), true)
    MQLib.set_heuristic(model, first(MQLib.heuristics()))
end