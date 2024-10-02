import MQLib
import MQLib: MOI, QUBODrivers

QUBODrivers.test(MQLib.Optimizer) do model
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MQLib.Heuristic(), first(MQLib.heuristics()))
end
