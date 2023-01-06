import MQLib

MQLib.test(; examples=true) do model
    MQLib.set_heuristic(model, first(MQLib.heuristics()))
end