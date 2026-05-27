#include "mqlib_c_api.h"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <utility>
#include <vector>

#include "heuristics/heuristic_factory.h"
#include "heuristics/maxcut/max_cut_simple_solution.h"
#include "heuristics/qubo/qubo_simple_solution.h"
#include "metrics/max_cut_metrics.h"
#include "problem/max_cut_instance.h"
#include "problem/qubo_instance.h"
#include "util/random.h"
#include "util/randomForest.h"

namespace {

bool has_text(const char *value) {
    return value != NULL && value[0] != '\0';
}

enum HyperheuristicProblem {
    HYPERHEURISTIC_MAXCUT,
    HYPERHEURISTIC_QUBO
};

struct HyperheuristicChoice {
    bool found;
    HyperheuristicProblem problem;
    std::string code;

    HyperheuristicChoice() :
        found(false),
        problem(HYPERHEURISTIC_MAXCUT),
        code() {}
};

bool parse_double_token(const std::string &text, double *value) {
    char *end = NULL;
    errno = 0;
    const double parsed = std::strtod(text.c_str(), &end);
    if (end == text.c_str() || *end != '\0' || errno == ERANGE) {
        return false;
    }
    *value = parsed;
    return true;
}

bool parse_history(
    const std::string &text,
    std::vector<std::pair<double, double> > *history
) {
    history->clear();
    if (text.size() < 2 || text[0] != '[' || text[text.size() - 1] != ']') {
        return false;
    }

    const std::string body = text.substr(1, text.size() - 2);
    if (body.empty()) {
        return true;
    }

    std::stringstream stream(body);
    std::string item;
    while (std::getline(stream, item, ';')) {
        const std::string::size_type separator = item.find(':');
        if (separator == std::string::npos) {
            return false;
        }

        double objective = 0.0;
        double time = 0.0;
        if (!parse_double_token(item.substr(0, separator), &objective) ||
            !parse_double_token(item.substr(separator + 1), &time)) {
            return false;
        }

        history->push_back(std::make_pair(objective, time));
    }

    return true;
}

std::string model_path(const std::string &data_dir, const std::string &code) {
    const std::string root = data_dir.empty() ? std::string("hhdata") : data_dir;
    return root + "/" + code + ".rf";
}

bool file_exists(const std::string &path) {
    std::ifstream file(path.c_str());
    return file.good();
}

double elapsed_seconds(const struct timeval &start) {
    struct timeval end;
    gettimeofday(&end, 0);
    return (end.tv_sec - start.tv_sec) +
        0.000001 * (end.tv_usec - start.tv_usec);
}

double remaining_runtime_limit(double runtime_limit_seconds, double elapsed) {
    const double remaining = runtime_limit_seconds - elapsed;
    return remaining > 0.0 ? remaining : 0.0;
}

void update_hyperheuristic_choice(
    const std::string &code,
    HyperheuristicProblem problem,
    const std::vector<double> &metrics,
    const std::string &data_dir,
    double *best_probability,
    int *num_best,
    HyperheuristicChoice *choice
) {
    const std::string path = model_path(data_dir, code);
    if (!file_exists(path)) {
        return;
    }

    RandomForest random_forest(path);
    const double probability = random_forest.Predict(metrics);
    if (probability > *best_probability) {
        *best_probability = probability;
        *num_best = 1;
        choice->found = true;
        choice->problem = problem;
        choice->code = code;
    } else if (probability == *best_probability &&
               Random::RandInt(0, *num_best) == *num_best) {
        ++(*num_best);
        choice->found = true;
        choice->problem = problem;
        choice->code = code;
    }
}

int select_hyperheuristic(
    HeuristicFactory *factory,
    const MaxCutInstance &mi,
    const char *hyperheuristic_data_dir,
    HyperheuristicChoice *choice
) {
    GraphMetrics graph_metrics(mi);
    std::vector<double> metrics;
    graph_metrics.AllMetrics(&metrics, NULL);

    const std::string data_dir = has_text(hyperheuristic_data_dir) ?
        std::string(hyperheuristic_data_dir) :
        std::string();
    double best_probability = -1.0;
    int num_best = 1;

    std::vector<std::string> codes;
    factory->MaxCutHeuristicCodes(&codes);
    for (std::vector<std::string>::const_iterator code = codes.begin();
         code != codes.end();
         ++code) {
        update_hyperheuristic_choice(
            *code,
            HYPERHEURISTIC_MAXCUT,
            metrics,
            data_dir,
            &best_probability,
            &num_best,
            choice
        );
    }

    factory->QUBOHeuristicCodes(&codes);
    for (std::vector<std::string>::const_iterator code = codes.begin();
         code != codes.end();
         ++code) {
        update_hyperheuristic_choice(
            *code,
            HYPERHEURISTIC_QUBO,
            metrics,
            data_dir,
            &best_probability,
            &num_best,
            choice
        );
    }

    if (!choice->found) {
        return MQLIB_STATUS_HYPERHEURISTIC_DATA_NOT_FOUND;
    }

    return MQLIB_STATUS_OK;
}

int validate_buffers(
    const MQLibCQUBOInput *input,
    MQLibCQUBOResult *result,
    int32_t *solution_capacity,
    int32_t *selected_capacity
) {
    if (input == NULL || result == NULL) {
        return MQLIB_STATUS_INVALID_ARGUMENT;
    }
    if (input->abi_version != MQLIB_C_ABI_VERSION ||
        result->abi_version != MQLIB_C_ABI_VERSION) {
        return MQLIB_STATUS_ABI_VERSION_MISMATCH;
    }

    *solution_capacity = result->solution_length;
    *selected_capacity = result->selected_heuristic_length;
    result->solution_length = input->dimension;
    result->history_length = 0;

    if (input->dimension <= 0 || input->linear == NULL) {
        return MQLIB_STATUS_INVALID_ARGUMENT;
    }
    if (input->quadratic_count < 0 ||
        input->quadratic_count > std::numeric_limits<int32_t>::max()) {
        return MQLIB_STATUS_INVALID_ARGUMENT;
    }
    if (input->quadratic_count > 0 &&
        (input->quadratic_i == NULL ||
         input->quadratic_j == NULL ||
         input->quadratic_value == NULL)) {
        return MQLIB_STATUS_INVALID_ARGUMENT;
    }
    if (input->index_base != MQLIB_C_INDEX_BASE_ZERO &&
        input->index_base != MQLIB_C_INDEX_BASE_ONE) {
        return MQLIB_STATUS_INVALID_ARGUMENT;
    }
    if (!std::isfinite(input->runtime_limit_seconds) ||
        input->runtime_limit_seconds < 0.0) {
        return MQLIB_STATUS_INVALID_ARGUMENT;
    }
    if (input->random_seed < 0 || input->random_seed > 65535) {
        return MQLIB_STATUS_INVALID_ARGUMENT;
    }

    if (result->solution == NULL ||
        *solution_capacity < input->dimension) {
        return MQLIB_STATUS_BUFFER_TOO_SMALL;
    }
    if (*selected_capacity < 0 ||
        (*selected_capacity > 0 && result->selected_heuristic == NULL)) {
        return MQLIB_STATUS_INVALID_ARGUMENT;
    }
    if (result->history_capacity < 0 ||
        (result->history_capacity > 0 &&
         (result->history_objective_values == NULL ||
          result->history_times_seconds == NULL))) {
        return MQLIB_STATUS_INVALID_ARGUMENT;
    }

    return MQLIB_STATUS_OK;
}

int build_instance_data(
    const MQLibCQUBOInput *input,
    std::vector<double> *linear,
    std::vector<Instance::InstanceTuple> *quadratic
) {
    linear->assign(input->linear, input->linear + input->dimension);
    for (std::vector<double>::const_iterator iter = linear->begin();
         iter != linear->end();
         ++iter) {
        if (!std::isfinite(*iter)) {
            return MQLIB_STATUS_INVALID_ARGUMENT;
        }
    }

    quadratic->clear();
    quadratic->reserve(static_cast<size_t>(input->quadratic_count));

    const int32_t lower = input->index_base;
    const int32_t upper = input->index_base + input->dimension - 1;
    for (int64_t k = 0; k < input->quadratic_count; ++k) {
        const int32_t i = input->quadratic_i[k];
        const int32_t j = input->quadratic_j[k];
        const double value = input->quadratic_value[k];
        if (i < lower || i > upper || j < lower || j > upper ||
            i == j || !std::isfinite(value)) {
            return MQLIB_STATUS_INVALID_ARGUMENT;
        }

        const int first = static_cast<int>(i - input->index_base + 1);
        const int second = static_cast<int>(j - input->index_base + 1);
        quadratic->push_back(
            Instance::InstanceTuple(std::make_pair(first, second), value)
        );
    }

    return MQLIB_STATUS_OK;
}

int copy_selected_heuristic(
    const std::string &selected,
    int32_t selected_capacity,
    MQLibCQUBOResult *result
) {
    const int32_t required =
        static_cast<int32_t>(selected.size()) + static_cast<int32_t>(1);
    result->selected_heuristic_length = required;

    if (selected_capacity == 0 || result->selected_heuristic == NULL) {
        return MQLIB_STATUS_OK;
    }
    if (selected_capacity < required) {
        return MQLIB_STATUS_BUFFER_TOO_SMALL;
    }

    std::memcpy(
        result->selected_heuristic,
        selected.c_str(),
        static_cast<size_t>(required)
    );
    return MQLIB_STATUS_OK;
}

int copy_history(
    const std::vector<std::pair<double, double> > &history,
    MQLibCQUBOResult *result
) {
    result->history_length = static_cast<int32_t>(history.size());

    if (result->history_capacity == 0) {
        return MQLIB_STATUS_OK;
    }

    const int32_t to_copy = std::min(
        result->history_capacity,
        result->history_length
    );
    for (int32_t k = 0; k < to_copy; ++k) {
        result->history_objective_values[k] = history[static_cast<size_t>(k)].first;
        result->history_times_seconds[k] = history[static_cast<size_t>(k)].second;
    }

    if (result->history_capacity < result->history_length) {
        return MQLIB_STATUS_BUFFER_TOO_SMALL;
    }

    return MQLIB_STATUS_OK;
}

int solve_qubo_impl(
    const MQLibCQUBOInput *input,
    MQLibCQUBOResult *result,
    int32_t selected_capacity
) {
    std::vector<double> linear;
    std::vector<Instance::InstanceTuple> quadratic;
    int status = build_instance_data(input, &linear, &quadratic);
    if (status != MQLIB_STATUS_OK) {
        return status;
    }

    std::srand(static_cast<unsigned int>(input->random_seed));

    QUBOInstance qi(quadratic, linear, static_cast<int>(input->dimension));
    HeuristicFactory factory;
    std::unique_ptr<MaxCutInstance> mi;
    std::unique_ptr<QUBOInstance> hyperheuristic_qi;
    std::unique_ptr<MaxCutHeuristic> maxcut_heuristic;
    std::unique_ptr<QUBOHeuristic> qubo_heuristic;
    Heuristic *heuristic = NULL;
    std::string selected;
    const bool validation = false;
    bool qubo_solution_uses_original_instance = true;
    double hyperheuristic_selection_seconds = 0.0;
    double hyperheuristic_runtime_seconds = -1.0;

    if (has_text(input->heuristic)) {
        const std::string requested(input->heuristic);
        if (factory.ValidQUBOHeuristicCode(requested)) {
            qubo_heuristic.reset(factory.RunQUBOHeuristic(
                requested,
                qi,
                input->runtime_limit_seconds,
                validation,
                NULL
            ));
            heuristic = qubo_heuristic.get();
            selected = requested;
        } else if (factory.ValidMaxCutHeuristicCode(requested)) {
            mi.reset(new MaxCutInstance(qi));
            maxcut_heuristic.reset(factory.RunMaxCutHeuristic(
                requested,
                *mi,
                input->runtime_limit_seconds,
                validation,
                NULL
            ));
            heuristic = maxcut_heuristic.get();
            selected = requested;
        } else {
            return MQLIB_STATUS_INVALID_HEURISTIC;
        }
    } else {
        mi.reset(new MaxCutInstance(qi));
        struct timeval hyperheuristic_start;
        gettimeofday(&hyperheuristic_start, 0);

        HyperheuristicChoice choice;
        status = select_hyperheuristic(
            &factory,
            *mi,
            input->hyperheuristic_data_dir,
            &choice
        );
        if (status != MQLIB_STATUS_OK) {
            return status;
        }

        hyperheuristic_selection_seconds =
            elapsed_seconds(hyperheuristic_start);
        const double selected_runtime_limit = remaining_runtime_limit(
            input->runtime_limit_seconds,
            hyperheuristic_selection_seconds
        );

        std::srand(static_cast<unsigned int>(input->random_seed));
        if (choice.problem == HYPERHEURISTIC_MAXCUT) {
            maxcut_heuristic.reset(factory.RunMaxCutHeuristic(
                choice.code,
                *mi,
                selected_runtime_limit,
                validation,
                NULL
            ));
            heuristic = maxcut_heuristic.get();
        } else {
            hyperheuristic_qi.reset(new QUBOInstance(*mi));
            qubo_heuristic.reset(factory.RunQUBOHeuristic(
                choice.code,
                *hyperheuristic_qi,
                selected_runtime_limit,
                validation,
                NULL
            ));
            heuristic = qubo_heuristic.get();
            qubo_solution_uses_original_instance = false;
        }
        selected = "HH_" + choice.code;
        hyperheuristic_runtime_seconds =
            elapsed_seconds(hyperheuristic_start);
    }

    if (heuristic == NULL) {
        return MQLIB_STATUS_INTERNAL_ERROR;
    }

    result->runtime_seconds = hyperheuristic_runtime_seconds >= 0.0 ?
        hyperheuristic_runtime_seconds :
        heuristic->Runtime();

    std::vector<int> assignments;
    if (qubo_heuristic.get() != NULL) {
        const QUBOSimpleSolution &solution = qubo_heuristic->get_best_solution();
        if (qubo_solution_uses_original_instance) {
            assignments = solution.get_assignments();
            result->objective_value = solution.get_weight();
        } else {
            MaxCutSimpleSolution maxcut_solution(solution, *mi, NULL);
            QUBOSimpleSolution qubo_solution(maxcut_solution, qi, NULL);
            assignments = qubo_solution.get_assignments();
            result->objective_value = qubo_solution.get_weight();
        }
    } else {
        QUBOSimpleSolution solution(
            maxcut_heuristic->get_best_solution(),
            qi,
            NULL
        );
        assignments = solution.get_assignments();
        result->objective_value = solution.get_weight();
    }

    if (assignments.size() != static_cast<size_t>(input->dimension)) {
        return MQLIB_STATUS_INTERNAL_ERROR;
    }
    for (int32_t k = 0; k < input->dimension; ++k) {
        result->solution[k] = assignments[static_cast<size_t>(k)];
    }

    std::vector<std::pair<double, double> > history;
    if (!parse_history(heuristic->History(), &history)) {
        return MQLIB_STATUS_INTERNAL_ERROR;
    }
    if (hyperheuristic_runtime_seconds >= 0.0) {
        for (std::vector<std::pair<double, double> >::iterator iter =
                 history.begin();
             iter != history.end();
             ++iter) {
            if (iter != history.begin()) {
                iter->second += hyperheuristic_selection_seconds;
            }
        }
    }

    const int selected_status =
        copy_selected_heuristic(selected, selected_capacity, result);
    const int history_status = copy_history(history, result);
    if (selected_status != MQLIB_STATUS_OK) {
        return selected_status;
    }
    if (history_status != MQLIB_STATUS_OK) {
        return history_status;
    }

    return MQLIB_STATUS_OK;
}

}  // namespace

extern "C" MQLIB_C_API int mqlib_c_abi_version(void) {
    return static_cast<int>(MQLIB_C_ABI_VERSION);
}

extern "C" MQLIB_C_API const char *mqlib_c_status_message(int status) {
    switch (status) {
    case MQLIB_STATUS_OK:
        return "ok";
    case MQLIB_STATUS_ABI_VERSION_MISMATCH:
        return "ABI version mismatch";
    case MQLIB_STATUS_INVALID_ARGUMENT:
        return "invalid argument";
    case MQLIB_STATUS_INVALID_HEURISTIC:
        return "invalid heuristic";
    case MQLIB_STATUS_BUFFER_TOO_SMALL:
        return "buffer too small";
    case MQLIB_STATUS_ALLOCATION_FAILED:
        return "allocation failed";
    case MQLIB_STATUS_INTERNAL_ERROR:
        return "internal error";
    case MQLIB_STATUS_HYPERHEURISTIC_DATA_NOT_FOUND:
        return "hyperheuristic data not found";
    default:
        return "unknown status";
    }
}

extern "C" MQLIB_C_API int mqlib_solve_qubo(
    const MQLibCQUBOInput *input,
    MQLibCQUBOResult *result
) {
    int32_t solution_capacity = 0;
    int32_t selected_capacity = 0;
    int status = validate_buffers(
        input,
        result,
        &solution_capacity,
        &selected_capacity
    );
    if (status != MQLIB_STATUS_OK) {
        return status;
    }

    (void)solution_capacity;

    try {
        return solve_qubo_impl(input, result, selected_capacity);
    } catch (const std::bad_alloc &) {
        return MQLIB_STATUS_ALLOCATION_FAILED;
    } catch (const std::exception &) {
        return MQLIB_STATUS_INTERNAL_ERROR;
    } catch (...) {
        return MQLIB_STATUS_INTERNAL_ERROR;
    }
}
