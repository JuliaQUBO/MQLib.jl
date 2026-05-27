#include <stdio.h>

#include "mqlib_c_api.h"

static int solve_with_heuristic(
    const char *heuristic,
    const char *hhdata_dir,
    double runtime_limit_seconds
) {
    const double linear[3] = {5.0, 3.0, 1.0};
    const int32_t quadratic_i[2] = {0, 1};
    const int32_t quadratic_j[2] = {1, 2};
    const double quadratic_value[2] = {-6.0, -1.0};

    int32_t solution[3] = {0, 0, 0};
    char selected_heuristic[64];
    double history_values[32];
    double history_times[32];

    MQLibCQUBOInput input = {
        MQLIB_C_ABI_VERSION,
        3,
        linear,
        2,
        quadratic_i,
        quadratic_j,
        quadratic_value,
        MQLIB_C_INDEX_BASE_ZERO,
        heuristic,
        runtime_limit_seconds,
        1234,
        hhdata_dir
    };

    MQLibCQUBOResult result = {
        MQLIB_C_ABI_VERSION,
        0.0,
        0.0,
        solution,
        3,
        selected_heuristic,
        (int32_t)sizeof(selected_heuristic),
        history_values,
        history_times,
        32,
        0
    };

    const int status = mqlib_solve_qubo(&input, &result);
    if (status != MQLIB_STATUS_OK) {
        fprintf(stderr, "MQLib failed: %s\n", mqlib_c_status_message(status));
        return status;
    }

    printf(
        "%s objective %.15g solution %d %d %d runtime %.6f seconds\n",
        result.selected_heuristic,
        result.objective_value,
        result.solution[0],
        result.solution[1],
        result.solution[2],
        result.runtime_seconds
    );
    return MQLIB_STATUS_OK;
}

int main(int argc, char **argv) {
    int status = solve_with_heuristic("ALKHAMIS1998", NULL, 0.01);
    if (status != MQLIB_STATUS_OK) {
        return status;
    }

    if (argc < 2) {
        fprintf(stderr, "Pass the path to MQLib's hhdata directory to run the hyperheuristic example.\n");
        return 2;
    }

    return solve_with_heuristic(NULL, argv[1], 0.75);
}
