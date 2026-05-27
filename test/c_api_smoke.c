#include <stdio.h>
#include <string.h>

#include "mqlib_c_api.h"

static void init_input(MQLibCQUBOInput *input, const char *heuristic, const char *hhdata_dir) {
    static const double linear[3] = {5.0, 3.0, 1.0};
    static const int32_t quadratic_i[2] = {0, 1};
    static const int32_t quadratic_j[2] = {1, 2};
    static const double quadratic_value[2] = {-6.0, -1.0};

    input->abi_version = MQLIB_C_ABI_VERSION;
    input->dimension = 3;
    input->linear = linear;
    input->quadratic_count = 2;
    input->quadratic_i = quadratic_i;
    input->quadratic_j = quadratic_j;
    input->quadratic_value = quadratic_value;
    input->index_base = MQLIB_C_INDEX_BASE_ZERO;
    input->heuristic = heuristic;
    input->runtime_limit_seconds = 0.01;
    input->random_seed = 1234;
    input->hyperheuristic_data_dir = hhdata_dir;
}

static void init_result(MQLibCQUBOResult *result, int32_t *solution, char *selected_heuristic) {
    result->abi_version = MQLIB_C_ABI_VERSION;
    result->objective_value = 0.0;
    result->runtime_seconds = 0.0;
    result->solution = solution;
    result->solution_length = 3;
    result->selected_heuristic = selected_heuristic;
    result->selected_heuristic_length = 64;
    result->history_objective_values = NULL;
    result->history_times_seconds = NULL;
    result->history_capacity = 0;
    result->history_length = 0;
}

static int run_success_case(const char *heuristic, const char *hhdata_dir, const char *expected_prefix) {
    MQLibCQUBOInput input;
    MQLibCQUBOResult result;
    int32_t solution[3] = {0, 0, 0};
    char selected_heuristic[64] = {0};

    init_input(&input, heuristic, hhdata_dir);
    init_result(&result, solution, selected_heuristic);

    const int status = mqlib_solve_qubo(&input, &result);
    if (status != MQLIB_STATUS_OK) {
        fprintf(stderr, "expected success, got %s\n", mqlib_c_status_message(status));
        return 1;
    }
    if (strncmp(result.selected_heuristic, expected_prefix, strlen(expected_prefix)) != 0) {
        fprintf(stderr, "unexpected selected heuristic: %s\n", result.selected_heuristic);
        return 1;
    }
    if (result.solution_length != 3) {
        fprintf(stderr, "unexpected solution length: %d\n", result.solution_length);
        return 1;
    }
    return 0;
}

static int run_missing_hhdata_case(void) {
    MQLibCQUBOInput input;
    MQLibCQUBOResult result;
    int32_t solution[3] = {0, 0, 0};
    char selected_heuristic[64] = {0};

    init_input(&input, NULL, NULL);
    init_result(&result, solution, selected_heuristic);

    const int status = mqlib_solve_qubo(&input, &result);
    if (status != MQLIB_STATUS_HYPERHEURISTIC_DATA_NOT_FOUND) {
        fprintf(stderr, "expected missing hyperheuristic data, got %s\n", mqlib_c_status_message(status));
        return 1;
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: %s /path/to/MQLib/hhdata\n", argv[0]);
        return 2;
    }

    if (run_success_case("ALKHAMIS1998", NULL, "ALKHAMIS1998") != 0) {
        return 1;
    }
    if (run_missing_hhdata_case() != 0) {
        return 1;
    }
    if (run_success_case(NULL, argv[1], "HH_") != 0) {
        return 1;
    }

    return 0;
}
