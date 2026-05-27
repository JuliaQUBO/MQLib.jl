#ifndef MQLIB_C_API_H
#define MQLIB_C_API_H

#include <stdint.h>

#define MQLIB_C_ABI_VERSION 1u

#define MQLIB_C_INDEX_BASE_ZERO 0
#define MQLIB_C_INDEX_BASE_ONE 1

#if defined(_WIN32) && defined(MQLIB_C_BUILD_SHARED)
#define MQLIB_C_API __declspec(dllexport)
#elif defined(_WIN32) && defined(MQLIB_C_USE_SHARED)
#define MQLIB_C_API __declspec(dllimport)
#else
#define MQLIB_C_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef enum MQLibCStatus {
    MQLIB_STATUS_OK = 0,
    MQLIB_STATUS_ABI_VERSION_MISMATCH = 1,
    MQLIB_STATUS_INVALID_ARGUMENT = 2,
    MQLIB_STATUS_INVALID_HEURISTIC = 3,
    MQLIB_STATUS_BUFFER_TOO_SMALL = 4,
    MQLIB_STATUS_ALLOCATION_FAILED = 5,
    MQLIB_STATUS_INTERNAL_ERROR = 6,
    MQLIB_STATUS_HYPERHEURISTIC_DATA_NOT_FOUND = 7
} MQLibCStatus;

typedef struct MQLibCQUBOInput {
    uint32_t abi_version;

    int32_t dimension;
    const double *linear;

    int64_t quadratic_count;
    const int32_t *quadratic_i;
    const int32_t *quadratic_j;
    const double *quadratic_value;
    int32_t index_base;

    /*
     * NULL or an empty string runs the MQLib hyperheuristic. Otherwise this
     * must be a valid MQLib QUBO or Max-Cut heuristic code.
     */
    const char *heuristic;

    double runtime_limit_seconds;
    int32_t random_seed;

    /*
     * Optional path to the directory containing the hyperheuristic .rf model
     * files. Used only when heuristic is NULL or empty. If this is NULL or
     * empty, mqlib_solve_qubo looks for ./hhdata relative to the caller's
     * current working directory.
     */
    const char *hyperheuristic_data_dir;
} MQLibCQUBOInput;

typedef struct MQLibCQUBOResult {
    uint32_t abi_version;

    double objective_value;
    double runtime_seconds;

    /*
     * Caller-owned buffer with capacity in solution_length on input. The
     * required length is written back on output.
     */
    int32_t *solution;
    int32_t solution_length;

    /*
     * Optional caller-owned buffer. selected_heuristic_length is the byte
     * capacity on input and the required byte count, including NUL, on output.
     */
    char *selected_heuristic;
    int32_t selected_heuristic_length;

    /*
     * Optional caller-owned buffers for objective history. history_capacity is
     * the input capacity of both arrays. history_length is the required entry
     * count on output. If the capacity is too small, the prefix that fits is
     * copied and MQLIB_STATUS_BUFFER_TOO_SMALL is returned.
     */
    double *history_objective_values;
    double *history_times_seconds;
    int32_t history_capacity;
    int32_t history_length;
} MQLibCQUBOResult;

MQLIB_C_API int mqlib_c_abi_version(void);
MQLIB_C_API const char *mqlib_c_status_message(int status);
MQLIB_C_API int mqlib_solve_qubo(
    const MQLibCQUBOInput *input,
    MQLibCQUBOResult *result
);

#ifdef __cplusplus
}
#endif

#endif
