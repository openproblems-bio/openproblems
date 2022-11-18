import json
import openproblems

_CORE_TEST_SUITES = ["(test_0_ or test_3_)"]
_TASK_TEST_SUITES = ["test_1_", "test_2_"]


def generate_matrix():
    suites = _CORE_TEST_SUITES.copy()
    for task in openproblems.TASKS:
        task_name = task.__name__.split(".")[-1]
        suites.extend([f"{suite} and {task_name}" for suite in _TASK_TEST_SUITES])
    return suites


def main():
    matrix = generate_matrix()
    print(json.dumps(matrix))


if __name__ == "__main__":
    main()
