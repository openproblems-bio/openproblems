import json
import openproblems

_CORE_TEST_SUITES = ["test_core"]
_TASK_TEST_SUITES = ["test_task"]


def generate_matrix():
    suites = _CORE_TEST_SUITES.copy()
    for task in openproblems.TASKS:
        task_name = openproblems.utils.get_member_id(task)
        suites.extend([f"{suite} and {task_name}" for suite in _TASK_TEST_SUITES])
    return suites


def main():
    matrix = generate_matrix()
    print(json.dumps(matrix))


if __name__ == "__main__":
    main()
