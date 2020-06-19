from openproblems import TASKS

def test():
    for task in TASKS:
        for dataset in task.DATASETS:
            X = dataset()
            for method in task.METHODS:
                y = method(X)
                for metric in task.METRICS:
                    m = metric(X, y)


