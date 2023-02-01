import openproblems
import parameterized
import unittest
import utils.name

TASK_SUMMARY_MINLEN = 40
TASK_SUMMARY_MAXLEN = 80

MODULE = type(openproblems)


@parameterized.parameterized_class(
    ("task",),
    [(task,) for task in openproblems.TASKS],
    class_name_func=utils.name.name_test,
)
class TestTask(unittest.TestCase):
    def test_class(self):
        """Test that task is a module"""
        assert isinstance(self.task, MODULE)

    def test_members(self):
        """Test that task has the required members"""
        assert hasattr(self.task, "_task_name")
        assert isinstance(self.task._task_name, str)
        assert hasattr(self.task, "_task_summary")
        assert isinstance(self.task._task_summary, str)
        assert len(self.task._task_summary) > TASK_SUMMARY_MINLEN
        assert len(self.task._task_summary) < TASK_SUMMARY_MAXLEN
        assert hasattr(self.task, "DEFAULT_LAYER")
        assert isinstance(self.task.DEFAULT_LAYER, str)
        assert self.task.DEFAULT_LAYER in ["counts", "log_normalized", "log_cpm"]
        assert hasattr(self.task, "api")
        assert isinstance(self.task.api, MODULE)
        for list_name in ["DATASETS", "METHODS", "METRICS"]:
            assert hasattr(self.task, list_name)
            method_list = getattr(self.task, list_name)
            assert isinstance(method_list, list)
            assert len(method_list) > 0
            for method in method_list:
                assert callable(method)

    def test_api_members(self):
        """Test that task.api has the required members"""
        assert hasattr(self.task.api, "check_dataset")
        assert hasattr(self.task.api, "check_method")
        assert hasattr(self.task.api, "sample_dataset")
        assert hasattr(self.task.api, "sample_method")
        assert callable(self.task.api.check_dataset)
        assert callable(self.task.api.check_method)
        assert callable(self.task.api.sample_dataset)
        assert callable(self.task.api.sample_method)
        assert hasattr(self.task.api.sample_dataset, "metadata")

    def test_api_is_consistent(self):
        """Test that a task's API is self-consistent"""
        adata = self.task.api.sample_dataset()
        assert self.task.api.check_dataset(adata)
        adata = self.task.api.sample_method(adata)
        assert self.task.api.check_method(adata)
