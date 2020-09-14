import gc


def garbage_collection():
    # runs gc.collect multiple times to free memory to OS
    # rather than just to Python
    gc.collect()
    gc.collect()
    gc.collect()
