def get_members(module):
    namespace = [attr for attr in dir(module) if not attr.startswith("_")]
    return [getattr(module, attr) for attr in namespace]

