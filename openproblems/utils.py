def get_members(module):
    namespace = [attr for attr in dir(module) if not attr.startswith("_")]
    return [getattr(module, attr) for attr in namespace]


def get_callable_members(module):
    return [member for member in get_members(module) if callable(member)]
