def get_tmm_dir():
    import os
    return os.path.dirname(__file__)


def get_config():
    conf = {}
    conf['TMM_DIR'] = get_tmm_dir()
    return conf
