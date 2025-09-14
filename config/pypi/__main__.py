if __name__ == "__main__":
    import sys
    if "--prefix" in sys.argv:
        from . import get_tmm_dir
        print(get_tmm_dir())
        del get_tmm_dir
    del sys
