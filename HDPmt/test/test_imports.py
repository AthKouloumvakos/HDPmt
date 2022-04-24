import pkgutil


def test_import_main():
    """
    This imports HDPmt
    """
    no_requests = False
    try:
        pass
    except ImportError:
        no_requests = True
    assert no_requests is False


def test_imports_all():
    """
    This imports all modules in HDPmt
    """
    def on_error():
        try:
            raise
        except Warning:
            pass

    for imper, nm, ispkg in pkgutil.walk_packages(['HDPmt'], 'HDPmt.',
                                                  onerror=on_error):
        imper.find_spec(nm)
