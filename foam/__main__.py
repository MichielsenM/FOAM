"""
Foam module interface.

Usage:
    foam version 
    foam test

"""

if __name__ == "__main__":

    import sys
    import foam

    if sys.argv[-1] == 'version':
        sys.exit(foam.__version__)
    
    if sys.argv[-1] == 'test':
        import pytest, os
        from pathlib import Path
        module_dir = str(Path(os.path.abspath(foam.__file__)).parent.parent)
        pytest.main([module_dir])

