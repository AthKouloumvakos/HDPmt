'''
When HDPmt is installed in an enviroment,
with this script you can directry call the package as an excecutable.
$ python3 -m pip install ./ --use-feature=in-tree-build
$ hdpmt streamlit
'''
import os

import click
import streamlit.cli


@click.group()
def main():
    pass


@main.command('streamlit')
def main_streamlit():
    '''
    Run the main program in browser.
    '''
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, 'PyThea_app.py')
    args = []
    streamlit.cli._main_run(filename, args)


@main.command('help')
@click.pass_context
def help(ctx):
    """Print this help message."""
    # For 'pythea --help' instead of 'pythea help'.
    import sys
    sys.argv[1] = '--help'
    main()


@main.command('version')
def main_version():
    '''
    Print PyThea's version number.
    '''
    print('Showing version is not implemented yet.')


@main.command('docs')
def docs():
    '''
    Show documents in browser
    '''
    print('Showing document page in browser is not implemented yet.')

# TODO Include the test here.


if __name__ == '__main__':
    main()
