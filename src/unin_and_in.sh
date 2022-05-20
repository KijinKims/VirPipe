pip uninstall -y virpipe
rm -rf virpipe.egg-info
rm -rf dist
python setup.py sdist
pip install dist/virpipe-1.0.tar.gz
