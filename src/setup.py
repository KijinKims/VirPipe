import setuptools

setuptools.setup(
    name='virpipe',
    version='1.0',
    scripts=['./scripts/virpipe'],
    author='Kijin Kim',
    description='virpipe',
    packages=['virpipe'],
    install_requires=[
        'setuptools'
    ],
    include_package_data=True
)
