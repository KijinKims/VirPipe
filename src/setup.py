from setuptools import setup, find_packages

setup_requires = [
    ]

install_requires = [
    'docker',
    'pyfastx',
    ]

setup(
    name='VirPipe',
    version='1.0.9',
    description='VirPipe is a easy-to-use computational pipeline to identify virus sequences from high-throughput sequencing reads.',
    author='Kijin Kim',
    author_email='skkujin@gmail.com',
    url='https://github.com/KijinKims/VirPipe',
    packages=find_packages(),
    install_requires=install_requires,
    setup_requires=setup_requires,
    license='MIT',
    entry_points={
        'console_scripts': [
            'virpipe = virpipe.__main__:main',
            ],
        },
    )
