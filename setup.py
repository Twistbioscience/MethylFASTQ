import subprocess
import sys

from setuptools import setup, find_packages
from setuptools.command.test import test


setup_requirements = [
    'pytest-runner==5.2',
]


tests_requirements = [
    'pytest==5.4.3',
    'coverage==4.5.4',
    'pytest-cov==2.7.1',
]


class PyTest(test):
    """Define what happens when we run python setup.py test"""
    user_options = [('pytest-args=', 'a', "Arguments to pass to pytest")]

    def initialize_options(self):
        test.initialize_options(self)
        self.pytest_args = 'tests/'

    def run_tests(self):
        # Import here, because outside the eggs aren't loaded
        import shlex
        import pytest
        err_number = pytest.main(shlex.split(self.pytest_args))
        sys.exit(err_number)


setup(
    name='bfx-ngs-te-simulate-methyl',
    version='0.0.1',
    description="Simulate methylation target enrichment sequencing",
    author="Twist Bioscience",
    author_email='bioinformatics@twistbioscience.com',
    url='https://github.com/Twistbioscience/MethylFASTQ/',
    packages=find_packages(exclude=('tests*', 'docs')),
    include_package_data=True,
    install_requires=open('requirements.txt').readlines(),
    setup_requires=setup_requirements,
    tests_require=tests_requirements,
    zip_safe=False,
    entry_points={
        'console_scripts': ['ngs-te-simulate-methyl=bfx_ngs_te_simulate_methyl.methylFASTQ:main'],
    },
    classifiers=[
        'Programming Language :: Python :: 3.7',
    ],
    test_suite='tests',
    cmdclass={'tests': PyTest},
)