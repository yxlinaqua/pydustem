from setuptools import setup, find_packages

version = '0.1'

setup(
    name='pydustem',
    version=version,
    description="DustEM wrapper for Python",
    long_description="""\
    DustEM wrapper for Python""",
    classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords='',
    author='Tomohiko NAKAMURA',
    author_email='tnakamura@astron.s.u-tokyo.ac.jp',
    url='',
    license='',
    packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
    include_package_data=True,
    zip_safe=False,
    install_requires = [
    ],
    entry_points = {
        'console_scripts' : [
        ],
    },
    test_suite='nose.collector',
    tests_require=['Nose'],
)
